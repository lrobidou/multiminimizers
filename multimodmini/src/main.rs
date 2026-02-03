use clap::Parser;
use serde::{Deserialize, Serialize};

#[derive(Debug, Parser)]
#[command(name = "ocmm")]
#[command(disable_help_subcommand = true)]
struct Args {
    #[arg(long)]
    size_read: usize,

    #[arg(long)]
    minimizer_size: usize,

    #[arg(long = "window", short)]
    w: usize,

    #[arg(long = "anchor-s")]
    anchor_s: Option<usize>,

    #[arg(short, default_value = "4")]
    r: usize,

    #[arg(long = "hashes", alias = "hash-count", default_value = "4")]
    hash_count: usize,

    #[arg(long = "seed")]
    base_seed: Option<u64>,

    #[arg(long = "verify")]
    verify: bool,
}

#[derive(Debug)]
pub struct Config {
    size_read: usize,
    k: usize,
    w: usize,
    t: usize,
    anchor_s: usize,
    hash_count: usize,
    base_seed: u64,
    verify: bool,
}

impl Config {
    pub fn new(
        size_read: usize,
        minimizer_size: usize,
        w: usize,
        r: usize,
        anchor_s: Option<usize>,
        hash_count: usize,
        base_seed: u64,
        verify: bool,
    ) -> Result<Self, String> {
        let k = minimizer_size;
        if k == 0 || w == 0 {
            return Err("Both --kmer and --window must be positive integers".into());
        }

        if r == 0 {
            return Err("--r must be at least 1".into());
        }
        if r > k {
            return Err("--r must not exceed --kmer".into());
        }

        if hash_count == 0 {
            return Err("--hashes must be at least 1".into());
        }

        let t = if k <= w {
            k
        } else {
            let rem = (k - r) % w;
            r + rem
        };
        if t == 0 {
            return Err("Derived anchor length t must be positive".into());
        }
        if t > k {
            return Err("Derived anchor length t can not exceed k".into());
        }

        let anchor_s = anchor_s.unwrap_or_else(|| std::cmp::max(1, t / 2));
        if anchor_s == 0 || anchor_s > t {
            return Err("--anchor-s must satisfy 1 <= anchor-s <= derived anchor length".into());
        }

        Ok(Self {
            size_read,
            k,
            w,
            t,
            anchor_s,
            hash_count,
            base_seed,
            verify,
        })
    }

    fn from_args(args: Args) -> Result<Self, String> {
        Self::new(
            args.size_read,
            args.minimizer_size,
            args.w,
            args.r,
            args.anchor_s,
            args.hash_count,
            args.base_seed.unwrap_or(0x9E37_79B1_85EB_CA87),
            args.verify,
        )
    }
}

fn main() {
    let args = Args::parse();
    let config = Config::from_args(args).expect("invalid values passed from the command line");
    let result = run(config).expect("error during execution");
    println!("{:?}", result);
}

fn random_dna_seq(len: usize) -> String {
    use rand::prelude::IndexedRandom;

    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::rng();
    let seq: Vec<u8> = (0..len).map(|_| *bases.choose(&mut rng).unwrap()).collect();
    String::from_utf8(seq).unwrap()
}

pub fn run(config: Config) -> Result<(f64, f64), String> {
    let seeds = build_seeds(config.base_seed, config.hash_count);
    let seq = random_dna_seq(config.size_read);
    process_sequence(&seq, &config, &seeds)
}

fn process_sequence(sequence: &str, config: &Config, seeds: &[Seed]) -> Result<(f64, f64), String> {
    let seq_len = sequence.len();
    let window_span = config.w + config.k - 1;
    if seq_len < window_span {
        eprintln!();
        return Err(format!(
            "Sequence '{}' shorter than required span ({} < {}), skipping.",
            sequence, seq_len, window_span
        ));
    }
    let num_kmers = seq_len - config.k + 1;
    let num_windows = seq_len - window_span + 1;
    let anchor_span = config.w + config.k - config.t;

    let mut segments_per_hash: Vec<Vec<Segment>> = Vec::with_capacity(seeds.len());
    for seed in seeds.iter() {
        let info = precompute_t_info(sequence.as_bytes(), config.t, config.anchor_s, seed);
        let segments = build_segments_for_hash(info, config, num_windows, anchor_span);
        segments_per_hash.push(segments);
    }

    let selections = select_segments(&segments_per_hash, num_kmers)?;

    if config.verify {
        verify_coverage(sequence, num_kmers, &selections)?;
    }

    let nb_superkmer = selections.len();
    let density = nb_superkmer as f64 / (config.size_read - config.k + 1) as f64;
    let bit_per_kmer_in_superkmer =
        get_bit_per_kmer_in_superkmer(seq_len, config.k, config.w, &selections).unwrap();
    Ok((density, bit_per_kmer_in_superkmer))
}

#[derive(Clone)]
struct Seed {
    order_seed: u64,
    s_seed: u64,
}

fn build_seeds(base_seed: u64, count: usize) -> Vec<Seed> {
    (0..count)
        .map(|idx| {
            let mut seed =
                base_seed.wrapping_add(0x9E37_79B1_85EB_CA87_u64.wrapping_mul(idx as u64 + 1));
            seed = mix64(seed);
            let order_seed = mix64(seed ^ 0xD1B5_4A32_D192_ED03);
            let s_seed = mix64(order_seed ^ 0x94D0_49BB_1331_11EB);
            Seed { order_seed, s_seed }
        })
        .collect()
}

#[derive(Clone)]
struct TMerInfo {
    order: u64,
    is_open: bool,
    is_closed: bool,
}

#[derive(Clone)]
struct Segment {
    minimizer_pos: usize,
    start: usize,
    end: usize,
}

#[derive(Clone)]
struct SegmentState {
    minimizer_pos: usize,
    start_window: usize,
    last_window: usize,
}

#[derive(Clone, Serialize, Deserialize)]
struct Selection {
    hash_idx: usize,
    minimizer_pos: usize,
    cover_start: usize,
    cover_end: usize,
}

fn precompute_t_info(sequence: &[u8], t: usize, s: usize, seed: &Seed) -> Vec<TMerInfo> {
    if sequence.len() < t || s == 0 || s > t {
        return Vec::new();
    }
    let num_tmers = sequence.len() - t + 1;
    let num_smers = sequence.len() - s + 1;
    let mut s_hashes = Vec::with_capacity(num_smers);
    for start in 0..num_smers {
        let slice = &sequence[start..start + s];
        s_hashes.push(hash64(seed.s_seed, slice));
    }
    let open_offset = (t - s) / 2;
    let max_offset = t - s;
    let mut infos = Vec::with_capacity(num_tmers);
    for pos in 0..num_tmers {
        let t_slice = &sequence[pos..pos + t];
        let order = hash64(seed.order_seed, t_slice);
        let mut min_hash = u64::MAX;
        let mut min_offset = 0usize;
        for rel in 0..=max_offset {
            let hash = s_hashes[pos + rel];
            if hash < min_hash || (hash == min_hash && rel < min_offset) {
                min_hash = hash;
                min_offset = rel;
            }
        }
        infos.push(TMerInfo {
            order,
            is_open: min_offset == open_offset,
            is_closed: min_offset == 0 || min_offset == max_offset,
        });
    }
    infos
}

fn build_segments_for_hash(
    infos: Vec<TMerInfo>,
    config: &Config,
    num_windows: usize,
    anchor_span: usize,
) -> Vec<Segment> {
    if num_windows == 0 || infos.is_empty() {
        return Vec::new();
    }
    let mut segments = Vec::new();
    let mut current: Option<SegmentState> = None;
    for window_start in 0..num_windows {
        let anchor_start = choose_anchor(&infos, window_start, anchor_span);
        let anchor_offset = anchor_start - window_start;
        let k_idx = anchor_offset % config.w;
        let k_pos = window_start + k_idx;

        match current.take() {
            Some(mut state) => {
                if state.minimizer_pos == k_pos {
                    state.last_window = window_start;
                    current = Some(state);
                } else {
                    segments.push(Segment {
                        minimizer_pos: state.minimizer_pos,
                        start: state.start_window,
                        end: state.last_window + config.w,
                    });
                    current = Some(SegmentState {
                        minimizer_pos: k_pos,
                        start_window: window_start,
                        last_window: window_start,
                    });
                }
            }
            None => {
                current = Some(SegmentState {
                    minimizer_pos: k_pos,
                    start_window: window_start,
                    last_window: window_start,
                });
            }
        }
    }
    if let Some(state) = current {
        segments.push(Segment {
            minimizer_pos: state.minimizer_pos,
            start: state.start_window,
            end: state.last_window + config.w,
        });
    }
    segments
}

fn select_segments(
    segments_per_hash: &[Vec<Segment>],
    num_kmers: usize,
) -> Result<Vec<Selection>, String> {
    let hash_count = segments_per_hash.len();
    if hash_count == 0 {
        return Ok(Vec::new());
    }
    let mut pointers = vec![0usize; hash_count];
    let mut coverage_idx = 0usize;
    let mut selections = Vec::new();

    while coverage_idx < num_kmers {
        let mut best: Option<(usize, &Segment)> = None;
        for hash_idx in 0..hash_count {
            let segments = &segments_per_hash[hash_idx];
            while let Some(seg) = segments.get(pointers[hash_idx]) {
                if coverage_idx >= seg.end {
                    pointers[hash_idx] += 1;
                } else {
                    break;
                }
            }
            let seg = match segments.get(pointers[hash_idx]) {
                Some(seg) if seg.start <= coverage_idx => seg,
                _ => continue,
            };
            let replace = match best {
                Some((_, best_seg)) => seg.end > best_seg.end,
                None => true,
            };
            if replace {
                best = Some((hash_idx, seg));
            }
        }
        let (hash_idx, seg) = best.ok_or_else(|| {
            format!(
                "No minimizer segment covers k-mer index {} (hashes evaluated: {})",
                coverage_idx, hash_count
            )
        })?;
        let cover_end = seg.end.min(num_kmers);
        selections.push(Selection {
            hash_idx,
            minimizer_pos: seg.minimizer_pos,
            cover_start: coverage_idx,
            cover_end,
        });
        coverage_idx = cover_end;
    }

    Ok(selections)
}

fn verify_coverage(
    sequence_name: &str,
    num_kmers: usize,
    selections: &[Selection],
) -> Result<(), String> {
    if num_kmers == 0 {
        return Ok(());
    }
    let mut covered = vec![false; num_kmers];
    for sel in selections {
        if sel.cover_start >= sel.cover_end || sel.cover_end > num_kmers {
            return Err(format!(
                "Invalid coverage range [{}..{}) for hash {} on sequence {}",
                sel.cover_start, sel.cover_end, sel.hash_idx, sequence_name
            ));
        }
        for idx in sel.cover_start..sel.cover_end {
            covered[idx] = true;
        }
    }
    if let Some(miss) = covered.iter().position(|covered| !covered) {
        return Err(format!(
            "Verification failed: k-mer {} in sequence '{}' is uncovered",
            miss, sequence_name
        ));
    }
    Ok(())
}

fn get_bit_per_kmer_in_superkmer(
    seq_len: usize,
    k: usize,
    w: usize,
    selections: &[Selection],
) -> Result<f64, String> {
    // Warning: I'm using a different notation here, because I'm more used to it
    let m = k;
    let k = w + m - 1;
    let mut bit_size = 0;
    let nb_kmers = seq_len - k + 1;
    for sel in selections {
        if sel.cover_start >= sel.cover_end || sel.cover_end > seq_len - m + 1 {
            return Err(format!(
                "Invalid coverage range [{}..{}) for hash {} on sequence of {} kmers",
                sel.cover_start, sel.cover_end, sel.hash_idx, nb_kmers
            ));
        }
        let nb_mini_covered = sel.cover_end - sel.cover_start;
        let size_superkmer_in_bases = nb_mini_covered + k - 1;
        let size_superkmer_in_bits = 2 * size_superkmer_in_bases;
        bit_size += size_superkmer_in_bits;
    }
    let bit_per_kmer_in_superkmer = bit_size as f64 / seq_len as f64;
    Ok(bit_per_kmer_in_superkmer)
}

fn choose_anchor(infos: &[TMerInfo], start: usize, length: usize) -> usize {
    let range_end = start + length;
    let mut best_open: Option<(usize, u64)> = None;
    let mut best_closed: Option<(usize, u64)> = None;
    let mut best_any: Option<(usize, u64)> = None;
    for pos in start..range_end {
        let info = match infos.get(pos) {
            Some(info) => info,
            None => continue,
        };
        if info.is_open {
            update_best(&mut best_open, pos, info.order);
        } else if info.is_closed {
            update_best(&mut best_closed, pos, info.order);
        }
        update_best(&mut best_any, pos, info.order);
    }
    best_open
        .or(best_closed)
        .or(best_any)
        .map(|(pos, _)| pos)
        .expect("Window must contain at least one anchor candidate")
}

fn update_best(candidate: &mut Option<(usize, u64)>, pos: usize, order: u64) {
    match candidate {
        Some((best_pos, best_order)) => {
            if order < *best_order || (order == *best_order && pos < *best_pos) {
                *candidate = Some((pos, order));
            }
        }
        None => *candidate = Some((pos, order)),
    }
}

fn hash64(mut seed: u64, bytes: &[u8]) -> u64 {
    const MULT: u64 = 0x9E37_79B1_85EB_CA87;
    for &b in bytes {
        seed ^= (b as u64).wrapping_mul(0x1000_0000_01B3);
        seed = seed
            .rotate_left(27)
            .wrapping_mul(MULT)
            .wrapping_add(0x94D0_49BB_1331_11EB);
    }
    mix64(seed ^ (bytes.len() as u64))
}

fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hash64_changes() {
        let a = hash64(123, b"ACGT");
        let b = hash64(123, b"TGCA");
        assert_ne!(a, b);
    }

    #[test]
    fn test_precompute_and_choose_anchor() {
        let sequence = b"ACGTACGTACGT".to_vec();
        let seed = Seed {
            order_seed: 1,
            s_seed: 2,
        };
        let infos = precompute_t_info(&sequence, 4, 2, &seed);
        assert!(!infos.is_empty());
        let anchor = choose_anchor(&infos, 0, 3);
        assert!(anchor < infos.len());
    }
}
