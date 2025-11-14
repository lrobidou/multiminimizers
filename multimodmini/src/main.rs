use std::borrow::Cow;
use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use multiminimizers::compute_superkmers_linear_streaming;


const SIZE_READ:usize = 50_000;

fn main() -> Result<(), Box<dyn Error>> {
    let config = Config::from_env()?;
    run(config)?;
    Ok(())
}

#[derive(Debug)]
struct Config {
    // input: String,
    emit_path: Option<String>,
    k: usize,
    w: usize,
    t: usize,
    anchor_s: usize,
    hash_count: usize,
    base_seed: u64,
    verify: bool,
    compare_classic: bool,
    classic_canonical: bool,
}

impl Config {
    fn from_env() -> Result<Self, Box<dyn Error>> {
        let mut args: Vec<String> = env::args().collect();
        // Skip program path.
        if args.len() == 1 {
            return Err(Box::from(Config::usage()));
        }
        // Remove the binary path.
        args.remove(0);
        let mut input = None;
        let mut emit_path = None;
        let mut k = None;
        let mut w = None;
        let mut anchor_s = None;
        let mut r = 4usize;
        let mut hash_count = 4usize;
        let mut base_seed: u64 = 0x9E37_79B1_85EB_CA87;
        let mut verify = false;
        let mut compare_classic = false;
        let mut classic_canonical = false;

        let mut iter = args.into_iter();
        while let Some(arg) = iter.next() {
            match arg.as_str() {
                "-i" | "--input" => {
                    input = Some(
                        iter.next()
                            .ok_or_else(|| format!("Missing value for {arg}"))?,
                    )
                }
                "-k" | "--kmer" => {
                    k = Some(
                        iter.next()
                            .ok_or_else(|| format!("Missing value for {arg}"))?
                            .parse::<usize>()
                            .map_err(|_| format!("Invalid integer for {arg}"))?,
                    );
                }
                "-w" | "--window" => {
                    w = Some(
                        iter.next()
                            .ok_or_else(|| format!("Missing value for {arg}"))?
                            .parse::<usize>()
                            .map_err(|_| format!("Invalid integer for {arg}"))?,
                    );
                }
                "--anchor-s" => {
                    anchor_s = Some(
                        iter.next()
                            .ok_or_else(|| "Missing value for --anchor-s".to_string())?
                            .parse::<usize>()
                            .map_err(|_| "Invalid integer for --anchor-s")?,
                    );
                }
                "--r" => {
                    r = iter
                        .next()
                        .ok_or_else(|| "Missing value for --r".to_string())?
                        .parse::<usize>()
                        .map_err(|_| "Invalid integer for --r")?;
                }
                "--hashes" | "--hash-count" => {
                    hash_count = iter
                        .next()
                        .ok_or_else(|| "Missing value for --hashes".to_string())?
                        .parse::<usize>()
                        .map_err(|_| "Invalid integer for --hashes")?;
                }
                "--seed" => {
                    base_seed = iter
                        .next()
                        .ok_or_else(|| "Missing value for --seed".to_string())?
                        .parse::<u64>()
                        .map_err(|_| "Invalid integer for --seed")?;
                }
                "--emit" => {
                    emit_path = Some(
                        iter.next()
                            .ok_or_else(|| "Missing value for --emit".to_string())?,
                    );
                }
                "--verify" => {
                    verify = true;
                }
                "--compare-classic" => {
                    compare_classic = true;
                }
                "--classic-canonical" => {
                    classic_canonical = true;
                }
                "--help" | "-h" => {
                    return Err(Box::from(Config::usage()));
                }
                other => {
                    return Err(Box::from(format!(
                        "Unknown argument '{other}'\n{}",
                        Config::usage()
                    )));
                }
            }
        }

        // let input = input.ok_or_else(|| Config::usage())?;
        let k = k.ok_or_else(|| "Missing required argument --kmer".to_string())?;
        let w = w.ok_or_else(|| "Missing required argument --window".to_string())?;
        if k == 0 || w == 0 {
            return Err(Box::from(
                "Both --kmer and --window must be positive integers",
            ));
        }
        if hash_count == 0 {
            return Err(Box::from("--hashes must be at least 1"));
        }
        if r == 0 {
            return Err(Box::from("--r must be at least 1"));
        }
        if r > k {
            return Err(Box::from("--r must not exceed --kmer"));
        }
        let t = if k <= w {
            k
        } else {
            let rem = (k - r) % w;
            r + rem
        };
        if t == 0 {
            return Err(Box::from("Derived anchor length t must be positive"));
        }
        if t > k {
            return Err(Box::from("Derived anchor length t can not exceed k"));
        }
        let anchor_s = anchor_s.unwrap_or_else(|| std::cmp::max(1, t / 2));
        if anchor_s == 0 || anchor_s > t {
            return Err(Box::from(
                "--anchor-s must satisfy 1 <= anchor-s <= derived anchor length",
            ));
        }

        Ok(Self {
            // input,
            emit_path,
            k,
            w,
            t,
            anchor_s,
            hash_count,
            base_seed,
            verify,
            compare_classic,
            classic_canonical,
        })
    }

    fn usage() -> String {
        r#"Usage: ocmm --input <fasta> --kmer <k> --window <w> [options]

Options:
  --anchor-s <s>        Length of the inner syncmer substring (default: max(1, t/2))
  --r <value>           Lower bound for anchor length t in the mod-minimizer (default: 4)
  --hashes <n>          Number of independent hash functions to combine (default: 4)
  --seed <u64>          Base seed used to derive per-hash permutations (default: 0x9E3779B185EBCA87)
  --emit <path|->       Emit tab-separated minimizers to the given file or stdout (-)
  --verify              After sampling, verify every k-mer is covered by the selected minimizers
  --compare-classic     Also run the classical multiminimizer with the same parameters (requires odd --window <= --kmer)
  --classic-canonical   Use canonical mode for the classical multiminimizer comparison
  -h, --help            Show this message
"#
        .to_string()
    }
}

#[derive(Debug)]
struct Sequence {
    name: String,
    data: Vec<u8>,
}

fn read_fasta<P: AsRef<Path>>(path: P) -> io::Result<Vec<Sequence>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut sequences = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();
    use std::mem::take;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if let Some(name) = current_name.take() {
                let data = take(&mut current_seq);
                sequences.push(Sequence { name, data });
            }
            let header = line[1..].trim().to_string();
            current_name = Some(if header.is_empty() {
                format!("seq{}", sequences.len() + 1)
            } else {
                header
            });
        } else {
            let trimmed = line.trim();
            if !trimmed.is_empty() {
                for byte in trimmed.as_bytes() {
                    let upper = byte.to_ascii_uppercase();
                    current_seq.push(upper);
                }
            }
        }
    }
    if let Some(name) = current_name {
        sequences.push(Sequence {
            name,
            data: current_seq,
        });
    }
    Ok(sequences)
}

fn random_dna_seq(len: usize) -> String {
    use rand::prelude::IndexedRandom;

    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::rng();
    let seq: Vec<u8> = (0..len).map(|_| *bases.choose(&mut rng).unwrap()).collect();
    String::from_utf8(seq).unwrap()
}

fn run(config: Config) -> Result<(), Box<dyn Error>> {
    // let sequences = read_fasta(&config.input)?;
    // if sequences.is_empty() {
    //     return Err(Box::from("Input FASTA does not contain any sequences"));
    // }
    let mut writer: Option<Box<dyn Write>> = match config.emit_path.as_deref() {
        Some("-") => Some(Box::new(BufWriter::new(io::stdout()))),
        Some(path) => Some(Box::new(BufWriter::new(File::create(path)?))),
        None => None,
    };
    if let Some(w) = writer.as_mut() {
        writeln!(
            w,
            "#sequence\tcover_start\tcover_end\tminimizer_pos\tkmer\thash_ix"
        )?;
    }

    let seeds = build_seeds(config.base_seed, config.hash_count);
    let mut summary = Summary::new(config.k);
    let seq = random_dna_seq(SIZE_READ);
    let writer_ref = writer.as_mut().map(|w| w.as_mut() as &mut dyn Write);
    process_sequence(&seq, &config, &seeds, writer_ref, &mut summary)?;
    if let Some(w) = writer.as_mut() {
        w.flush()?;
    }
    let classic_stats =
    
    // if config.compare_classic {
        // Some(
            // run_classic_comparison(&sequences, &config)
                // .map_err(|msg| Box::<dyn Error>::from(msg))?,
        // )
    // } else {
        None
    // }
    ;
    summary.print(classic_stats.as_ref());
    Ok(())
}

fn process_sequence(
    sequence: &str,
    config: &Config,
    seeds: &[Seed],
    mut writer: Option<&mut dyn Write>,
    summary: &mut Summary,
) -> io::Result<()> {
    let seq_len = sequence.len();
    let window_span = config.w + config.k - 1;
    if seq_len < window_span {
        eprintln!(
            "Sequence '{}' shorter than required span ({} < {}), skipping.",
            sequence, seq_len, window_span
        );
        return Ok(());
    }
    let num_kmers = seq_len - config.k + 1;
    let num_windows = seq_len - window_span + 1;
    let anchor_span = config.w + config.k - config.t;

    let mut segments_per_hash: Vec<Vec<Segment>> = Vec::with_capacity(seeds.len());
    for seed in seeds.iter() {
        let info = precompute_t_info(&sequence.as_bytes(), config.t, config.anchor_s, seed);
        let segments = build_segments_for_hash(info, config, num_windows, anchor_span);
        segments_per_hash.push(segments);
    }

    let selections = select_segments(&segments_per_hash, num_kmers).map_err(|msg| {
        io::Error::new(io::ErrorKind::Other, format!("{} ({})", msg, sequence))
    })?;

    if config.verify {
        verify_coverage(&sequence, num_kmers, &selections)
            .map_err(|msg| io::Error::new(io::ErrorKind::Other, msg))?;
    }

    summary.total_kmers += num_kmers;
    summary.total_windows += num_windows;
    summary.sequences += 1;
    summary.selected_minimizers += selections.len();

    if let Some(w) = writer.as_mut() {
        for sel in selections {
            let kmer = &sequence.as_bytes()[sel.minimizer_pos..sel.minimizer_pos + config.k];
            let kmer_str = String::from_utf8_lossy(kmer);
            writeln!(
                *w,
                "{}\t{}\t{}\t{}\t{}\t{}",
                sequence,
                sel.cover_start,
                sel.cover_end,
                sel.minimizer_pos,
                kmer_str,
                sel.hash_idx
            )?;
        }
    }
    Ok(())
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

#[derive(Clone)]
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

const MAX_CLASSIC_HASHES: usize = 16;

struct ClassicStats {
    minimizers: usize,
}

fn run_classic_comparison(sequences: &[Sequence], config: &Config) -> Result<ClassicStats, String> {
    let minimizer_len = derive_classic_minimizer_len(config)?;
    let window_span = config.w + config.k - 1;
    let minimizers = if config.classic_canonical {
        dispatch_classic_comparison::<true>(sequences, config, minimizer_len, window_span)?
    } else {
        dispatch_classic_comparison::<false>(sequences, config, minimizer_len, window_span)?
    };
    Ok(ClassicStats { minimizers })
}

fn derive_classic_minimizer_len(config: &Config) -> Result<usize, String> {
    if config.w % 2 == 0 {
        return Err(
            "Classical comparison requires an odd --window (sticky minimizers need odd width)"
                .to_string(),
        );
    }
    let k_minus_w = config
        .k
        .checked_sub(config.w)
        .ok_or_else(|| "Classical comparison requires --window <= --kmer".to_string())?;
    Ok(k_minus_w + 1)
}

fn dispatch_classic_comparison<const CANONICAL: bool>(
    sequences: &[Sequence],
    config: &Config,
    minimizer_len: usize,
    window_span: usize,
) -> Result<usize, String> {
    let k = config.k;
    match config.hash_count {
        1 => Ok(classic_total::<CANONICAL, 1>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        2 => Ok(classic_total::<CANONICAL, 2>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        3 => Ok(classic_total::<CANONICAL, 3>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        4 => Ok(classic_total::<CANONICAL, 4>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        5 => Ok(classic_total::<CANONICAL, 5>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        6 => Ok(classic_total::<CANONICAL, 6>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        7 => Ok(classic_total::<CANONICAL, 7>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        8 => Ok(classic_total::<CANONICAL, 8>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        9 => Ok(classic_total::<CANONICAL, 9>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        10 => Ok(classic_total::<CANONICAL, 10>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        11 => Ok(classic_total::<CANONICAL, 11>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        12 => Ok(classic_total::<CANONICAL, 12>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        13 => Ok(classic_total::<CANONICAL, 13>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        14 => Ok(classic_total::<CANONICAL, 14>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        15 => Ok(classic_total::<CANONICAL, 15>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        16 => Ok(classic_total::<CANONICAL, 16>(
            sequences,
            k,
            minimizer_len,
            window_span,
        )),
        other => Err(format!(
            "Classical comparison currently supports up to {MAX_CLASSIC_HASHES} hash functions (requested {other})"
        )),
    }
}

fn classic_total<const CANONICAL: bool, const N: usize>(
    sequences: &[Sequence],
    k: usize,
    minimizer_len: usize,
    window_span: usize,
) -> usize {
    sequences
        .iter()
        .filter(|seq| seq.data.len() >= window_span)
        .map(|seq| classic_sequence_superkmers::<CANONICAL, N>(seq, k, minimizer_len))
        .sum()
}

fn classic_sequence_superkmers<const CANONICAL: bool, const N: usize>(
    sequence: &Sequence,
    k: usize,
    minimizer_len: usize,
) -> usize {
    let sanitized = sanitize_sequence_for_classic(&sequence.data);
    compute_superkmers_linear_streaming::<N, CANONICAL>(sanitized.as_ref(), k, minimizer_len)
        .map(|iter| iter.count())
        .unwrap_or(0)
}

fn sanitize_sequence_for_classic(data: &[u8]) -> Cow<'_, [u8]> {
    if data.iter().any(|base| *base == b'N') {
        let mut owned = data.to_vec();
        for base in &mut owned {
            if *base == b'N' {
                *base = b'A';
            }
        }
        Cow::Owned(owned)
    } else {
        Cow::Borrowed(data)
    }
}

struct Summary {
    k:usize,
    total_kmers: usize,
    total_windows: usize,
    selected_minimizers: usize,
    sequences: usize,
}

impl Summary {
    fn new(k:usize) -> Self {
        Self {
            k,
            total_kmers: 0,
            total_windows: 0,
            selected_minimizers: 0,
            sequences: 0,
        }
    }

    fn print(&self, comparison: Option<&ClassicStats>) {
        if self.total_kmers == 0 {
            eprintln!("No windows long enough for the requested parameters.");
            return;
        }
        if let Some(classic) = comparison {
            println!("ocmm\t{}", self.selected_minimizers);
            println!("classic\t{}", classic.minimizers);
            let delta = self.selected_minimizers as i128 - classic.minimizers as i128;
            println!("delta\t{}", delta);
            if classic.minimizers > 0 {
                println!(
                    "ratio\t{:.6}",
                    self.selected_minimizers as f64 / classic.minimizers as f64
                );
            }
        } else {
            // print density
            println!("{}", self.selected_minimizers as f64 / (SIZE_READ-self.k+1) as f64);
        }
    }
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
