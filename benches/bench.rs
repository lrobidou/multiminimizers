use criterion::{Criterion, criterion_group, criterion_main};

use multiminimizers::compute_superkmers_linear_streaming;
use multiminimizers::needless_modules::{Index, LinesIter};
use simd_minimizers::seeded::canonical_minimizer_and_superkmer_positions;
use std::hint::black_box;

use std::fs::{self, File};
use std::io::Write;
use uuid::Uuid;

// TODO stop duplicating this
fn random_dna_seq(len: usize) -> String {
    use rand::prelude::IndexedRandom;

    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::rng();
    let seq: Vec<u8> = (0..len).map(|_| *bases.choose(&mut rng).unwrap()).collect();
    String::from_utf8(seq).unwrap()
}

fn original_iterator(seq: &[u8], k: usize, m: usize) -> usize {
    let w = k - m + 1;
    let mut min_pos_vec = vec![];
    let mut sks_pos_vec = vec![];
    canonical_minimizer_and_superkmer_positions(seq, m, w, 0, &mut min_pos_vec, &mut sks_pos_vec);

    sks_pos_vec.len()
}

fn multi_mini_iterator<const N: usize>(seq: &[u8], k: usize, m: usize) -> usize {
    let superkmer_iter = compute_superkmers_linear_streaming::<N, true>(seq, k, m);
    superkmer_iter.unwrap().count()
}

fn bench_sk_iter(c: &mut Criterion) {
    let read_size = 100000;
    let k = 31;
    let m = 19;
    let read = random_dna_seq(read_size);

    let mut group = c.benchmark_group("Superkmer iterator");

    group.bench_function("original iterator", |b| {
        b.iter(|| original_iterator(black_box(read.as_bytes()), black_box(k), black_box(m)))
    });
    group.bench_function("multi_mini_iterator::<1>", |b| {
        b.iter(|| multi_mini_iterator::<1>(black_box(read.as_bytes()), black_box(k), black_box(m)))
    });
    group.bench_function("multi_mini_iterator::<2>", |b| {
        b.iter(|| multi_mini_iterator::<2>(black_box(read.as_bytes()), black_box(k), black_box(m)))
    });
    group.bench_function("multi_mini_iterator::<3>", |b| {
        b.iter(|| multi_mini_iterator::<3>(black_box(read.as_bytes()), black_box(k), black_box(m)))
    });
    group.bench_function("multi_mini_iterator::<4>", |b| {
        b.iter(|| multi_mini_iterator::<4>(black_box(read.as_bytes()), black_box(k), black_box(m)))
    });

    group.finish();
}

fn build_index<const NB_HASH: usize>(input: &str, k: usize, m: usize) {
    let sequences = LinesIter::from_path(input);

    let mut index: Index<NB_HASH> = Index::new(k, m);
    index.add_reads(sequences);
    let output = format!("{}.needless", input);

    index
        .dump(output)
        .expect("should have been able to dump the index");
}

fn bench_index_hash(c: &mut Criterion) {
    // write input
    let folder_name = Uuid::new_v4().to_string();
    fs::create_dir(&folder_name).unwrap();
    let read_size = 100000;
    let file_path = format!("{}/index.fa", folder_name);
    let mut file = File::create(&file_path).unwrap();
    let dna_sequence = random_dna_seq(read_size);
    writeln!(file, ">header").unwrap();
    writeln!(file, "{}", dna_sequence).unwrap();

    let k = 31;
    let m = 19;

    let mut group = c.benchmark_group("Index");

    group.bench_function("1 hash function", |b| {
        b.iter(|| build_index::<1>(black_box(&file_path), black_box(k), black_box(m)))
    });
    group.bench_function("2 hash function", |b| {
        b.iter(|| build_index::<2>(black_box(&file_path), black_box(k), black_box(m)))
    });
    group.bench_function("3 hash function", |b| {
        b.iter(|| build_index::<3>(black_box(&file_path), black_box(k), black_box(m)))
    });
    group.bench_function("4 hash function", |b| {
        b.iter(|| build_index::<4>(black_box(&file_path), black_box(k), black_box(m)))
    });
    group.bench_function("5 hash function", |b| {
        b.iter(|| build_index::<5>(black_box(&file_path), black_box(k), black_box(m)))
    });
    group.bench_function("6 hash function", |b| {
        b.iter(|| build_index::<6>(black_box(&file_path), black_box(k), black_box(m)))
    });
    group.bench_function("7 hash function", |b| {
        b.iter(|| build_index::<7>(black_box(&file_path), black_box(k), black_box(m)))
    });
    group.bench_function("8 hash function", |b| {
        b.iter(|| build_index::<8>(black_box(&file_path), black_box(k), black_box(m)))
    });

    group.finish();

    let dest = "bench_input";
    let _ = fs::remove_dir_all(dest);
    fs::rename(&folder_name, dest).unwrap();
}

criterion_group!(benches, bench_sk_iter, bench_index_hash);
criterion_main!(benches);
