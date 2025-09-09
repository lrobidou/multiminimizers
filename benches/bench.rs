use criterion::{Criterion, criterion_group, criterion_main};

use simd_minimizers::seeded::canonical_minimizer_and_superkmer_positions;
use std::hint::black_box;

use sticky_mini::compute_superkmers_linear_streaming;

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

fn bench_f(c: &mut Criterion) {
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
    // c.bench_function("f::<10>", |b| b.iter(|| f::<10>()));
    // c.bench_function("f::<256>", |b| b.iter(|| f::<256>()));
    // c.bench_function("f::<4096>", |b| b.iter(|| f::<4096>()));
}

criterion_group!(benches, bench_f);
criterion_main!(benches);
