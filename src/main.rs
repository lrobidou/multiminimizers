use serde::{Deserialize, Serialize};
use std::fs;
use sticky_mini::{
    superkmer::{NoAnchor, Superkmer},
    superkmers_computation::compute_superkmers_linear_streaming,
};

mod simd_canonical_sticky_minimizer_iterator;
mod subsequence;
mod superkmer;
mod superkmers_computation;
mod two_bits;

type Minimizer = u64;

use itertools::Itertools;

#[derive(Serialize, Deserialize)]
struct DataToPlot {
    size_read: usize,
    k: usize,
    m: usize,
    limit_overhead: i32,
    overhead: Vec<f64>,
    limit_avg_superkmer_size: usize,
    avg_superkmer_size: Vec<f64>,
    limit_nb_superkmer: f64,
    nb_superkmer: Vec<usize>,
}

impl DataToPlot {
    pub fn new(
        size_read: usize,
        k: usize,
        m: usize,
        overhead: Vec<f64>,
        avg_superkmer_size: Vec<f64>,
        nb_superkmer: Vec<usize>,
    ) -> Self {
        let limit_overhead = 2;
        let limit_avg_superkmer_size = 2 * k - m;
        let limit_nb_superkmer = size_read as f64 / k as f64;
        Self {
            size_read,
            k,
            m,
            limit_overhead,
            overhead,
            limit_avg_superkmer_size,
            avg_superkmer_size,
            limit_nb_superkmer,
            nb_superkmer,
        }
    }

    pub fn from_array(size_read: usize, k: usize, m: usize, arr: [(f64, f64, usize); 16]) -> Self {
        let mut overheads = vec![];
        let mut avgs_superkmer_size = vec![];
        let mut nbs_superkmer = vec![];

        for (overhead, avg_superkmer_size, nb_superkmer) in arr {
            overheads.push(overhead);
            avgs_superkmer_size.push(avg_superkmer_size);
            nbs_superkmer.push(nb_superkmer);
        }
        Self::new(
            size_read,
            k,
            m,
            overheads,
            avgs_superkmer_size,
            nbs_superkmer,
        )
    }
}

fn random_dna_seq(len: usize) -> String {
    use rand::prelude::IndexedRandom;

    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::rng();
    let seq: Vec<u8> = (0..len).map(|_| *bases.choose(&mut rng).unwrap()).collect();
    String::from_utf8(seq).unwrap()
}

fn nb_base_in_representation(v: &[Superkmer<'_, NoAnchor>]) -> usize {
    v.iter().map(|sk| sk.superkmer.len()).sum()
}

fn canonical_get_overhead_and_avg_size_and_nb_of_sk<const N: usize>(
    read: &str,
    k: usize,
    m: usize,
) -> (f64, f64, usize) {
    let superkmer_iter = compute_superkmers_linear_streaming::<N, true>(read.as_bytes(), k, m);
    let superkmers: Vec<Superkmer<'_, NoAnchor>> = superkmer_iter.unwrap().collect_vec();

    let nb_base = nb_base_in_representation(&superkmers);
    let overhead = nb_base as f64 / read.len() as f64;
    let avg_size = nb_base as f64 / superkmers.len() as f64;
    let nb_sk = superkmers.len();

    (overhead, avg_size, nb_sk)
}

fn non_canonical_get_overhead_and_avg_size_and_nb_of_sk<const N: usize>(
    read: &str,
    k: usize,
    m: usize,
) -> (f64, f64, usize) {
    let superkmer_iter = compute_superkmers_linear_streaming::<N, false>(read.as_bytes(), k, m);
    let superkmers: Vec<Superkmer<'_, NoAnchor>> = superkmer_iter.unwrap().collect_vec();

    let nb_base = nb_base_in_representation(&superkmers);
    let overhead = nb_base as f64 / read.len() as f64;
    let avg_size = nb_base as f64 / superkmers.len() as f64;
    let nb_sk = superkmers.len();

    (overhead, avg_size, nb_sk)
}

fn main() {
    let size_read = 10_000_000;
    let k = 301;
    let m = 21;
    let read = random_dna_seq(size_read);

    let arr: [(f64, f64, usize); 16] = collect_macro::collect!((
        16,
        canonical_get_overhead_and_avg_size_and_nb_of_sk,
        &read,
        k,
        m
    ));
    let data = DataToPlot::from_array(size_read, k, m, arr);
    let data_str = serde_json::to_string(&data).unwrap();
    fs::write("data_canonical.txt", data_str).expect("Should be able to write to `data/txt`");

    let arr = collect_macro::collect!((
        16,
        non_canonical_get_overhead_and_avg_size_and_nb_of_sk,
        &read,
        k,
        m
    ));
    let data = DataToPlot::from_array(size_read, k, m, arr);
    let data_str = serde_json::to_string(&data).unwrap();
    fs::write("data_non_canonical.txt", data_str).expect("Should be able to write to `data/txt`");

    println!("output written to data_canonical.txt and data_non_canonical.txt");
}

// #[cfg(test)]
// mod tests {
//     use super::*;

// fn calculate_hash<T: std::hash::Hash>(t: &T, seed: u64) -> u64 {
//     let mut s = DefaultHashBuilder::with_seed(seed);
//     t.hash(&mut s);
//     s.finish()
// }

// Wrapper around a Hasher that prints when hashing.
// struct DebugHasher<H> {
//     inner: H,
// }

// impl<H: Hasher> Hasher for DebugHasher<H> {
//     fn finish(&self) -> u64 {
//         let f = self.inner.finish();
//         println!("finish = {f}");
//         f
//     }

//     fn write(&mut self, bytes: &[u8]) {
//         // println!("Hashing bytes: {:?}", bytes);
//         self.inner.write(bytes);
//     }

//     // override other write_* methods similarly if needed
// }

// struct DebugBuildHasher<B> {
//     inner: B,
// }

// impl<B> DebugBuildHasher<B> {
//     fn new(inner: B) -> Self {
//         DebugBuildHasher { inner }
//     }
// }

// impl<B: BuildHasher> BuildHasher for DebugBuildHasher<B> {
//     type Hasher = DebugHasher<B::Hasher>;

//     fn build_hasher(&self) -> Self::Hasher {
//         DebugHasher {
//             inner: self.inner.build_hasher(),
//         }
//     }
// }
//     #[test]
//     fn test_seed_effect() {
//         // let read = random_dna_seq(31);
//         let read = String::from("GTAGAGGATCAGCCCCTGCACTCCCTTCCTTGTTTTACAGGCTCGAATTTGTT");
//         let read = read.as_bytes();
//         let k = 31;
//         let m = 11;

//         // seed 0
//         let hasher = DebugBuildHasher::new(DefaultHashBuilder::with_seed(0));

//         let hasher = ahash::RandomState::with_seeds(
//             0x0123_4567_89AB_ABAB,
//             0x89AB_CDEF_BB23_4567,
//             0xFEDC_BA98_76C4_3210,
//             0x7654_3210_FEDC_BA98,
//         );

//         // let minimizer_iter = MinimizerBuilder::<u64>::new()
//         //     .minimizer_size(m)
//         //     .width((k - m + 1) as u16)
//         //     .hasher(hasher)
//         //     // .seed(0)
//         //     .iter(read);
//         println!("=== {}", hasher.hash_one(1));

//         let minimizer_iter: CanonicalMinimizerPosIterator<'_, u64, ahash::RandomState> =
//             CanonicalMinimizerPosIterator::new(read, m, (k - m + 1) as u16, hasher, ENCODING);
//         let minimizers = minimizer_iter.collect_vec();

//         println!("");

//         // seed 1
//         // let hasher = DebugBuildHasher::new(DefaultHashBuilder::with_seed(1));
//         let hasher = ahash::RandomState::with_seeds(
//             0x0123_4567_89AB_CDEF,
//             0x89AB_CDEF_0123_4567,
//             0xFEDC_BA98_7654_3211,
//             0x7654_3210_FEDC_BA98,
//         );
//         println!("=== {}", hasher.hash_one(1));
//         let minimizer_iter_seed: CanonicalMinimizerPosIterator<'_, u64, ahash::RandomState> =
//             CanonicalMinimizerPosIterator::new(read, m, (k - m + 1) as u16, hasher, ENCODING);
//         let minimizers_seed = minimizer_iter_seed.collect_vec();

//         // same value
//         assert!(minimizers != minimizers_seed);
//     }
// }
