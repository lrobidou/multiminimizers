use crate::{
    superkmer::{NoAnchor, Superkmer},
    superkmers_computation::compute_superkmers_linear_streaming,
};
use std::fs;

use itertools::Itertools;

mod my_minimizer_iter;
mod subsequence;
mod superkmer;
mod superkmers_computation;
mod two_bits;

type Minimizer = u64;

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

fn get_overhead_and_avg_size<const N: usize>(read: &str, k: usize, m: usize) -> (f64, f64) {
    let superkmer_iter = compute_superkmers_linear_streaming::<N>(read.as_bytes(), k, m);
    let superkmers: Vec<Superkmer<'_, NoAnchor>> = superkmer_iter.unwrap().collect_vec();

    let nb_base = nb_base_in_representation(&superkmers);
    let overhead = nb_base as f64 / read.len() as f64;
    let avg_size = nb_base as f64 / superkmers.len() as f64;

    (overhead, avg_size)
}

fn main() {
    let size_read = 1000000;
    let k = 201;
    let m = 17;
    let read = random_dna_seq(size_read);

    let arr = collect_macro::collect!((16, get_overhead_and_avg_size, &read, k, m));

    let data = format!("{:?}", arr);
    fs::write("data.txt", data).expect("Should be able to write to `data/txt`");

    // println!("{:?}", arr);
    println!("output written to data.txt");
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
