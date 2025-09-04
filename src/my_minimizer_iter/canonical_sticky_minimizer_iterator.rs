// use super::CanonicalMinimizerIterator;
use minimizer_iter::iterator::CanonicalMinimizerIterator;
// use simd_minimizers::seed::canonical_minimizer_and_superkmer_positions;

use core::hash::Hash;
use num_traits::{AsPrimitive, PrimInt};
use std::fmt::Debug;
use std::iter::Peekable;

use xxhash_rust::xxh3::Xxh3Builder;

// pub struct Xxh3BuildHasher {
//     seed: u64,
// }

// impl Xxh3BuildHasher {
//     pub fn new(seed: u64) -> Self {
//         Self { seed }
//     }
// }

// impl BuildHasher for Xxh3BuildHasher {
//     type Hasher = Xxh3;

//     fn build_hasher(&self) -> Self::Hasher {
//         Xxh3::with_seed(self.seed)
//     }
// }

/// Sets the binary encoding of the bases.
pub const fn get_encoding(a: u8, c: u8, g: u8, t: u8) -> [u8; 256] {
    let mut encoding = [0; 256];
    encoding[b'A' as usize] = a;
    encoding[b'a' as usize] = a;
    encoding[b'C' as usize] = c;
    encoding[b'c' as usize] = c;
    encoding[b'G' as usize] = g;
    encoding[b'g' as usize] = g;
    encoding[b'T' as usize] = t;
    encoding[b't' as usize] = t;
    encoding
}

pub const ENCODING: [u8; 256] = get_encoding(0b00000000, 0b00000001, 0b00000011, 0b00000010);

// pub struct Xxh3BuildHasher {
//     seed: u64,
// }

// impl Xxh3BuildHasher {
//     pub fn new(seed: u64) -> Self {
//         Self { seed }
//     }
// }

// impl BuildHasher for Xxh3BuildHasher {
//     type Hasher = Xxh3;

//     fn build_hasher(&self) -> Self::Hasher {
//         Xxh3::with_seed(self.seed)
//     }
// }

struct MiniIterator<I, T: PrimInt + Hash = u64>
where
    I: Iterator<Item = (T, usize, bool)>,
{
    minimizer_positions: I,
    prev: Option<(T, usize, bool)>,
    prev_sk_start: usize,
    width: u16,
}

impl<I, T> MiniIterator<I, T>
where
    I: Iterator<Item = (T, usize, bool)>,
    T: PrimInt + Hash,
{
    fn new(mut iter: I, width: u16) -> Self {
        let prev = iter.next();
        Self {
            minimizer_positions: iter,
            prev,
            prev_sk_start: 0,
            width,
        }
    }
}

impl<I, T> Iterator for MiniIterator<I, T>
where
    I: Iterator<Item = (T, usize, bool)>,
    T: PrimInt + Hash + Debug, // TODO remove debug
{
    type Item = (T, usize, usize, bool); // (val, minimizer_start, superkmer_start, canonicity)

    fn next(&mut self) -> Option<Self::Item> {
        // println!("[MiniIterator] iter is {:?}", self.prev);

        let (val, prev_mini_start, c) = self.prev?;
        let prev_sk_start = self.prev_sk_start;

        self.prev = self.minimizer_positions.next();
        println!(
            "will return {:?}",
            Some((val, prev_mini_start, prev_sk_start, c))
        );

        if let Some((next_val, next_prev_mini_start, _next_c)) = self.prev {
            println!(
                "just got {:?} {} \n ",
                Some((next_val, next_prev_mini_start, _next_c)),
                self.width
            );
            self.prev_sk_start = if next_val < val {
                // new minimizer encountererd by sliding the window
                next_prev_mini_start - self.width as usize - 1 // TODO check
            } else {
                // previous mini is now out of the window
                prev_mini_start + 1
            };
        }

        // println!("[MiniIterator] next iter will be {:?}", self.prev);
        // =;
        Some((val, prev_mini_start, prev_sk_start, c))
    }
}

/// An iterator over `N` `CanonicalMinimizerIterator` of a sequence.
/// It requires an odd width to break ties between multiple minimizers.
pub struct CanonicalStickyMinimizerIterator<
    'a,
    const N: usize,
    T: PrimInt + Hash + Debug + 'static = u64, // TODO remove debug
> where
    u8: AsPrimitive<T>,
{
    minimizer_iters: [Peekable<MiniIterator<CanonicalMinimizerIterator<'a, T, Xxh3Builder>, T>>; N],
    minimizers: [Option<(T, usize, usize, bool)>; N],
    current_mini: Option<usize>,
}

/// Gets the further minimizer if it exists, else returns None.
fn get_max_pos<T, const N: usize>(arr: [Option<(T, usize, usize, bool)>; N]) -> Option<usize> {
    arr.iter()
        .enumerate()
        .filter_map(|(i, opt)| opt.as_ref().map(|(_, second, _, _)| (i, second)))
        .max_by_key(|&(_, second)| second)
        .map(|(i, _)| i)
    // arr.into_iter().flaltten().max_by_key(|x| x.1)
}

impl<'a, const N: usize, T: PrimInt + Debug + Hash + 'static>
    CanonicalStickyMinimizerIterator<'a, N, T>
where
    u8: AsPrimitive<T>,
{
    pub fn new(seq: &'a [u8], minimizer_size: usize, width: u16) -> Self {
        assert_eq!(
            width % 2,
            1,
            "width must be odd to break ties between multiple minimizers"
        );
        // OPTIMIZE make an unsafe constructor for `CanonicalMinimizerIterator` that doesn't check if the width is odd ?
        let mut minimizer_iters: [Peekable<
            MiniIterator<CanonicalMinimizerIterator<'_, T, Xxh3Builder>, T>,
        >; N] = std::array::from_fn(|i| {
            let build_hasher = Xxh3Builder::new().with_seed(i as u64);
            let inner: CanonicalMinimizerIterator<'_, T, Xxh3Builder> =
                CanonicalMinimizerIterator::new(seq, minimizer_size, width, build_hasher, ENCODING);
            let iter = MiniIterator::new(inner, width);
            iter.peekable()
        });
        let minimizers = std::array::from_fn(|i| minimizer_iters[i].next());
        // println!("minimizers at the start = {:?}", minimizers);
        let current_mini = get_max_pos(minimizers);

        // the current minimizer is self.minimizers[self.current_mini?]?
        // invariants:
        // 1. the current minimizer is the farthest in the read
        // 2. all minimizers'superkmer covers the first kmer of the superkmer denoted by the current minimizer
        // 3. the next minimizer of all minimizers does not cover said kmer

        Self {
            minimizer_iters,
            minimizers,
            current_mini,
        }
    }

    fn ensure_all_iters_cover_sk_starting_pos(&mut self, start_sk: usize) {
        // println!("start");
        for i in 0..self.minimizer_iters.len() {
            // println!("check {i}");
            while let Some((_, _, next_start_sk, _)) = self.minimizer_iters[i].peek()
                && *next_start_sk <= start_sk
            {
                // println!("increase");
                self.minimizers[i] = self.minimizer_iters[i].next()
            }
        }
        // println!("stop");
        // println!("minimizers = {:?}", self.minimizers);
    }
}

impl<'a, const N: usize, T: PrimInt + Hash + Debug + 'static> Iterator
    for CanonicalStickyMinimizerIterator<'a, N, T>
where
    u8: AsPrimitive<T>,
{
    type Item = (usize, usize, bool);

    fn next(&mut self) -> Option<Self::Item> {
        // early exits if nothing more to return
        let current_mini = self.current_mini?;
        let (_val, sk_start, mini_start, c) = self.minimizers[current_mini]?;
        // move the current minimizer iterator

        self.minimizers[current_mini] = self.minimizer_iters[current_mini].next();

        // invariants might be broken, let's fix them
        if let Some((_, _, start_sk, _)) = self.minimizers[current_mini] {
            // fix invariants 2 and 3
            self.ensure_all_iters_cover_sk_starting_pos(start_sk);
            // fix invarient 1
            self.current_mini = get_max_pos(self.minimizers);
        }

        Some((sk_start, mini_start, c))
    }
}

// #[cfg(test)]
// mod tests {
//     use itertools::Itertools;
//     use minimizer_iter::{DefaultHashBuilder, iterator::CanonicalMinimizerPosIterator};

//     use super::*;

//     #[test]
//     fn test_iter_single_kmer() {
//         let kmer = "TGATGAGTACGTAGCGAAAAAAAAAAGGGTACGTGCATGCAGTGACG".as_bytes();
//         let k = kmer.len();
//         let m = 11;

//         let build_hasher = Xxh3Builder::new().with_seed(0 as u64);
//         let minimizer_iter: CanonicalMinimizerIterator<u64, _> =
//             CanonicalMinimizerIterator::new(kmer, m, (k - m + 1) as u16, build_hasher, ENCODING);
//         let minimizer_iter = minimizer_iter.collect_vec();
//         assert_eq!(minimizer_iter.len(), 1);
//     }

//     fn random_dna_seq(len: usize) -> String {
//         use rand::prelude::IndexedRandom;

//         let bases = [b'A', b'C', b'G', b'T'];
//         let mut rng = rand::rng();
//         let seq: Vec<u8> = (0..len).map(|_| *bases.choose(&mut rng).unwrap()).collect();
//         String::from_utf8(seq).unwrap()
//     }

//     #[test]
//     fn test_diff_seed() {
//         let read = random_dna_seq(1000);
//         let read = read.as_bytes();
//         let k = 31;
//         let m = 11;
//         let hasher = Xxh3Builder::new().with_seed(0 as u64);
//         let minimizer_iter: CanonicalMinimizerIterator<'_, u64, Xxh3Builder> =
//             CanonicalMinimizerIterator::new(read, m, (k - m + 1) as u16, hasher, ENCODING);
//         let minimizers = minimizer_iter.collect_vec();
//         for seed in 1..100 {
//             let hasher = Xxh3Builder::new().with_seed(seed as u64);
//             let minimizer_iter_seed: CanonicalMinimizerIterator<'_, u64, Xxh3Builder> =
//                 CanonicalMinimizerIterator::new(read, m, (k - m + 1) as u16, hasher, ENCODING);
//             let minimizers_seed = minimizer_iter_seed.collect_vec();
//             assert!(minimizers != minimizers_seed);
//         }
//     }

//     // #[test]
//     // fn todo_remove() {
//     //     // TODO remove
//     //     let read = random_dna_seq(1000);
//     //     let read = read.as_bytes();
//     //     let k = 31;
//     //     let m = 11;
//     //     let hasher = DefaultHashBuilder::with_seed(0);
//     //     let minimizer_iter: CanonicalMinimizerPosIterator<'_> =
//     //         CanonicalMinimizerPosIterator::new(read, m, (k - m + 1) as u16, hasher, ENCODING);
//     //     let minimizers = minimizer_iter.collect_vec();
//     //     // println!("{:?}", minimizers);
//     //     // for seed in 1..100 {
//     //     let hasher = DefaultHashBuilder::with_seed(1);
//     //     let minimizer_iter_seed: CanonicalMinimizerPosIterator<'_, u64> =
//     //         CanonicalMinimizerPosIterator::new(read, m, (k - m + 1) as u16, hasher, ENCODING);
//     //     let minimizers_seed = minimizer_iter_seed.collect_vec();
//     //     assert!(minimizers != minimizers_seed);
//     //     // }
//     // }
// }
