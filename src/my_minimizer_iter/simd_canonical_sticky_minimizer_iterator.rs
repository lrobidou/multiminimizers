use minimizer_iter::iterator::CanonicalMinimizerIterator;
use simd_minimizers::seed::canonical_minimizer_and_superkmer_positions;

use core::hash::Hash;
use num_traits::{AsPrimitive, PrimInt};
use std::fmt::Debug;
use std::iter::Peekable;

use xxhash_rust::xxh3::Xxh3Builder;

struct MiniIterator {
    pos: usize,
    min_pos_vec: Vec<u32>,
    sks_pos_vec: Vec<u32>,
}

impl MiniIterator {
    fn new(min_pos_vec: Vec<u32>, sks_pos_vec: Vec<u32>) -> Option<Self> {
        if min_pos_vec.len() != sks_pos_vec.len() {
            return None;
        }
        Some(Self {
            pos: 0,
            min_pos_vec,
            sks_pos_vec,
        })
    }
}

impl Iterator for MiniIterator {
    type Item = (u32, u32); // (minimizer_start, superkmer_start)

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.min_pos_vec.len() {
            return None;
        }

        let res = (self.min_pos_vec[self.pos], self.sks_pos_vec[self.pos]);
        self.pos += 1;
        Some(res)
    }
}

/// Gets the further minimizer if it exists, else returns None.
fn get_max_pos<const N: usize>(arr: [Option<(u32, u32)>; N]) -> Option<usize> {
    arr.iter()
        .enumerate()
        .filter_map(|(i, opt)| opt.as_ref().map(|(a, _b)| (i, a)))
        .max_by_key(|&(_, second)| second)
        .map(|(i, _)| i)
    // arr.into_iter().flaltten().max_by_key(|x| x.1)
}

/// A SIMD iterator over `N` `CanonicalMinimizerIterator` of a sequence.
/// It requires an odd width to break ties between multiple minimizers.
pub struct CanonicalStickyMinimizerIteratorSIMD<const N: usize> {
    minimizer_iters: [Peekable<MiniIterator>; N],
    minimizers: [Option<(u32, u32)>; N],
    current_mini: Option<usize>,
}

impl<const N: usize> CanonicalStickyMinimizerIteratorSIMD<N> {
    pub fn new(seq: &[u8], minimizer_size: usize, width: u16) -> Self {
        assert_eq!(
            width % 2,
            1,
            "width must be odd to break ties between multiple minimizers"
        );
        // OPTIMIZE make an unsafe constructor for `CanonicalMinimizerIterator` that doesn't check if the width is odd ?
        let mut minimizer_iters: [Peekable<MiniIterator>; N] = std::array::from_fn(|i| {
            let mut min_pos_vec = vec![];
            let mut sks_pos_vec = vec![];
            // TODO
            canonical_minimizer_and_superkmer_positions(
                seq,
                minimizer_size, // TODO ask
                width as usize,
                i as u32,
                &mut min_pos_vec,
                &mut sks_pos_vec,
            );
            let minimizer_iters = MiniIterator::new(min_pos_vec, sks_pos_vec).unwrap();
            minimizer_iters.peekable()
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

    fn ensure_all_iters_cover_sk_starting_pos(&mut self, start_sk: u32) {
        // println!("start");
        for i in 0..self.minimizer_iters.len() {
            // println!("check {i}");
            while let Some((_next_start_mini, next_start_sk)) = self.minimizer_iters[i].peek()
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

impl<const N: usize> Iterator for CanonicalStickyMinimizerIteratorSIMD<N> {
    type Item = (u32, u32);

    fn next(&mut self) -> Option<Self::Item> {
        // early exits if nothing more to return
        let current_mini = self.current_mini?;
        let (mini_start, sk_start) = self.minimizers[current_mini]?;
        // move the current minimizer iterator

        self.minimizers[current_mini] = self.minimizer_iters[current_mini].next();

        // invariants might be broken, let's fix them
        if let Some((_start_mini, start_sk)) = self.minimizers[current_mini] {
            // fix invariants 2 and 3
            self.ensure_all_iters_cover_sk_starting_pos(start_sk);
            // fix invarient 1
            self.current_mini = get_max_pos(self.minimizers);
        }
        assert!(sk_start <= mini_start);
        Some((mini_start, sk_start))
    }
}

// #[cfg(test)]
// mod tests {
//     use itertools::Itertools;
//     use minimizer_iter::{DefaultHashBuilder, iterator::CanonicalMinimizerPosIterator};

//     use super::*;

//     #[test]
//     fn test_iter_single_kmer() {

//     }
// }
