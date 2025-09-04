use num_traits::PrimInt;

use super::superkmer::Superkmer;
use crate::simd_canonical_sticky_minimizer_iterator::CanonicalStickyMinimizerIteratorSIMD;
use crate::superkmer::{AnchorInfos, NoAnchor, REVCOMP_TAB};
use std::cmp::Ordering;
use std::iter::{Copied, Map, Rev};
use std::marker::PhantomData;
use std::slice::Iter;

// Branch prediction hint. This is currently only available on nightly but it
// consistently improves performance by 10-15%.
#[cfg(not(feature = "nightly"))]
use core::convert::identity as unlikely;
#[cfg(feature = "nightly")]
use core::intrinsics::unlikely;

pub struct SuperkmerIterator<'a, const N: usize> {
    sequence: &'a [u8],
    minimizer_iter: std::iter::Peekable<CanonicalStickyMinimizerIteratorSIMD<N>>,
    k: usize,
    m: usize,
    previous_minimizer: Option<(u32, u32 /*, bool*/)>,
}

impl<'a, const N: usize> SuperkmerIterator<'a, N> {
    fn new(
        sequence: &'a [u8],
        minimizer_iter: CanonicalStickyMinimizerIteratorSIMD<N>,
        k: usize,
        m: usize,
    ) -> Self {
        let mut minimizer_iter = minimizer_iter.peekable();
        let previous_minimizer: Option<(u32, u32 /*, bool*/)> = minimizer_iter.next();
        Self {
            sequence,
            minimizer_iter,
            k,
            m,
            previous_minimizer,
        }
    }
}

// Canonical windows have >half TG characters.
fn is_canonical(sequence: &[u8], start_mini: usize, m: usize) -> bool {
    let is_t_g = |base: &&u8| **base == b'C' || **base == b'G';
    let mmer = &sequence[start_mini..start_mini + m];
    let count_t_g = mmer.iter().filter(is_t_g).count();
    count_t_g * 2 > m
}

impl<'a, const N: usize> Iterator for SuperkmerIterator<'a, N> {
    type Item = Superkmer<'a, NoAnchor>;

    fn next(&mut self) -> Option<Self::Item> {
        let previous_minimizer = self.previous_minimizer?;
        let start_mini = previous_minimizer.0 as usize;
        let start_sk = previous_minimizer.1 as usize;

        match self.minimizer_iter.next() {
            Some(current_minimizer) => {
                // println!(
                //     "[SuperkmerIterator] current minimizer is {:?}",
                //     previous_minimizer
                // );

                // println!(
                //     "[SuperkmerIterator] FYI, next minimizer will be  {:?}",
                //     current_minimizer
                // );

                let sk = Superkmer::new(
                    self.sequence,
                    start_mini,
                    start_mini + self.m,
                    start_sk,
                    current_minimizer.1 as usize + self.k - 1,
                    is_canonical(self.sequence, start_mini, self.m),
                );
                self.previous_minimizer = Some(current_minimizer);
                // println!("[SuperkmerIterator] will return  {:?}", sk);
                Some(sk)
            }
            None => {
                // println!("[SuperkmerIterator] no other minimizer");
                let sk = Superkmer::new(
                    self.sequence,
                    start_mini,
                    start_mini + self.m,
                    start_sk,
                    self.sequence.len(),
                    is_canonical(self.sequence, start_mini, self.m),
                );
                self.previous_minimizer = None;
                Some(sk)
            }
        }
    }
}

pub fn compute_superkmers_linear_streaming<'a, const N: usize>(
    sequence: &'a [u8],
    k: usize,
    m: usize,
) -> Option<SuperkmerIterator<'a, N>> {
    const {
        assert!(N >= 1, "At least one hash function is required");
    }
    if sequence.len() < k {
        None
    } else {
        let minimizer_iter: CanonicalStickyMinimizerIteratorSIMD<N> =
            CanonicalStickyMinimizerIteratorSIMD::new(sequence, m, (k - m + 1) as u16);
        let superkmer_iter = SuperkmerIterator::new(sequence, minimizer_iter, k, m);
        Some(superkmer_iter)
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::superkmer::reverse_complement_byte_iter;

    use super::*;

    #[test]
    fn test_reverse_complement() {
        let revcomp: Vec<u8> = reverse_complement_byte_iter("ACTGTGCAGTGCA".as_bytes()).collect();
        assert_eq!(revcomp, b"TGCACTGCACAGT");
    }

    #[test]
    fn test_reverse_complement_n() {
        let revcomp: Vec<u8> =
            reverse_complement_byte_iter("ACTGTGCAGTNNGNCA".as_bytes()).collect();
        // assert_eq!(revcomp, b"TG\0C\0\0ACTGCACAGT");
        assert_eq!(revcomp, b"TGTCTTACTGCACAGT");
    }

    #[test]
    fn test_compute_superkmers_sequence_too_short() {
        // there are no superkmers for sequence < k
        if let Some(_x) =
            compute_superkmers_linear_streaming::<1>("AGCAGCTAGCATTTT".as_bytes(), 16, 5)
        {
            panic!()
        }
    }
}
