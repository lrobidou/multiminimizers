use num_traits::PrimInt;

use super::superkmer::Superkmer;
use crate::simd_canonical_sticky_minimizer_iterator::{
    CanonicalStickyMinimizerIteratorSIMD, SKInfos,
};
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

pub struct SuperkmerIterator<'a, const N: usize, const CANONICAL: bool> {
    sequence: &'a [u8],
    minimizer_iter: std::iter::Peekable<CanonicalStickyMinimizerIteratorSIMD<N, CANONICAL>>,
    k: usize,
    m: usize,
    previous_minimizer: Option<SKInfos>,
}

impl<'a, const N: usize, const CANONICAL: bool> SuperkmerIterator<'a, N, CANONICAL> {
    fn new(
        sequence: &'a [u8],
        minimizer_iter: CanonicalStickyMinimizerIteratorSIMD<N, CANONICAL>,
        k: usize,
        m: usize,
    ) -> Self {
        let mut minimizer_iter = minimizer_iter.peekable();
        let previous_minimizer = minimizer_iter.next();
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

impl<'a, const N: usize, const CANONICAL: bool> Iterator for SuperkmerIterator<'a, N, CANONICAL> {
    type Item = Superkmer<'a, NoAnchor>;

    fn next(&mut self) -> Option<Self::Item> {
        let previous_minimizer = self.previous_minimizer?;
        let start_mini = previous_minimizer.minimizer_start;
        let start_sk = previous_minimizer.superkmer_start;
        let end_sk = previous_minimizer.superkmer_end;

        let sk = Superkmer::new(
            self.sequence,
            start_mini,
            start_mini + self.m,
            start_sk,
            end_sk,
            is_canonical(self.sequence, start_mini, self.m),
        );

        self.previous_minimizer = self.minimizer_iter.next();
        Some(sk)
    }
}

pub fn compute_superkmers_linear_streaming<'a, const N: usize, const CANONICAL: bool>(
    sequence: &'a [u8],
    k: usize,
    m: usize,
) -> Option<SuperkmerIterator<'a, N, CANONICAL>> {
    const {
        assert!(N >= 1, "At least one hash function is required");
    }
    if sequence.len() < k {
        None
    } else {
        let minimizer_iter: CanonicalStickyMinimizerIteratorSIMD<N, CANONICAL> =
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
        if compute_superkmers_linear_streaming::<1, true>("AGCAGCTAGCATTTT".as_bytes(), 16, 5)
            .is_some()
        {
            panic!()
        }

        if compute_superkmers_linear_streaming::<1, false>("AGCAGCTAGCATTTT".as_bytes(), 16, 5)
            .is_some()
        {
            panic!()
        }
    }

    fn random_dna_seq(len: usize) -> String {
        use rand::prelude::IndexedRandom;

        let bases = [b'A', b'C', b'G', b'T'];
        let mut rng = rand::rng();
        let seq: Vec<u8> = (0..len).map(|_| *bases.choose(&mut rng).unwrap()).collect();
        String::from_utf8(seq).unwrap()
    }

    fn test_all_kmers_covered_inner<const CANONICAL: bool>() {
        for _ in 0..100 {
            let size_read = 1000;
            let k = 31;
            let m = 17;
            let read = random_dna_seq(size_read);

            // used to break my algo computing the end of superkmer
            // let read = "AGAACAGTTATCGTCTGATTCGCGAGGTACGGTTGGACACGGTAGGCTTGGTCGAAA";

            let sk_iter =
                compute_superkmers_linear_streaming::<2, CANONICAL>(read.as_bytes(), k, m).unwrap();

            let sks = sk_iter.collect_vec();
            let mut curr_sk = 0;
            for i in 0..(read.len() - k + 1) {
                let sk = sks[curr_sk];
                assert!(sk.superkmer.start() < sk.superkmer.end() - k + 1);
                if (sk.superkmer.start()..sk.superkmer.end() - k + 1).contains(&i) {
                    // ok
                } else {
                    let next_sk = sks[curr_sk + 1];
                    assert!(
                        (next_sk.superkmer.start()..next_sk.superkmer.end() - k + 1).contains(&i)
                    );
                    curr_sk += 1;
                }
            }
        }
    }

    #[test]
    fn test_all_kmers_covered() {
        test_all_kmers_covered_inner::<false>();
        test_all_kmers_covered_inner::<true>();
    }

    fn check_density<'a>(m: usize, sk_iter: impl IntoIterator<Item = Superkmer<'a, NoAnchor>>) {
        let sks = sk_iter.into_iter().collect_vec();
        let read_size = sks[0].get_read().len();
        let nb_minimizer = sks.len();
        let sum_of_gaps: usize = sks
            .iter()
            .tuple_windows()
            .map(|(a, b)| b.start_of_minimizer() - a.start_of_minimizer())
            .sum();

        let avg_sum_of_gaps = sum_of_gaps as f64 / sks.len() as f64;
        let density = nb_minimizer as f64 / (read_size - m) as f64;

        // let's check the product of the two are roughly 1
        let fraction = density * avg_sum_of_gaps;
        let target = 1.0;
        let epsilon = 0.01;
        assert!((target - epsilon..target + epsilon).contains(&fraction));
    }
    #[test]
    fn test_product_density() {
        let size_read = 10000;
        let k = 31;
        let m = 19;
        let read = random_dna_seq(size_read);

        let sk_iter =
            compute_superkmers_linear_streaming::<1, true>(read.as_bytes(), k, m).unwrap();
        check_density(m, sk_iter);
    }
}
