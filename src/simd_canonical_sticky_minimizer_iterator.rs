use simd_minimizers::seeded::canonical_minimizer_and_superkmer_positions;
use simd_minimizers::seeded::minimizer_and_superkmer_positions;

use std::fmt::Debug;
use std::iter::Peekable;

#[derive(Copy, Clone, Debug)]
pub struct SKInfos {
    pub minimizer_start: usize,
    pub superkmer_start: usize,
    pub superkmer_end: usize,
}

struct MiniIterator {
    pos: usize,
    k: usize,
    size_read: usize,
    min_pos_vec: Vec<u32>,
    sks_pos_vec: Vec<u32>,
}

impl MiniIterator {
    fn new(
        k: usize,
        size_read: usize,
        min_pos_vec: Vec<u32>,
        sks_pos_vec: Vec<u32>,
    ) -> Option<Self> {
        if min_pos_vec.len() != sks_pos_vec.len() {
            return None;
        }

        Some(Self {
            k,
            size_read,
            pos: 0,
            min_pos_vec,
            sks_pos_vec,
        })
    }
}

impl Iterator for MiniIterator {
    type Item = SKInfos;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.min_pos_vec.len() {
            return None;
        }

        let infos = if self.pos + 1 < self.min_pos_vec.len() {
            SKInfos {
                minimizer_start: self.min_pos_vec[self.pos] as usize,
                superkmer_start: self.sks_pos_vec[self.pos] as usize,
                superkmer_end: self.sks_pos_vec[self.pos + 1] as usize + self.k - 1,
            }
        } else {
            SKInfos {
                minimizer_start: self.min_pos_vec[self.pos] as usize,
                superkmer_start: self.sks_pos_vec[self.pos] as usize,
                superkmer_end: self.size_read,
            }
        };

        self.pos += 1;
        Some(infos)
    }
}

/// Gets the further minimizer if it exists, else returns None.
fn get_max_pos<const N: usize>(arr: &[Option<SKInfos>; N]) -> Option<usize> {
    arr.iter()
        .enumerate()
        .filter_map(|(i, opt)| opt.map(|infos| (i, infos.minimizer_start)))
        .max_by_key(|&(_i, second)| second)
        .map(|(i, _second)| i)
}

/// A SIMD iterator over `N` `CanonicalMinimizerIterator` of a sequence.
/// It requires an odd width to break ties between multiple minimizers.
pub struct CanonicalStickyMinimizerIteratorSIMD<const N: usize, const CANONICAL: bool> {
    minimizer_iters: [Peekable<MiniIterator>; N],
    minimizers: [Option<SKInfos>; N],
    current_mini: Option<usize>,
    k: usize,
}

impl<const N: usize, const CANONICAL: bool> CanonicalStickyMinimizerIteratorSIMD<N, CANONICAL> {
    pub fn new(seq: &[u8], minimizer_size: usize, width: u16) -> Self {
        assert_eq!(
            width % 2,
            1,
            "width must be odd to break ties between multiple minimizers"
        );
        let k = minimizer_size + width as usize - 1;

        let mut minimizer_iters: [Peekable<MiniIterator>; N] = std::array::from_fn(|i| {
            let mut min_pos_vec = vec![];
            let mut sks_pos_vec = vec![];
            if CANONICAL {
                canonical_minimizer_and_superkmer_positions(
                    seq,
                    minimizer_size, // TODO ask
                    width as usize,
                    i as u32,
                    &mut min_pos_vec,
                    &mut sks_pos_vec,
                );
            } else {
                minimizer_and_superkmer_positions(
                    seq,
                    minimizer_size, // TODO ask
                    width as usize,
                    i as u32,
                    &mut min_pos_vec,
                    &mut sks_pos_vec,
                );
            }

            let minimizer_iters =
                MiniIterator::new(k, seq.len(), min_pos_vec, sks_pos_vec).unwrap();
            minimizer_iters.peekable()
        });

        let minimizers = std::array::from_fn(|i| minimizer_iters[i].next());
        let current_mini = get_max_pos(&minimizers);

        // the current minimizer is self.minimizers[self.current_mini?]?
        // invariants:
        // 1. the current minimizer is the farthest in the read
        // 2. all minimizers'superkmer covers the first kmer of the superkmer denoted by the current minimizer
        // 3. the next minimizer of all minimizers does not cover said kmer
        Self {
            minimizer_iters,
            minimizers,
            current_mini,
            k,
        }
    }

    fn ensure_all_iters_cover_sk_starting_pos(&mut self, start_sk: usize) {
        for i in 0..self.minimizer_iters.len() {
            while let Some(infos) = self.minimizer_iters[i].peek()
                && infos.superkmer_start <= start_sk
            {
                self.minimizers[i] = self.minimizer_iters[i].next()
            }
        }
    }
}

impl<const N: usize, const CANONICAL: bool> Iterator
    for CanonicalStickyMinimizerIteratorSIMD<N, CANONICAL>
{
    type Item = SKInfos;

    fn next(&mut self) -> Option<Self::Item> {
        // early exits if nothing more to return
        let current_mini = self.current_mini?;
        let infos = &self.minimizers[current_mini]?;

        // move the current minimizer iterator

        self.minimizers[current_mini] = self.minimizer_iters[current_mini].next();

        // invariants might be broken, let's fix them
        if self.minimizers[current_mini].is_some() {
            // fix invariants 2 and 3
            self.ensure_all_iters_cover_sk_starting_pos(infos.superkmer_end - self.k + 1);
            // fix invariant 1
            self.current_mini = get_max_pos(&self.minimizers);
        }
        Some(*infos)
    }
}
