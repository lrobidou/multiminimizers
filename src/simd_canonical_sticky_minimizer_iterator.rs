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

/// Gets the furthest superkmer if it exists, else returns None.
fn get_max_pos<const N: usize>(arr: &[Option<SKInfos>; N]) -> Option<usize> {
    arr.iter()
        .enumerate()
        .filter_map(|(i, opt)| opt.map(|infos| (i, infos.superkmer_end)))
        .max_by_key(|&(_i, second)| second)
        .map(|(i, _second)| i)
}

// /// Stores the `SKInfos` skipped when iterating over multiminimizers
// pub struct NoSPSS<const N: usize> {
//     /// Data for all hash functions
//     data: [Vec<SKInfos>; N],
// }

// /// skipped `SKInfos` are useless when iterating over SPSS => the SPSS struct is empty
// pub struct SPSS {}

/// A SIMD iterator over `N` `CanonicalMinimizerIterator` of a sequence.
/// It requires an odd width to break ties between multiple minimizers.
pub struct CanonicalStickyMinimizerIteratorSIMD<const N: usize, const CANONICAL: bool> {
    minimizer_iters: [Peekable<MiniIterator>; N],
    minimizers: [Option<SKInfos>; N],
    current_mini: Option<usize>,
    k: usize,
}

// trait WithCapacity {
//     fn with_capacity(capacity: usize) -> Self;
// }

// impl<const N: usize> WithCapacity for NoSPSS<N> {
//     /// Initializes the data for all hash functions with the given capacity.
//     fn with_capacity(capacity: usize) -> Self {
//         Self {
//             data: std::array::from_fn(|_| Vec::with_capacity(capacity)),
//         }
//     }
// }

// impl WithCapacity for SPSS {
//     fn with_capacity(_capacity: usize) -> Self {
//         Self {}
//     }
// }

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
            // skipped: S::with_capacity(k - minimizer_size),
            // skipped_old: S::with_capacity(k - minimizer_size),
        }
    }
}

impl<const N: usize, const CANONICAL: bool> CanonicalStickyMinimizerIteratorSIMD<N, CANONICAL> {
    fn ensure_all_iters_cover_sk_starting_pos(&mut self, start_sk: usize) {
        for i in 0..self.minimizer_iters.len() {
            while let Some(next_infos) = self.minimizer_iters[i].peek()
                && next_infos.superkmer_start <= start_sk
            {
                self.minimizers[i] = self.minimizer_iters[i].next()
            }
        }
    }
}

// impl<const N: usize, const CANONICAL: bool>
//     CanonicalStickyMinimizerIteratorSIMD<N, CANONICAL, NoSPSS<N>>
// {
//     fn ensure_all_iters_cover_starting_pos(
//         &mut self,
//         start_sk: usize,
//         // hash_i: usize,
//     ) {
//         // clear all data from previous iteration
//         self.skipped
//             .data
//             .iter_mut()
//             .for_each(|hash_vec| hash_vec.clear());

//         // let mut first = true;
//         // skip minimizer before the one selected
//         for i in 0..self.minimizer_iters.len() {
//             while let Some(next_infos) = self.minimizer_iters[i].peek()
//                 && next_infos.superkmer_start <= start_sk
//             {
//                 // if hash_i == i && first {
//                 //     first = false;
//                 // } else {

//                 // TODO should awlays be none ? unwrap ? perf cost ?
//                 // copy the skipped minimizers to be able to backtrack
//                 self.skipped.data[i].push(self.minimizers[i].unwrap());
//                 self.minimizers[i] = self.minimizer_iters[i].next()

//                 // }
//             }
//         }
//     }
// }

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

// impl<'a, const N: usize, const CANONICAL: bool>
//     CanonicalStickyMinimizerIteratorSIMD<N, CANONICAL, NoSPSS<N>>
// {
//     fn next(&'a mut self) -> Option<(SKInfos, &'a [Vec<SKInfos>; N])> {
//         // early exits if nothing more to return
//         let current_mini = self.current_mini?;
//         let infos = &self.minimizers[current_mini]?;
//         self.skipped_old.data = self.skipped.data.clone();

//         // move the current minimizer iterator
//         self.minimizers[current_mini] = self.minimizer_iters[current_mini].next();

//         // invariants might be broken, let's fix them
//         if self.minimizers[current_mini].is_some() {
//             // fix invariants 2 and 3

//             self.ensure_all_iters_cover_starting_pos(infos.superkmer_end - self.k + 1);
//             // fix invariant 1
//             self.current_mini = get_max_pos(&self.minimizers);
//         }
//         Some((*infos, &self.skipped_old.data))
//     }
// }

/// A SIMD iterator over `N` `CanonicalMinimizerIterator` of a sequence.
/// It requires an odd width to break ties between multiple minimizers.
/// Returns all the possible minimizers.
/// # Warning
/// This iterator can go back in the read, and a given minimizer can be covered multiple times.
pub struct AllCanonicalStickyMinimizerIteratorSIMD<const N: usize, const CANONICAL: bool> {
    minimizer_iters: [Peekable<MiniIterator>; N],
    current_mini: usize,
    end: [bool; N],
}

impl<const N: usize, const CANONICAL: bool> AllCanonicalStickyMinimizerIteratorSIMD<N, CANONICAL> {
    pub fn new(seq: &[u8], minimizer_size: usize, width: u16) -> Self {
        assert_eq!(
            width % 2,
            1,
            "width must be odd to break ties between multiple minimizers"
        );
        let k = minimizer_size + width as usize - 1;

        let minimizer_iters: [Peekable<MiniIterator>; N] = std::array::from_fn(|i| {
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

        let current_mini = 0;

        // the current minimizer is self.minimizers[self.current_mini?]?
        // invariants:
        // 1. the current minimizer is the farthest in the read
        // 2. all minimizers'superkmer covers the first kmer of the superkmer denoted by the current minimizer
        // 3. the next minimizer of all minimizers does not cover said kmer
        Self {
            minimizer_iters,
            current_mini,
            end: [false; N],
        }
    }
}

impl<const N: usize, const CANONICAL: bool> Iterator
    for AllCanonicalStickyMinimizerIteratorSIMD<N, CANONICAL>
{
    type Item = SKInfos;
    fn next(&mut self) -> Option<Self::Item> {
        let n = self.minimizer_iters.len();
        if n == 0 {
            return None;
        }

        let mut checked = 0;

        while checked < n {
            let i = self.current_mini;
            self.current_mini = (self.current_mini + 1) % n;

            if self.end[i] {
                checked += 1;
                continue;
            }

            match self.minimizer_iters[i].next() {
                Some(item) => return Some(item),
                None => {
                    self.end[i] = true; // mark iterator finished
                }
            }

            checked += 1;
        }

        None
    }
}
