use itertools::Itertools;
use xxhash_rust::const_xxh3::xxh3_64;

// TODO test

use crate::superkmer::{
    reverse_complement, reverse_complement_ascii_to_ascii, reverse_complement_no_copy,
};

#[cfg(debug_assertions)]
use super::superkmer::is_canonical;

// Branch prediction hint. This is currently only available on nightly but it
// consistently improves performance by 10-15%.
#[cfg(debug_assertions)]
#[cfg(not(feature = "nightly"))]
use core::convert::identity as unlikely;
#[cfg(debug_assertions)]
#[cfg(feature = "nightly")]
use core::intrinsics::unlikely;

// states of Subsequence
#[derive(Debug, Clone, PartialEq)]
pub struct BitPacked<'a> {
    total_base_in_sequence: usize,
    read: &'a [u64],
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NoBitPacked<'a> {
    read: &'a [u8],
}

/// Represents a subsequence, possibly in reverse
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Subsequence<Packing> {
    start: usize,
    end: usize,
    same_orientation: bool,
    packing: Packing,
}

impl<'a> Subsequence<NoBitPacked<'a>> {
    pub fn new(read: &'a [u8], start: usize, end: usize, same_orientation: bool) -> Self {
        debug_assert!(start <= read.len());
        debug_assert!(end <= read.len());
        debug_assert!(start <= end);
        Self {
            start, // in base
            end,   // in base
            same_orientation,
            packing: NoBitPacked { read },
        }
    }

    pub fn get_read(&self) -> &[u8] {
        self.packing.read
    }

    pub fn get_subsequence_as_in_read(&self) -> &[u8] {
        &self.packing.read[self.start..self.end]
    }

    #[cfg(debug_assertions)]
    pub fn is_canonical(&self) -> bool {
        let subsequence = &self.packing.read[self.start..self.end];
        let is_original_subsequence_canonical = is_canonical(subsequence);
        self.same_orientation == is_original_subsequence_canonical
    }

    #[cfg(debug_assertions)]
    pub fn is_equal_to_its_revcomp(&self) -> bool {
        use crate::superkmer::is_equal_to_its_revcomp;

        let subsequence = &self.packing.read[self.start..self.end];
        is_equal_to_its_revcomp(subsequence)
    }

    #[cfg(debug_assertions)]
    pub fn to_canonical(self) -> Self {
        let subsequence = &self.packing.read[self.start..self.end];

        if is_canonical(subsequence) == self.same_orientation {
            self
        } else {
            Self {
                start: self.start,
                end: self.end,
                same_orientation: !self.same_orientation,
                packing: self.packing,
            }
        }
    }

    pub fn hash(&self) -> u64 {
        let subsequence = &self.packing.read[self.start..self.end];
        if self.same_orientation {
            xxh3_64(subsequence)
        } else {
            let revcomp = reverse_complement_no_copy(subsequence.iter().copied());
            // TODO copy into a vec
            xxh3_64(&revcomp.collect_vec())
        }
    }
}

impl std::fmt::Display for Subsequence<NoBitPacked<'_>> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let subsequence = &self.packing.read[self.start..self.end];
        let string = if self.same_orientation {
            String::from_utf8(Vec::from(subsequence)).unwrap()
        } else {
            reverse_complement(subsequence)
        };
        write!(f, "{}", string)
    }
}

impl Subsequence<NoBitPacked<'_>> {
    pub fn as_vec(&self) -> Vec<u8> {
        let subsequence = &self.packing.read[self.start..self.end];
        if self.same_orientation {
            Vec::from(subsequence)
        } else {
            reverse_complement_ascii_to_ascii(subsequence)
        }
    }
}

impl<'a> Subsequence<BitPacked<'a>> {
    pub fn whole_bitpacked(read: &'a [u64], nb_bases: usize) -> Self {
        #[cfg(debug_assertions)]
        {
            // let read = lock.get_slice_from_internal_id(id);
            debug_assert!((nb_bases / 32) + ((nb_bases % 32 != 0) as usize) == read.len());
        }
        Self {
            start: 0,
            end: nb_bases,
            same_orientation: true,
            packing: BitPacked {
                total_base_in_sequence: nb_bases,
                read,
            },
        }
    }

    // pub fn decode_2bits(&self) -> Vec<u8> {
    //     let decoder = Decoder::new(self.packing.read, self.packing.total_base_in_sequence);
    //     let ascii_whole = decoder.collect_vec();
    //     let ascii = ascii_whole[self.start..self.end]
    //         .iter()
    //         .copied()
    //         .collect_vec();
    //     ascii
    // }
}

impl<Packing> Subsequence<Packing> {
    pub fn change_orientation(self) -> Self {
        Self {
            start: self.start,
            end: self.end,
            same_orientation: !self.same_orientation,
            packing: self.packing,
        }
    }

    pub fn change_orientation_if(self, cond: bool) -> Self {
        if cond {
            Self {
                start: self.start,
                end: self.end,
                same_orientation: !self.same_orientation,
                packing: self.packing,
            }
        } else {
            self
        }
    }

    pub fn start(&self) -> usize {
        self.start
    }

    pub fn end(&self) -> usize {
        self.end
    }

    pub fn same_orientation(&self) -> bool {
        self.same_orientation
    }

    pub fn len(&self) -> usize {
        self.end - self.start
    }

    /// Extract a subsequence
    /// When not bit packed, equivalent to:
    /// ```
    /// Self::whole_string(self.to_string()[start..end])
    /// ```
    pub fn subsequence(self, start: usize, end: usize) -> Subsequence<Packing> {
        if self.same_orientation {
            Self {
                start: self.start + start,
                end: self.start + end,
                same_orientation: self.same_orientation,
                packing: self.packing,
            }
        } else {
            Self {
                start: self.end - end,
                end: self.end - start,
                same_orientation: self.same_orientation,
                packing: self.packing,
            }
        }
    }
}

#[cfg(debug_assertions)]
fn prefix_str(a: &str, b: &str) -> usize {
    let a = a.as_bytes();
    let b = b.as_bytes();
    iter_prefix_len(a.iter().copied(), b.iter().copied())
}

#[cfg(debug_assertions)]
fn suffix_str(a: &str, b: &str) -> usize {
    let a = a.as_bytes();
    let b = b.as_bytes();
    iter_suffix_len(&mut a.iter().copied(), &mut b.iter().copied())
}

#[cfg(debug_assertions)]
fn iter_prefix_len(mut x: impl Iterator<Item = u8>, mut y: impl Iterator<Item = u8>) -> usize {
    let mut length = 0;
    while let (Some(xc), Some(yc)) = (x.next(), y.next()) {
        // N is treated like a A
        let xc = if unlikely(xc == b'N') { b'A' } else { xc };
        let yc = if unlikely(yc == b'N') { b'A' } else { yc };

        if xc == yc {
            length += 1;
        } else {
            break;
        }
    }
    length
}

#[cfg(debug_assertions)]
fn iter_suffix_len(
    x: &mut impl DoubleEndedIterator<Item = u8>,
    y: &mut impl DoubleEndedIterator<Item = u8>,
) -> usize {
    let mut x = x.rev();
    let mut y = y.rev();
    let mut length = 0;
    while let (Some(xc), Some(yc)) = (x.next(), y.next()) {
        // N is treated like a A
        let xc = if unlikely(xc == b'N') { b'A' } else { xc };
        let yc = if unlikely(yc == b'N') { b'A' } else { yc };

        if xc == yc {
            length += 1;
        } else {
            break;
        }
    }
    length
}
