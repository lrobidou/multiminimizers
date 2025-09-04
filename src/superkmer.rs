use xxhash_rust::const_xxh3::xxh3_64;

use crate::Minimizer;
use crate::subsequence::NoBitPacked;
use crate::subsequence::Subsequence;
use crate::two_bits;

use num_traits::PrimInt;
use std::cmp::min;
#[cfg(debug_assertions)]
use std::iter::Copied;
use std::iter::Map;
use std::iter::Rev;
use std::marker::PhantomData;

// Branch prediction hint. This is currently only available on nightly but it
// consistently improves performance by 10-15%.
#[cfg(not(feature = "nightly"))]
use core::convert::identity as unlikely;
#[cfg(feature = "nightly")]
use core::intrinsics::unlikely;
#[cfg(debug_assertions)]
use std::slice::Iter;

pub const REVCOMP_TAB_CHAR: [char; 255] = {
    let mut tab = ['A'; 255];
    tab[b'A' as usize] = 'T';
    tab[b'T' as usize] = 'A';
    tab[b'C' as usize] = 'G';
    tab[b'G' as usize] = 'C';
    tab[b'N' as usize] = 'T';
    tab
};

pub fn reverse_complement(seq: &[u8]) -> String {
    seq.iter()
        .rev()
        .map(|base| unsafe { *REVCOMP_TAB_CHAR.get_unchecked(*base as usize) })
        .collect()
}

pub fn reverse_complement_ascii_to_ascii(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|base| unsafe { *REVCOMP_TAB_CHAR.get_unchecked(*base as usize) } as u8)
        .collect()
}

#[cfg(debug_assertions)]
pub fn same_orientation(seq: &[u8]) -> Copied<Iter<'_, u8>> {
    seq.iter().copied()
}

#[cfg(debug_assertions)]
pub fn is_canonical(seq: &[u8]) -> bool {
    let mut orientation_1 = same_orientation(seq);
    let mut orientation_2 = reverse_complement_byte_iter(seq);
    while let (Some(xc), Some(yc)) = (orientation_1.next(), orientation_2.next()) {
        use std::cmp::Ordering;

        let xc = if unlikely(xc == b'N') { b'A' } else { xc };

        match xc.cmp(&yc) {
            Ordering::Less => return true,
            Ordering::Greater => return false,
            Ordering::Equal => {}
        }
    }
    // in case of palindrome, prefer saying the sequence is canonical
    true
}

// TODO should we do something else instead of using this complex type?
#[allow(clippy::type_complexity)]
#[cfg(debug_assertions)]
pub fn reverse_complement_byte_iter<'a>(seq: &'a [u8]) -> Map<Rev<Iter<'a, u8>>, fn(&'a u8) -> u8> {
    seq.iter()
        .rev()
        .map(|base| unsafe { *REVCOMP_TAB.get_unchecked(*base as usize) })
}

#[cfg(debug_assertions)]
pub fn is_equal_to_its_revcomp(seq: &[u8]) -> bool {
    let mut orientation_1 = same_orientation(seq);
    let mut orientation_2 = reverse_complement_byte_iter(seq);
    while let (Some(xc), Some(yc)) = (orientation_1.next(), orientation_2.next()) {
        let xc = if unlikely(xc == b'N') { b'A' } else { xc };

        if xc != yc {
            return false;
        }
    }
    true
}

// TODO duplication: move to a module
pub const REVCOMP_TAB: [u8; 255] = {
    let mut tab = [0; 255];
    tab[b'A' as usize] = b'T';
    tab[b'T' as usize] = b'A';
    tab[b'C' as usize] = b'G';
    tab[b'G' as usize] = b'C';
    tab[b'N' as usize] = b'T'; // N is A, so its complement is T
    tab
};

pub fn revcomp_ascii(ascii_base: u8) -> u8 {
    // TODO document unsafe
    // safe as u8 is smaller or equal to 255
    unsafe { *REVCOMP_TAB.get_unchecked(ascii_base as usize) }
}

// TODO "discuss" is there no copy here ? What is the cost of moving references ?
pub fn reverse_complement_no_copy(
    seq: impl DoubleEndedIterator<Item = u8>,
) -> Map<Rev<impl DoubleEndedIterator<Item = u8>>, fn(u8) -> u8> {
    seq.rev()
        .map(|base| unsafe { *REVCOMP_TAB.get_unchecked(base as usize) })
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Superkmer<'a, T: Copy> {
    minimizer: u64,
    start_mini: usize,
    end_mini: usize,
    anchor_infos: T,
    pub superkmer: Subsequence<NoBitPacked<'a>>,
}

#[derive(Copy, Clone, Debug)]
pub struct NoAnchor;

#[derive(Copy, Clone, Debug)]
pub struct AnchorInfos<T: PrimInt, const SIZE_HALF_G_BITS: usize> {
    _phantom_data: PhantomData<T>,
}

impl<T: PrimInt, const SIZE_HALF_G_BITS: usize> AnchorInfos<T, SIZE_HALF_G_BITS> {
    // const PASS: () = assert!(G <= size_of::<T>() * 8);
    // TODO find a way to make that fail to compile
    const PASS: () = assert!(false);
}

impl<'a> Superkmer<'a, NoAnchor> {
    pub fn new(
        read: &'a [u8],
        start_mini: usize,
        end_mini: usize,
        start_sk: usize,
        end_sk: usize,
        same_orientation: bool,
    ) -> Self {
        let minizer_subsequence = &read[start_mini..end_mini];
        let minimizer = if same_orientation {
            two_bits::encode_minimizer(minizer_subsequence.iter().copied())
        } else {
            two_bits::encode_minimizer(reverse_complement_no_copy(
                minizer_subsequence.iter().copied(),
            ))
        };

        let minimizer = xxh3_64(&minimizer.to_le_bytes());

        Self {
            minimizer,
            start_mini,
            end_mini,
            anchor_infos: NoAnchor,
            superkmer: Subsequence::<NoBitPacked>::new(read, start_sk, end_sk, same_orientation),
        }
    }
}

impl<'a, T: PrimInt, const SIZE_HALF_G_BITS: usize>
    Superkmer<'a, AnchorInfos<T, SIZE_HALF_G_BITS>>
{
    pub fn new_with_anchor(
        read: &'a [u8],
        start_mini: usize,
        end_mini: usize,
        start_sk: usize,
        end_sk: usize,
        same_orientation: bool,
    ) -> Self {
        const {
            assert!(SIZE_HALF_G_BITS <= size_of::<T>() * 8); // TODO redondant ?
        }
        let minizer_subsequence = &read[start_mini..end_mini];
        let minimizer = if same_orientation {
            two_bits::encode_minimizer(minizer_subsequence.iter().copied())
        } else {
            two_bits::encode_minimizer(reverse_complement_no_copy(
                minizer_subsequence.iter().copied(),
            ))
        };

        let minimizer = xxh3_64(&minimizer.to_le_bytes());

        Self {
            minimizer,
            start_mini,
            end_mini,
            anchor_infos: AnchorInfos {
                _phantom_data: PhantomData::<T>,
            },
            superkmer: Subsequence::<NoBitPacked>::new(read, start_sk, end_sk, same_orientation),
        }
    }
}

impl<'a, T: Copy> Superkmer<'a, T> {
    pub fn hash_superkmer(&self) -> u64 {
        self.superkmer.hash()
    }

    pub fn get_minimizer(&self) -> Minimizer {
        self.minimizer
    }

    pub fn start_of_minimizer(&self) -> usize {
        self.start_mini
    }

    pub fn end_of_minimizer(&self) -> usize {
        self.end_mini
    }

    pub fn is_canonical_in_the_read(&self) -> bool {
        self.superkmer.same_orientation()
    }

    pub fn get_read(&self) -> &[u8] {
        self.superkmer.get_read()
    }

    // #[cfg(any(debug_assertions, test))]
    pub fn minimizer_string(&self) -> String {
        use itertools::Itertools;
        if self.is_canonical_in_the_read() {
            String::from_utf8(
                self.get_read()[self.start_of_minimizer()..self.end_of_minimizer()]
                    .iter()
                    .copied()
                    .collect_vec(),
            )
            .unwrap()
        } else {
            reverse_complement(&self.get_read()[self.start_of_minimizer()..self.end_of_minimizer()])
        }
    }
}

/// An incomplete part of a minimizer
#[derive(Debug, PartialEq)]
pub enum IncompleteSide {
    Left { known_size: u8 },
    Right { known_size: u8 },
}

/// A minimizer reconstructed around an anchor.
/// It has:
/// - a base pointed by the anchor
/// - bases at the left of it
/// - bases at the right of it
/// - since left or right be incomplete, an enum indicating if that if the case
/// - the known size of the unknown part
#[derive(Debug, PartialEq)]
pub struct AnchorBasedMinimizer<T> {
    pub left: T,
    pub right: T,
    pub middle: u8,
    pub unknown_side: Option<IncompleteSide>,
}

const _: () = {
    assert!(size_of::<AnchorBasedMinimizer<u8>>() == 5);
    assert!(size_of::<AnchorBasedMinimizer<u16>>() == 8);
};

/// Splits the minimizer into "anchors context part"
/// Returns the left, middle and right part of the context
impl<T, const SIZE_HALF_G_BITS: usize> Superkmer<'_, AnchorInfos<T, SIZE_HALF_G_BITS>>
where
    T: PrimInt + num_traits::AsPrimitive<u8>,
{
    pub fn simulate_anchor(&self) -> AnchorBasedMinimizer<T> {
        // Strategy: we are given the size of the context in bits
        // - convert that to bases
        // - try to extract the bases from the superkmer
        // - remove one bit if the number of bit was odd

        let m = self.end_mini - self.start_mini;
        let half_m = (m - 1) / 2;

        // check i'm not too tired
        debug_assert_eq!(self.start_mini + half_m + 1 + half_m, self.end_mini);

        // check there is enough space to place g
        let (half_all_g, all_g, size_bits_dest) = {
            // number of plain bases in a part of the context
            let half_plain_g: usize = SIZE_HALF_G_BITS / 2;
            // number of plain bases of the context
            // let plain_g = half_plain_g * 2 + 1;

            // number of bases (incl. non plain) in a part of  the context
            let half_all_g = half_plain_g + (SIZE_HALF_G_BITS % 2);
            // number of bases (incl. non plain) in the context
            let all_g = half_all_g * 2 + 1;

            let size_bits_dest = std::mem::size_of::<T>() * 8;
            debug_assert!(size_bits_dest.is_multiple_of(2));
            debug_assert!(size_bits_dest >= (all_g - 1));
            assert!(size_bits_dest >= SIZE_HALF_G_BITS);
            (half_all_g, all_g, size_bits_dest)
        };

        // TODO use unsafe at the constructor to remove the checks ?
        debug_assert!((m % 2) == 1, "minimizers' size should be odd");
        assert!(m <= all_g, "g should be greater than m");

        let read = self.get_read();
        let anchor_pos_in_read = self.start_mini + half_m;

        let (middle_base, left, left_size_bits, right, right_size_bits) = if self
            .is_canonical_in_the_read()
        {
            let middle_base = ((read[anchor_pos_in_read]) >> 1) & 0b00000011;
            let left_size_all_base = min(anchor_pos_in_read - self.superkmer.start(), half_all_g);
            let right_size_all_base =
                min(self.superkmer.end() - anchor_pos_in_read - 1, half_all_g);

            let mut left = T::zero();
            // OPTIMIZE on peut sortir des termes ici
            (anchor_pos_in_read - left_size_all_base..anchor_pos_in_read)
                .rev()
                .enumerate()
                .for_each(|(i, position_in_read)| {
                    let base_2_bit: u8 = (read[position_in_read] & 0b00000110) >> 1;
                    let shift = size_bits_dest - 2 * (i + 1);
                    let shifted: T = T::from(base_2_bit).unwrap() << shift;
                    left = left + shifted; // TODO chelou
                });
            left = left >> (size_bits_dest - half_all_g * 2);

            let mut right = T::zero();
            (anchor_pos_in_read + 1..anchor_pos_in_read + 1 + right_size_all_base)
                .enumerate()
                .for_each(|(i, position_in_read)| {
                    let base_2_bit: u8 = (read[position_in_read] & 0b00000110) >> 1;
                    let shift = size_bits_dest - 2 * (i + 1);
                    let shifted = T::from(base_2_bit).unwrap() << shift;
                    right = right + shifted; // TODO chelou
                });
            right = right >> (size_bits_dest - half_all_g * 2);

            if !SIZE_HALF_G_BITS.is_multiple_of(2) {
                // if left_size_all_base == half_all_g {
                left = left >> 1;
                // }
                // if right_size_all_base == half_all_g {
                right = right >> 1;
                // }
            }

            let size_to_remove_left =
                if !SIZE_HALF_G_BITS.is_multiple_of(2) && left_size_all_base == half_all_g {
                    1
                } else {
                    0
                };
            let size_to_remove_right =
                if !SIZE_HALF_G_BITS.is_multiple_of(2) && right_size_all_base == half_all_g {
                    1
                } else {
                    0
                };
            (
                middle_base,
                left,
                left_size_all_base * 2 - size_to_remove_left,
                right,
                right_size_all_base * 2 - size_to_remove_right,
            )
        } else {
            let middle_base = (revcomp_ascii(read[anchor_pos_in_read]) >> 1) & 0b00000011;
            let right_size_all_base =
                std::cmp::min(anchor_pos_in_read - self.superkmer.start(), half_all_g);
            let left_size_all_base =
                std::cmp::min(self.superkmer.end() - anchor_pos_in_read - 1, half_all_g);

            let mut left = T::zero();
            (anchor_pos_in_read + 1..anchor_pos_in_read + 1 + left_size_all_base)
                .enumerate()
                .for_each(|(i, position_in_read)| {
                    // fetch the base and compute its revcomp (2 bit representation)
                    let base_2_bit: u8 = (revcomp_ascii(read[position_in_read]) & 0b00000110) >> 1;
                    // compute the shift needed to place the base in `left`
                    let shift = size_bits_dest - 2 * (i + 1);
                    // update left
                    // TODO test with .sum() on the iterator and check genereted code
                    let shifted = T::from(base_2_bit).unwrap() << shift;

                    left = left + shifted;
                });
            left = left >> (size_bits_dest - half_all_g * 2);

            let mut right = T::zero();
            (anchor_pos_in_read - right_size_all_base..anchor_pos_in_read)
                .rev()
                .enumerate()
                .for_each(|(i, position_in_read)| {
                    let base_2_bit: u8 = (revcomp_ascii(read[position_in_read]) & 0b00000110) >> 1;
                    let shift = size_bits_dest - 2 * (i + 1);
                    let shifted = T::from(base_2_bit).unwrap() << shift;
                    right = right + shifted;
                });
            right = right >> (size_bits_dest - half_all_g * 2);

            if !SIZE_HALF_G_BITS.is_multiple_of(2) {
                // if left_size_all_base == half_all_g {
                left = left >> 1;
                // }
                // if right_size_all_base == half_all_g {
                right = right >> 1;
                // }
            }

            let size_to_remove_left =
                if !SIZE_HALF_G_BITS.is_multiple_of(2) && left_size_all_base == half_all_g {
                    1
                } else {
                    0
                };
            let size_to_remove_right =
                if !SIZE_HALF_G_BITS.is_multiple_of(2) && right_size_all_base == half_all_g {
                    1
                } else {
                    0
                };
            (
                middle_base,
                left,
                left_size_all_base * 2 - size_to_remove_left,
                right,
                right_size_all_base * 2 - size_to_remove_right,
            )
        };

        let left_is_complete = left_size_bits == SIZE_HALF_G_BITS;
        let right_is_complete = right_size_bits == SIZE_HALF_G_BITS;

        let unknown_side = if !right_is_complete {
            Some(IncompleteSide::Right {
                known_size: right_size_bits as u8,
            })
        } else if !left_is_complete {
            Some(IncompleteSide::Left {
                known_size: left_size_bits as u8,
            })
        } else {
            None
        };

        // TODO only in debug
        if !right_is_complete && !left_is_complete {
            unreachable!();
        }

        AnchorBasedMinimizer {
            left,
            right,
            unknown_side,
            middle: middle_base,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::superkmer::IncompleteSide;

    #[test]
    fn test_reverse_complement() {
        let read = "ACTGAGCTA";
        let bytes = read.as_bytes();
        let revcomp = reverse_complement(bytes);
        assert_eq!(revcomp, String::from("TAGCTCAGT"));
    }

    #[test]
    fn test_reverse_complement_ascii_to_ascii() {
        let read = "ACTGAGCTA";
        let bytes = read.as_bytes();
        let revcomp = reverse_complement_ascii_to_ascii(bytes);
        assert_eq!(revcomp, String::from("TAGCTCAGT").as_bytes());
    }

    #[test]
    fn test_get_minimizer() {
        let read = "ACTGAGCTA";
        let bytes = read.as_bytes();
        let sk = Superkmer::new(bytes, 1, 4, 1, 6, true);
        assert_eq!(sk.minimizer_string(), String::from("CTG"));
        let sk = Superkmer::new(bytes, 1, 4, 1, 6, false);
        assert_eq!(sk.minimizer_string(), String::from("CAG"));

        // TODO weird case: the minimizer is allowed to be outside of the superkmer
        let sk = Superkmer::new(bytes, 1, 4, 2, 6, false);
        assert_eq!(sk.minimizer_string(), String::from("CAG"));
    }

    // #[test]
    // fn test_brisk_mask() {
    //     let read = "AGCAGCTAGCATGCATGCAACT".as_bytes();
    //     let sk: Superkmer<'_, AnchorInfos<u8, 5>> =
    //         Superkmer::new_with_anchor(read, 6, 9, 0, 22, true);
    //     // AGCAG {minimizer} ATGCATGCAACT
    //     // => brisk encoding = A  G  T  A  G  C  C  G  A  A  T  /  G  /  C  /  // A  /  A  /  C  /  T  /
    //     // =>                  00 11 10 00 11 01 01 11 00 00 10 00 11 00 01 00 // 00 00 00 00 01 00 10 00

    //     let expected_brisk_mask_32_bases: u64 =
    //         0b11111111_11111111_11111100_11001100_11001100_11001100_00000000_00000000;
    //     let encoding = sk.brisk_mask_32_bases();
    //     assert_eq!(encoding, expected_brisk_mask_32_bases);
    // }

    #[test]
    fn test_simulate_anchor_10_complete() {
        let read = "AGCAGCTAGCATGCATGCAACT".as_bytes();
        let sk: Superkmer<'_, AnchorInfos<u16, 10>> =
            Superkmer::new_with_anchor(read, 5, 10, 0, 22, true);
        // left       middle            right
        // AGCAGCT      A      GCATGCATGCAACT
        // left = CAGCT = 01 00 11 01 10 -> rev 1001110001
        // right = GCATG = 11 01 00 10 11 -> 1101001011
        let expected = AnchorBasedMinimizer {
            left: 0b0000001001110001,
            right: 0b0000001101001011,
            middle: 0,
            unknown_side: None,
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_10_left_incomplete() {
        let read = "AGCAGCTAGCATGCATGCAACT".as_bytes();
        let left_shift = 3;

        // shift data from previous test
        let len = read.len();
        let read = &read[left_shift..len];
        let sk: Superkmer<'_, AnchorInfos<u16, 10>> = Superkmer::new_with_anchor(
            read,
            5 - left_shift,
            10 - left_shift,
            0,
            22 - left_shift,
            true,
        );
        // left       middle            right
        // AGCT      A      GCATGCATGCAACT
        // left = AGCT = 00 11 01 10 -> rev 10011100
        // right = GCATG = 11 01 00 10 11 -> 1101001011
        let expected = AnchorBasedMinimizer {
            left: 0b0000001001110000,
            right: 0b0000001101001011,
            middle: 0,
            unknown_side: Some(IncompleteSide::Left { known_size: 8 }),
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_10_right_incomplete() {
        let read = "AGCAGCTAGCATGCATGCAACT".as_bytes();
        let right_shift = 12;

        // shift data from previous test
        let len = read.len();
        let read = &read[0..len - right_shift];
        let sk: Superkmer<'_, AnchorInfos<u16, 10>> =
            Superkmer::new_with_anchor(read, 5, 10, 0, 22 - right_shift, true);
        // left       middle            right
        // AGCAGCT      A      GC
        // left = CAGCT = 01 00 11 01 10 -> rev 1001110001
        // right = GC = 11 01 -> 1101
        let expected = AnchorBasedMinimizer {
            left: 0b0000001001110001,
            right: 0b0000001101000000,
            middle: 0,
            unknown_side: Some(IncompleteSide::Right { known_size: 4 }),
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_10_complete_revcomp() {
        let read = "AGTTGCATGCATGCTAGCTGCT".as_bytes();
        let sk: Superkmer<'_, AnchorInfos<u16, 10>> =
            Superkmer::new_with_anchor(read, 22 - 10, 22 - 5, 0, 22, false);
        // original read:
        // left       middle            right
        // AGCAGCT      A      GCATGCATGCAACT
        // in revcomp (this test)
        // left          middle    right
        // AGTTGCATGCATGC  T     AGCTGCT

        // left = CAGCT = 01 00 11 01 10 -> rev 1001110001
        // right = GCATG = 11 01 00 10 11 -> 1101001011
        let expected = AnchorBasedMinimizer {
            left: 0b0000001001110001,
            right: 0b0000001101001011,
            middle: 0,
            unknown_side: None,
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_10_left_incomplete_revcomp() {
        let read = "AGTTGCATGCATGCTAGCT".as_bytes();
        let sk: Superkmer<'_, AnchorInfos<u16, 10>> =
            Superkmer::new_with_anchor(read, 22 - 10, 22 - 5, 0, 22 - 3, false);
        // original read:
        // left       middle            right
        // AGCT      A      GCATGCATGCAACT
        // in revcomp (this test)
        // left          middle    right
        // AGTTGCATGCATGC  T     AGCT

        // revcomp read (this test)
        // AGTTGCATGCATGC T AGCTGCT

        // left = XAGCT = X 00 11 01 10 -> rev 1001110000
        // right = GCATG = 11 01 00 10 11 -> 1101001011
        let expected = AnchorBasedMinimizer {
            left: 0b0000001001110000,
            right: 0b0000001101001011,
            middle: 0,
            unknown_side: Some(IncompleteSide::Left { known_size: 8 }),
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_10_right_incomplete_revcomp() {
        let read = "AGTTGCATGCATGCTAGCTGCT".as_bytes();
        let right_shift = 12;

        // shift data from previous test
        let len = read.len();
        let read = &read[right_shift..len];
        let sk: Superkmer<'_, AnchorInfos<u16, 10>> = Superkmer::new_with_anchor(
            read,
            22 - right_shift - 10,
            22 - right_shift - 5,
            0,
            22 - right_shift,
            false,
        );
        // original read
        // left       middle            right
        // AGCAGCT      A      GC
        // in revcomp (this test)
        // left
        // GC T AGCTGCT

        // left = CAGCT = 01 00 11 01 10 -> rev 1001110001
        // right = GC = 11 01 -> 1101
        let expected = AnchorBasedMinimizer {
            left: 0b0000001001110001,
            right: 0b0000001101000000,
            middle: 0,
            unknown_side: Some(IncompleteSide::Right { known_size: 4 }),
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_7_complete() {
        let read = "AGCAGCTAGCATGCATGCAACT".as_bytes();
        let sk: Superkmer<'_, AnchorInfos<u8, 7>> =
            Superkmer::new_with_anchor(read, 5, 10, 0, 22, true);
        // left       middle            right
        // AGCAGCT      A      GCATGCATGCAACT
        // left = AGCT = 00 11 01 10 -> rev 10011100
        // right = GCAT = 11 01 00 10 -> 11010010
        let expected = AnchorBasedMinimizer {
            left: 0b01001110,
            right: 0b01101001,
            middle: 0,
            unknown_side: None,
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_7_left_incomplete() {
        let read = "AGCAGCTAGCATGCATGCAACT".as_bytes();
        let left_shift = 4;

        // shift data from previous test
        let len = read.len();
        let read = &read[left_shift..len];
        let sk: Superkmer<'_, AnchorInfos<u8, 7>> = Superkmer::new_with_anchor(
            read,
            5 - left_shift,
            10 - left_shift,
            0,
            22 - left_shift,
            true,
        );
        // left       middle            right
        // AGCT      A      GCATGCATGCAACT
        // left = XGCT = XX 11 01 10 -> rev 100111XX
        // right = GCATG = 11 01 00 10 11 -> 1101001011
        let expected = AnchorBasedMinimizer {
            left: 0b01001110,
            right: 0b01101001,
            middle: 0,
            unknown_side: Some(IncompleteSide::Left { known_size: 6 }),
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_7_right_incomplete() {
        let read = "AGCAGCTAGCATGCATGCAACT".as_bytes();
        let right_shift = 12;

        // shift data from previous test
        let len = read.len();
        let read = &read[0..len - right_shift];
        let sk: Superkmer<'_, AnchorInfos<u8, 7>> =
            Superkmer::new_with_anchor(read, 5, 10, 0, 22 - right_shift, true);
        // left       middle            right
        // AGCAGCT      A      GC
        // left = CAGCT = 01 00 11 01 10 -> rev 1001110001
        // right = GC = 11 01 -> 1101
        let expected = AnchorBasedMinimizer {
            left: 0b01001110,
            right: 0b01101000,
            middle: 0,
            unknown_side: Some(IncompleteSide::Right { known_size: 4 }),
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_7_complete_revcomp() {
        let read = "AGTTGCATGCATGCTAGCTGCT".as_bytes();
        let sk: Superkmer<'_, AnchorInfos<u8, 7>> =
            Superkmer::new_with_anchor(read, 22 - 10, 22 - 5, 0, 22, false);
        // original read:
        // left       middle            right
        // AGCAGCT      A      GCATGCATGCAACT
        // in revcomp (this test)
        // left          middle    right
        // AGTTGCATGCATGC  T     AGCTGCT

        // left = CAGCT = 01 00 11 01 10 -> rev 1001110001
        // right = GCATG = 11 01 00 10 11 -> 1101001011
        let expected = AnchorBasedMinimizer {
            left: 0b01001110,
            right: 0b01101001,
            middle: 0,
            unknown_side: None,
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_7_left_incomplete_revcomp() {
        let read = "AGTTGCATGCATGCTAGCT".as_bytes();
        let sk: Superkmer<'_, AnchorInfos<u8, 7>> =
            Superkmer::new_with_anchor(read, 22 - 10, 22 - 5, 0, 22 - 4, false);
        // original read:
        // left       middle            right
        // AGCT      A      GCATGCATGCAACT
        // in revcomp (this test)
        // left          middle    right
        // AGTTGCATGCATGC  T     AGCT

        // revcomp read (this test)
        // AGTTGCATGCATGC T AGCTGCT

        // left = XGCT = XX 11 01 10 -> rev 100111XX
        // right = GCAT = 11 01 00 10 -> 11010010
        let expected = AnchorBasedMinimizer {
            left: 0b01001110,
            right: 0b01101001,
            middle: 0,
            unknown_side: Some(IncompleteSide::Left { known_size: 6 }),
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_7_right_incomplete_revcomp() {
        let read = "AGTTGCATGCATGCTAGCTGCT".as_bytes();
        let right_shift = 12;

        // shift data from previous test
        let len = read.len();
        let read = &read[right_shift..len];
        let sk: Superkmer<'_, AnchorInfos<u8, 7>> = Superkmer::new_with_anchor(
            read,
            22 - right_shift - 10,
            22 - right_shift - 5,
            0,
            22 - right_shift,
            false,
        );
        // original read
        // left       middle            right
        // AGCAGCT      A      GC
        // in revcomp (this test)
        // left
        // GC T AGCTGCT

        // left = AGCT = 00 11 01 10 -> rev 10011100
        // right = GC = 11 01 -> 1101
        let expected = AnchorBasedMinimizer {
            left: 0b01001110,
            right: 0b01101000,
            middle: 0,
            unknown_side: Some(IncompleteSide::Right { known_size: 4 }),
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_8_complete() {
        let read = "AGCAGCTAGCATGCATGCAACT".as_bytes();
        let sk: Superkmer<'_, AnchorInfos<u8, 8>> =
            Superkmer::new_with_anchor(read, 5, 10, 0, 22, true);
        // left       middle            right
        // AGCAGCT      A      GCATGCATGCAACT
        // left = AGCT = 00 11 01 10 -> rev 10011100
        // right = GCAT = 11 01 00 10 -> 11010010
        let expected = AnchorBasedMinimizer {
            left: 0b10011100,
            right: 0b11010010,
            middle: 0,
            unknown_side: None,
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_8_left_incomplete() {
        let read = "AGCAGCTAGCATGCATGCAACT".as_bytes();
        let left_shift = 4;

        // shift data from previous test
        let len = read.len();
        let read = &read[left_shift..len];
        let sk: Superkmer<'_, AnchorInfos<u8, 8>> = Superkmer::new_with_anchor(
            read,
            5 - left_shift,
            10 - left_shift,
            0,
            22 - left_shift,
            true,
        );
        // left       middle            right
        // AGCT      A      GCATGCATGCAACT
        // left = AGCT = 00 11 01 10 -> rev 10011100
        // right = GCATG = 11 01 00 10 11 -> 1101001011
        let expected = AnchorBasedMinimizer {
            left: 0b10011100,
            right: 0b11010010,
            middle: 0,
            unknown_side: Some(IncompleteSide::Left { known_size: 6 }),
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_8_right_incomplete() {
        let read = "AGCAGCTAGCATGCATGCAACT".as_bytes();
        let right_shift = 12;

        // shift data from previous test
        let len = read.len();
        let read = &read[0..len - right_shift];
        let sk: Superkmer<'_, AnchorInfos<u8, 8>> =
            Superkmer::new_with_anchor(read, 5, 10, 0, 22 - right_shift, true);
        // left       middle            right
        // AGCAGCT      A      GC
        // left = CAGCT = 01 00 11 01 10 -> rev 1001110001
        // right = GC = 11 01 -> 1101
        let expected = AnchorBasedMinimizer {
            left: 0b10011100,
            right: 0b11010000,
            middle: 0,
            unknown_side: Some(IncompleteSide::Right { known_size: 4 }),
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_8_complete_revcomp() {
        let read = "AGTTGCATGCATGCTAGCTGCT".as_bytes();
        let sk: Superkmer<'_, AnchorInfos<u8, 8>> =
            Superkmer::new_with_anchor(read, 22 - 10, 22 - 5, 0, 22, false);
        // original read:
        // left       middle            right
        // AGCAGCT      A      GCATGCATGCAACT
        // in revcomp (this test)
        // left          middle    right
        // AGTTGCATGCATGC  T     AGCTGCT

        // left = CAGCT = 01 00 11 01 10 -> rev 1001110001
        // right = GCATG = 11 01 00 10 11 -> 1101001011
        let expected = AnchorBasedMinimizer {
            left: 0b10011100,
            right: 0b11010010,
            middle: 0,
            unknown_side: None,
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_8_left_incomplete_revcomp() {
        let read = "AGTTGCATGCATGCTAGCT".as_bytes();
        let sk: Superkmer<'_, AnchorInfos<u8, 8>> =
            Superkmer::new_with_anchor(read, 22 - 10, 22 - 5, 0, 22 - 4, false);
        // original read:
        // left       middle            right
        // AGCT      A      GCATGCATGCAACT
        // in revcomp (this test)
        // left          middle    right
        // AGTTGCATGCATGC  T     AGCT

        // revcomp read (this test)
        // AGTTGCATGCATGC T AGCTGCT

        // left = XGCT = X X 11 01 10 -> rev 10011100
        // right = GCAT = 11 01 00 10 -> 11010010
        let expected = AnchorBasedMinimizer {
            left: 0b10011100,
            right: 0b11010010,
            middle: 0,
            unknown_side: Some(IncompleteSide::Left { known_size: 6 }),
        };
        assert_eq!(expected, sk.simulate_anchor());
    }

    #[test]
    fn test_simulate_anchor_8_right_incomplete_revcomp() {
        let read = "AGTTGCATGCATGCTAGCTGCT".as_bytes();
        let right_shift = 12;

        // shift data from previous test
        let len = read.len();
        let read = &read[right_shift..len];
        let sk: Superkmer<'_, AnchorInfos<u8, 8>> = Superkmer::new_with_anchor(
            read,
            22 - right_shift - 10,
            22 - right_shift - 5,
            0,
            22 - right_shift,
            false,
        );
        // original read
        // left       middle            right
        // AGCAGCT      A      GC
        // in revcomp (this test)
        // left
        // GC T AGCTGCT

        // left = AGCT = 00 11 01 10 -> rev 10011100
        // right = GC = 11 01 -> 1101
        let expected = AnchorBasedMinimizer {
            left: 0b10011100,
            right: 0b11010000,
            middle: 0,
            unknown_side: Some(IncompleteSide::Right { known_size: 4 }),
        };
        assert_eq!(expected, sk.simulate_anchor());
    }
}
