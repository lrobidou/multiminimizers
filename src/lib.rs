pub mod needless_modules;
mod simd_canonical_sticky_minimizer_iterator;
mod subsequence;
pub mod superkmer;
pub mod superkmers_computation;
mod two_bits;

type Minimizer = u64;

pub use crate::superkmers_computation::compute_superkmers_linear_streaming;

// pub use lines_iter::LinesIter;
