// mod canonical_minimizer_iterator;
mod canonical_sticky_minimizer_iterator;
mod simd_canonical_sticky_minimizer_iterator;

// use canonical_minimizer_iterator::CanonicalMinimizerIterator;
// pub use canonical_sticky_minimizer_iterator::CanonicalStickyMinimizerIterator;
pub use simd_canonical_sticky_minimizer_iterator::CanonicalStickyMinimizerIteratorSIMD;

pub use canonical_sticky_minimizer_iterator::ENCODING;
