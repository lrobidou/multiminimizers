use std::collections::HashSet;

const INDEX_BITS: usize = 10;
const INDEX_SIZE: usize = 1 << INDEX_BITS;

pub struct QuotientSet {
    buckets: Vec<HashSet<u32>>,
}

impl QuotientSet {
    pub fn new() -> Self {
        Self {
            buckets: vec![HashSet::new(); INDEX_SIZE],
        }
    }

    pub fn insert(&mut self, x: u64) {
        let q = (x >> (64 - 32)) as u32;
        let i = ((x >> (64 - 42)) & 0x3FF) as usize;
        debug_assert!(
            (x & ((1u64 << (64 - 42)) - 1)) == 0,
            "Ignored bits are not zero"
        );
        self.buckets[i].insert(q);
    }

    pub fn contains(&self, x: u64) -> bool {
        let q = (x >> (64 - 32)) as u32;
        let i = ((x >> (64 - 42)) & 0x3FF) as usize;
        self.buckets[i].contains(&q)
    }
}
