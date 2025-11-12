use std::collections::HashSet;

use serde::{Deserialize, Serialize};

const INDEX_BITS: usize = 10;
const INDEX_SIZE: usize = 1 << INDEX_BITS;

#[derive(Serialize, Deserialize)]
pub struct QuotientSet {
    buckets: Vec<HashSet<u64>>,
}

impl QuotientSet {
    pub fn new() -> Self {
        Self {
            buckets: vec![HashSet::new(); INDEX_SIZE],
        }
    }

    pub fn insert(&mut self, x: u64) {
        let q = (x >> (64 - INDEX_BITS)) as usize;
        let r = (x << INDEX_BITS) >> INDEX_BITS;
        self.buckets[q].insert(r);
    }

    pub fn contains(&self, x: u64) -> bool {
        let q = (x >> (64 - INDEX_BITS)) as usize;
        let r = (x << INDEX_BITS) >> INDEX_BITS;
        self.buckets[q].contains(&r)
    }
}
