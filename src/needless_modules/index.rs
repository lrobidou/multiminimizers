use itertools::Itertools;
use needletail::sequence;
use serde::Deserialize;
use serde::Serialize;
use std::fs::File;
use std::io::BufReader;
use std::io::BufWriter;
use std::path::Path;

use super::{LinesIter, QuotientSet};

use crate::compute_all_superkmers_linear_streaming as compute_all_superkmers;
use crate::compute_superkmers_linear_streaming as compute_superkmers;

const CANONICAL: bool = true;
const NB_HASH: usize = 8;

pub fn replace_n(sequences: &mut Vec<Vec<u8>>) {
    for sequence in sequences {
        for char in sequence {
            if char == &b'N' {
                *char = b'A';
            }
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct Index<const N: usize> {
    k: usize,
    m: usize,
    data: QuotientSet,
}

impl<const N: usize> Index<N> {
    pub fn new(k: usize, m: usize) -> Self {
        Self {
            k,
            m,
            data: QuotientSet::new(),
        }
    }

    pub fn add_reads(&mut self, sequences: LinesIter) {
        let mut lines = sequences.into_iter().collect_vec(); //TODO copy
        replace_n(&mut lines);
        for sequence in lines {
            let superkmers = match compute_superkmers::<N, CANONICAL>(&sequence, self.k, self.m) {
                Some(superkmers_iter) => superkmers_iter,
                None => continue,
            };
            superkmers.for_each(|sk| {
                self.data.insert(sk.get_minimizer_no_hashed());
            })
        }
    }

    pub fn search_ascii_read(&self, mut read: Vec<u8>) -> Vec<bool> {
        for char in read.iter_mut() {
            if char == &b'N' {
                *char = b'A';
            }
        }

        let mut response = vec![false; read.len() - self.k + 1];
        let superkmers = match compute_all_superkmers::<N, CANONICAL>(&read, self.k, self.m) {
            Some(superkmers_iter) => superkmers_iter,
            None => {
                return response;
            }
        };
        superkmers.for_each(|sk| {
            if self.data.contains(sk.get_minimizer_no_hashed()) {
                // hit: update the covered k-mers
                let start = sk.superkmer.start();
                let end = sk.superkmer.end() - self.k + 1;
                response[start..end].fill(true);
            }
        });
        response
    }

    pub fn dump<P: AsRef<Path>>(&self, filename: P) -> Result<(), Box<dyn std::error::Error>> {
        let file = File::create(filename)?;
        let buffer = BufWriter::with_capacity(9000000, file);
        bincode::serialize_into(buffer, &self)?;
        Ok(())
    }

    pub fn load<P: AsRef<Path>>(filename: P) -> Result<Self, Box<dyn std::error::Error>> {
        let file = File::open(filename)?;
        let buffer = BufReader::with_capacity(9000000, file);
        let decoded: Self = bincode::deserialize_from(buffer)?;
        Ok(decoded)
    }
}
