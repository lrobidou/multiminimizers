use clap::Parser;
use itertools::Itertools;
use std::collections::HashSet;
use std::path::Path;

use sticky_mini::superkmers_computation::compute_all_superkmers_linear_streaming as compute_all_superkmers;
use sticky_mini::superkmers_computation::compute_superkmers_linear_streaming as compute_superkmers;

mod needless_modules;
use crate::needless_modules::LinesIter;
use crate::needless_modules::QuotientSet;

fn replace_n(sequences: &mut Vec<Vec<u8>>) {
    for sequence in sequences {
        for char in sequence {
            if char == &b'N' {
                *char = b'A';
            }
        }
    }
}

const CANONICAL: bool = true;

/// Simple program to parse a list of integers
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    #[clap(short, long)]
    index: String,
    #[clap(short, long)]
    query: String,
    #[clap(short, long)]
    k: usize,
    #[clap(short, long)]
    m: usize,
}

struct Index<const N: usize> {
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

    pub fn search_ascii_read(&mut self, mut read: Vec<u8>) -> Vec<bool> {
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
}

fn main() {
    let args = Args::parse();
    let sequences = LinesIter::from_path(args.index);

    let mut index: Index<8> = Index::new(args.k, args.m);
    index.add_reads(sequences);

    LinesIter::from_path(args.query).for_each(|read| {
        let response = index.search_ascii_read(read);
        println!("{:?}", response);
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index_same_query() {
        let to_index = "tests/unit_test_data/index.fa";

        let k = 31;
        let m = 21;
        let mut index: Index<8> = Index::new(k, m);

        // index
        let sequences = LinesIter::from_path(to_index);
        index.add_reads(sequences);

        // query
        let responses = LinesIter::from_path(to_index)
            .map(|read| index.search_ascii_read(read))
            .collect_vec();

        // only one response and of the correct length
        assert_eq!(responses.len(), 1);
        assert_eq!(responses[0].len(), 100 - k + 1);

        // most new bases are not found
        assert!(responses[0].iter().all(|x| *x));
    }

    #[test]
    fn test_index_different_query() {
        let to_index = "tests/unit_test_data/index.fa";
        let to_query = "tests/unit_test_data/query.fa";

        let k = 31;
        let m = 21;
        let mut index: Index<8> = Index::new(k, m);

        // index
        let sequences = LinesIter::from_path(to_index);
        index.add_reads(sequences);

        // query
        let responses = LinesIter::from_path(to_query)
            .map(|read| index.search_ascii_read(read))
            .collect_vec();

        // only one response and of the correct length
        assert_eq!(responses.len(), 1);
        assert_eq!(responses[0].len(), 110 - k + 1);

        // most new bases are not found
        assert!(responses[0].iter().take(10).filter(|x| **x).count() < 5);
        assert!(responses[0].iter().skip(10).all(|x| *x));
    }
}
