use clap::{Args, Parser, Subcommand};

use itertools::Itertools;
use serde::Deserialize;
use serde::Serialize;
use std::fs::File;
use std::io::BufReader;
use std::io::BufWriter;
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
const NB_HASH: usize = 8;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Build an index indexing the k-mers of a FASTA/Q file
    #[clap(alias = "i")]
    Index(IndexArgs),
    /// Query an index
    #[clap(alias = "q")]
    Query(QueryArgs),
}

#[derive(Args, Debug)]
struct IndexArgs {
    /// K-mer size
    #[arg(short)]
    k: usize,
    /// Minimizer size
    #[arg(short, default_value_t = 21)]
    m: usize,
    /// Input file (FASTA/Q, possibly gzipped)
    #[arg(short, long)]
    input: String,
    /// Output file ({input}.needless by default)
    #[arg(short, long)]
    output: Option<String>,
}

#[derive(Args, Debug)]
struct QueryArgs {
    /// Input index
    #[arg(short, long)]
    index: String,
    /// Input fasta file
    #[arg(short, long)]
    fasta: String,
}

#[derive(Serialize, Deserialize)]
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

fn main() {
    let args = Cli::parse();
    match args.command {
        Command::Index(IndexArgs {
            k,
            m,
            input,
            output,
        }) => {
            let sequences = LinesIter::from_path(&input);

            let mut index: Index<NB_HASH> = Index::new(k, m);
            index.add_reads(sequences);
            let output = match output {
                Some(x) => x,
                None => format!("{}.needless", input),
            };
            index
                .dump(output)
                .expect("should have been able to dump the index");
        }
        Command::Query(QueryArgs { index, fasta }) => {
            let index: Index<NB_HASH> =
                Index::load(index).expect("should have been able to load index");

            LinesIter::from_path(fasta).for_each(|read| {
                let response = index.search_ascii_read(read);
                println!("{:?}", response);
            });
        }
    }
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
