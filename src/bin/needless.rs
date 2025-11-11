use clap::{Args, Parser, Subcommand};

use multiminimizers::needless_modules::{Index, LinesIter};

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
    use itertools::Itertools;

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
