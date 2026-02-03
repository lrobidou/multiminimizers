use clap::{Args, Parser, Subcommand};

use itertools::Itertools;
use multiminimizers::compute_superkmers_linear_streaming as compute_superkmers;
use multiminimizers::needless_modules::{LinesIter, replace_n};

pub mod constants {
    include!(concat!(env!("OUT_DIR"), "/constants.rs"));
}
use constants::NB_HASH;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Estimate the index size k-mers of a FASTA/Q file
    #[clap(alias = "i")]
    Estimate(EstimateArgs),
}

#[derive(Args, Debug)]
struct EstimateArgs {
    /// K-mer size
    #[arg(short)]
    k: usize,
    /// Minimizer size
    #[arg(short, default_value_t = 21)]
    m: usize,
    /// Input file (FASTA/Q, possibly gzipped)
    #[arg(short, long)]
    input: String,
}
fn main() {
    let args = Cli::parse();
    match args.command {
        Command::Estimate(EstimateArgs { k, m, input }) => {
            println!("using {NB_HASH} hash functions");
            let mut nb_hyperkmers = 0;

            let sequences = LinesIter::from_path(&input);

            let mut lines = sequences.into_iter().collect_vec(); //TODO copy
            replace_n(&mut lines);

            for sequence in lines {
                nb_hyperkmers += 1; // end of read 
                let superkmers = match compute_superkmers::<NB_HASH, true>(&sequence, k, m) {
                    Some(superkmers_iter) => superkmers_iter,
                    None => continue,
                };
                nb_hyperkmers += superkmers.count();
            }
            println!("estimation of number of hyperkmers: {nb_hyperkmers}");
            println!(
                "estimation of number of bases in hyperkmers: {}",
                nb_hyperkmers * k
            );
        }
    }
}
