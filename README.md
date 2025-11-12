![tests](https://github.com/lrobidou/multiminimizers/workflows/tests/badge.svg)
[![license](https://img.shields.io/badge/license-AGPL-purple)](https://github.com/lrobidou/multiminimizers//blob/main/LICENSE)

<!-- ## Plot the graph
```bash
cargo run  # write output to `data.txt`
python analyse.py  # plots the graph

``` -->

# Multiminimizers

Computes the multiminimizers scheme on a random minimizer scheme.

# How to use

## Cargo.toml
```toml
[dependencies]
multiminimizers = { git = "https://github.com/lrobidou/multiminimizers" }
```

## Public functions

The two main functions of this library are:
```rust
//  computes the superkmers selected by the multiminimizer scheme
use multiminimizers::compute_superkmers_linear_streaming;
//  computes all the candidates superkmers
use multiminimizers::compute_all_superkmers_linear_streaming;
```

## Usage

```rust
use multiminimizers::compute_superkmers_linear_streaming;

const CANONICAL: bool = true;
const NB_HASH: usize = 8;

fn f() {
    let sequence = String::from("ACTG");
    let bytes = sequence.as_bytes();
    let superkmers = match compute_superkmers::<NB_HASH, CANONICAL>(bytes, self.k, self.m) {
        Some(superkmers_iter) => superkmers_iter,
        None => /* handle error (no superkmer, is your sequence too short ?) */,
    };

    // now we can iterate over the superkmer of the sequence
    for superkmer in superkmers {
        // ...
    }
}
```

# Reproducing the figure in the paper

## Time plot
```bash
cargo criterion --message-format=json > results.jsonl  # generates JSON Lines
python3 benches/plot_bench.py results.jsonl # shows the graph "time wrt nb hash functions"
```

## Density plot
```bash
cargo run --release  # generates "data_fixed_w_{w}.json
cd density
cargo r -r -- eval -o density_4.json  # generates "density_4.json"
python3 plot-density.py  # loads every json files and save the plot to "density_4_w_{w}.pdf
```

<!-- ## Index and query
```bash
cargo run --bin needless -- index --input tests/unit_test_data/index.fa --output "index.needless" -k 31 -m 21
cargo run --bin needless -- query --index index.needless --fasta tests/unit_test_data/query.fa  # vec of booleans, one per k-mer in the query
`` -->