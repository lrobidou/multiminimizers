![tests](https://github.com/lrobidou/sticky-mini/workflows/tests/badge.svg)
[![license](https://img.shields.io/badge/license-AGPL-purple)](https://github.com/lrobidou/sticky-mini//blob/main/LICENSE)


# What

Computes the "sticky minimizers". TODO change name

# How

## Plot the graph
```bash
cargo run  # write output to `data.txt`
python analyse.py  # plots the graph
```

## Index and query
```bash
cargo run --bin needless -- index --input tests/unit_test_data/index.fa --output "index.needless" -k 31 -m 21
cargo run --bin needless -- query --index index.needless --fasta tests/unit_test_data/query.fa  # vec of booleans, one per k-mer in the query
``