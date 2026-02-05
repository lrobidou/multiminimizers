use multiminimizers::CanonicalStickyMinimizerIteratorSIMD;

use ahash::{HashMap, HashSet};
use serde::{Deserialize, Serialize};
use std::fs;

use rand::distr::weighted::WeightedIndex;
use rand::prelude::*;
use rand_distr::Geometric;

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

type Dict = HashMap<usize, Vec<f64>>;

#[derive(Serialize, Deserialize)]
struct DataToPlot {
    m: usize,
    w: u16,
    error_rates: Vec<f64>,
    data: Dict,
}

impl DataToPlot {
    pub fn new(m: usize, w: u16, data: Dict, error_rates: Vec<f64>) -> Self {
        Self {
            m,
            w,
            data,
            error_rates,
        }
    }

    pub fn write_to_file(&self, filename: &str) {
        let data_str = serde_json::to_string(&self).unwrap();
        fs::write(filename, data_str)
            .unwrap_or_else(|_| panic!("Should be able to write to `{filename}`"));
    }
}

fn random_seq(len: usize) -> Vec<u8> {
    let mut rng = rand::rng();
    let seq: Vec<u8> = (0..len).map(|_| *BASES.choose(&mut rng).unwrap()).collect();
    seq
}

fn mutate(seq: &[u8], subst: f64, ins: f64, del: f64) -> Vec<u8> {
    let mut res = Vec::with_capacity(seq.len() * 11 / 10);
    let mut rng = rand::rng();
    let p = subst + ins + del;
    let geom = Geometric::new(p).unwrap();
    let weighted = WeightedIndex::new([subst, ins, del]).unwrap();
    let mut i = 0;
    while i < seq.len() {
        let gap = geom.sample(&mut rng) as usize;
        let mut j = (i + gap).min(seq.len());
        res.extend_from_slice(&seq[i..j]);
        match weighted.sample(&mut rng) {
            0 => {
                res.push(*BASES.choose(&mut rng).unwrap());
                j += 1;
            }
            1 => {
                res.push(*BASES.choose(&mut rng).unwrap());
            }
            _ => {
                j += 1;
            }
        }
        i = j;
    }
    res
}

fn jaccard(a: &HashSet<&[u8]>, b: &HashSet<&[u8]>) -> f64 {
    let sa = a.len();
    let sb = b.len();
    let mut si = 0usize;
    a.iter().for_each(|&x| {
        if b.contains(x) {
            si += 1;
        }
    });
    si as f64 / (sa + sb - si) as f64
}

fn append_conservation<const N: usize>(m: usize, w: u16, data: &mut Dict, error_rates: &[f64]) {
    const CANONICAL: bool = true;
    let seq_len = 5_000_000;
    let seq = random_seq(seq_len);

    let iter = CanonicalStickyMinimizerIteratorSIMD::<N, CANONICAL>::new(&seq, m, w);
    let minis: HashSet<_> = iter
        .map(|info| {
            let i = info.minimizer_start;
            &seq[i..(i + m)]
        })
        .collect();

    let mut res = Vec::with_capacity(error_rates.len());
    for error_rate in error_rates {
        let seq = mutate(&seq, error_rate / 2., error_rate / 4., error_rate / 4.);
        let iter = CanonicalStickyMinimizerIteratorSIMD::<N, CANONICAL>::new(&seq, m, w);
        let minis_mut: HashSet<_> = iter
            .map(|info| {
                let i = info.minimizer_start;
                &seq[i..(i + m)]
            })
            .collect();
        res.push(jaccard(&minis, &minis_mut));
    }
    data.insert(N, res);
}

fn main() {
    let m = 21;
    let w = 15;
    let error_rates = vec![0.001, 0.0025, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035];
    let mut data = Dict::default();
    append_conservation::<1>(m, w, &mut data, &error_rates);
    append_conservation::<2>(m, w, &mut data, &error_rates);
    append_conservation::<3>(m, w, &mut data, &error_rates);
    append_conservation::<4>(m, w, &mut data, &error_rates);
    DataToPlot::new(m, w, data, error_rates).write_to_file("data_conservation.json");
}
