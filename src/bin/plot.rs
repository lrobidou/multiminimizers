use multiminimizers::{
    superkmer::{NoAnchor, Superkmer},
    superkmers_computation::compute_superkmers_linear_streaming,
};
use serde::{Deserialize, Serialize};
use std::{collections::HashMap, fs};

use itertools::Itertools;

#[derive(Serialize, Deserialize)]
struct DataToPlotForAK {
    size_read: usize,
    k: usize,
    m: usize,
    limit_overhead: f64,
    overhead: Vec<f64>,
    limit_avg_superkmer_size: usize,
    avg_superkmer_size: Vec<f64>,
    limit_nb_superkmer: f64,
    nb_superkmer: Vec<usize>,
}

impl DataToPlotForAK {
    pub fn new(
        size_read: usize,
        k: usize,
        m: usize,
        overhead: Vec<f64>,
        avg_superkmer_size: Vec<f64>,
        nb_superkmer: Vec<usize>,
    ) -> Self {
        let limit_avg_superkmer_size = 2 * k - m;
        let limit_nb_superkmer = size_read as f64 / (k - m + 1) as f64;

        let limit_overhead = (2.0 * k as f64 - m as f64) / (k - m + 1) as f64;

        Self {
            size_read,
            k,
            m,
            limit_overhead,
            overhead,
            limit_avg_superkmer_size,
            avg_superkmer_size,
            limit_nb_superkmer,
            nb_superkmer,
        }
    }

    pub fn from_vec(size_read: usize, k: usize, m: usize, arr: Vec<(f64, f64, usize)>) -> Self {
        let mut overheads = vec![];
        let mut avgs_superkmer_size = vec![];
        let mut nbs_superkmer = vec![];

        for (overhead, avg_superkmer_size, nb_superkmer) in arr {
            overheads.push(overhead);
            avgs_superkmer_size.push(avg_superkmer_size);
            nbs_superkmer.push(nb_superkmer);
        }
        Self::new(
            size_read,
            k,
            m,
            overheads,
            avgs_superkmer_size,
            nbs_superkmer,
        )
    }
}

#[derive(Serialize, Deserialize)]
struct DataToPlot {
    data: HashMap<usize, DataToPlotForAK>,
    ns: Vec<usize>,
}

impl DataToPlot {
    pub fn new(data: HashMap<usize, DataToPlotForAK>, ns: Vec<usize>) -> Self {
        Self { data, ns }
    }

    pub fn compute_canonical_data_fixed_w(
        m_iter: impl Iterator<Item = usize>,
        read: &str,
        w: usize,
    ) -> Self {
        let ns = vec![2, 3, 4, 8, 16, 32];
        let raw_data_for_multiple_ks = m_iter
            .map(|m| {
                // w = k - m + 1
                let k = w + m - 1;
                let arr = vec![
                    canonical_get_stats::<2>(read, k, m),
                    canonical_get_stats::<3>(read, k, m),
                    canonical_get_stats::<4>(read, k, m),
                    canonical_get_stats::<8>(read, k, m),
                    canonical_get_stats::<16>(read, k, m),
                    canonical_get_stats::<32>(read, k, m),
                ];
                (k, DataToPlotForAK::from_vec(read.len(), k, m, arr))
            })
            .collect();

        Self::new(raw_data_for_multiple_ks, ns)
    }

    pub fn write_to_file(&self, filename: &str) {
        let data_str = serde_json::to_string(&self).unwrap();
        fs::write(filename, data_str)
            .unwrap_or_else(|_| panic!("Should be able to write to `{filename}`"));
    }
}

fn random_dna_seq(len: usize) -> String {
    use rand::prelude::IndexedRandom;

    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::rng();
    let seq: Vec<u8> = (0..len).map(|_| *bases.choose(&mut rng).unwrap()).collect();
    String::from_utf8(seq).unwrap()
}

fn nb_base_in_representation(v: &[Superkmer<'_, NoAnchor>]) -> usize {
    v.iter().map(|sk| sk.superkmer.len()).sum()
}

fn canonical_get_stats<const N: usize>(read: &str, k: usize, m: usize) -> (f64, f64, usize) {
    let superkmer_iter = compute_superkmers_linear_streaming::<N, true>(read.as_bytes(), k, m);
    let superkmers: Vec<Superkmer<'_, NoAnchor>> = superkmer_iter.unwrap().collect_vec();

    let nb_base = nb_base_in_representation(&superkmers);
    let overhead = nb_base as f64 / read.len() as f64;
    let avg_size = nb_base as f64 / superkmers.len() as f64;
    let nb_sk = superkmers.len();

    (overhead, avg_size, nb_sk)
}

fn main() {
    let size_read = 5_000_000;
    let read = random_dna_seq(size_read);
    let m_start = 5;
    let m_stop = 81;
    let step = 2;

    for w in [3, 7, 15, 25, 31, 63] {
        let ms = (m_start..=m_stop).step_by(step);
        let filename_fixed_w = format!("data_fixed_w_{w}.json");
        let data = DataToPlot::compute_canonical_data_fixed_w(ms, &read, w);
        data.write_to_file(&filename_fixed_w);
    }
}
