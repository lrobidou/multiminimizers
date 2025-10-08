use needletail::{FastxReader, parse_fastx_file};
use std::path::Path;

pub struct LinesIter {
    data: Box<dyn FastxReader>,
}

impl LinesIter {
    pub fn new(data: Box<dyn FastxReader>) -> Self {
        Self { data }
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> Self {
        let sequences = parse_fastx_file(path).unwrap();
        Self::new(sequences)
    }
}

impl Iterator for LinesIter {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        let record = self.data.next()?.unwrap();
        let sequence = record.seq();
        let sequence_vec = sequence.iter().copied().collect(); // TODO copy
        Some(sequence_vec)
    }
}
