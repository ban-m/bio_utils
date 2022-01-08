//! A very simple fastq reader.
use serde::{Deserialize, Serialize};
use std::io;
use std::io::{BufRead, BufReader};
use std::path::Path;
#[derive(Debug)]
pub struct Reader<R: io::Read> {
    reader: BufReader<R>,
    line: Vec<u8>,
}

impl Reader<std::fs::File> {
    pub fn from_file<P: AsRef<Path>>(file: P) -> std::io::Result<Self> {
        let reader = std::fs::File::open(file).map(BufReader::new)?;
        let line = Vec::new();
        Ok(Self { reader, line })
    }
}

impl<R: io::Read> Reader<R> {
    pub fn new(reader: R) -> Self {
        let line = Vec::new();
        let reader = BufReader::new(reader);
        Self { reader, line }
    }
    pub fn read(&mut self, record: &mut Record) -> std::io::Result<usize> {
        // Note that the fastq file is four lines each.
        self.line.clear();
        record.clear();
        self.reader.read_until(b'\n', &mut self.line)?;
        if self.line.is_empty() {
            return Ok(1);
        }
        if !self.line[0] == b'@' {
            return Err(std::io::Error::from(std::io::ErrorKind::Other));
        }
        self.line.pop().unwrap();
        record.id = String::from_utf8_lossy(&self.line[1..]).to_string();
        // Base
        self.line.clear();
        self.reader.read_until(b'\n', &mut self.line)?;
        assert!(record.seq.is_empty());
        self.line.pop().unwrap();
        record.seq.extend_from_slice(&self.line);
        // Empty
        self.line.clear();
        self.reader.read_until(b'\n', &mut self.line)?;
        // Quality.
        self.line.clear();
        self.reader.read_until(b'\n', &mut self.line)?;
        assert!(record.qual.is_empty());
        self.line.pop().unwrap();
        record.qual.extend_from_slice(&self.line);
        Ok(1)
    }
    pub fn records(self) -> Records<R> {
        Records { inner: self }
    }
}

#[derive(Debug)]
pub struct Records<R: io::Read> {
    inner: Reader<R>,
}

impl<R: io::Read> Iterator for Records<R> {
    type Item = std::io::Result<Record>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut record = Record::default();
        let result = self.inner.read(&mut record);
        match result {
            Ok(_) if record.is_empty() => None,
            Ok(_) => Some(Ok(record)),
            Err(why) => Some(Err(why)),
        }
    }
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct Record {
    id: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl Record {
    fn clear(&mut self) {
        self.id.clear();
        self.seq.clear();
        self.qual.clear();
    }
    pub fn new() -> Self {
        Self::default()
    }
    pub fn with_data(id: &str, seq: &[u8], qual: &[u8]) -> Self {
        let id = id.to_string();
        let seq = seq.to_vec();
        let qual = qual.to_vec();
        Self { id, seq, qual }
    }
    pub fn is_empty(&self) -> bool {
        self.id.is_empty() && self.seq.is_empty() && self.qual.is_empty()
    }
    pub fn len(&self) -> usize {
        self.seq().len()
    }
    pub fn id(&self) -> &str {
        &self.id
    }
    pub fn seq(&self) -> &[u8] {
        &self.seq
    }
    pub fn quality(&self) -> &[u8] {
        &self.qual
    }
}

impl std::fmt::Display for Record {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let seq = String::from_utf8_lossy(&self.seq);
        let qual = String::from_utf8_lossy(&self.qual);
        write!(f, "@{}\n{}\n+\n{}", self.id, seq, qual)
    }
}

/// Fastest method to open and parse fasta file.
pub fn parse_into_vec<P: AsRef<Path>>(file: P) -> std::io::Result<Vec<Record>> {
    let reader = std::fs::File::open(file).map(std::io::BufReader::new)?;
    parse_into_vec_from(reader)
}

pub fn parse_into_vec_from<R: io::Read>(reader: R) -> std::io::Result<Vec<Record>> {
    Ok(Reader::new(reader)
        .records()
        .filter_map(|e| e.ok())
        .fuse()
        .collect())
}
