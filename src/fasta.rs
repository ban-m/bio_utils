
use serde::{Serialize,Deserialize};
use std::io;
use std::io::{BufRead, BufReader};
use std::path::Path;
#[derive(Debug)]
pub struct Reader<R: io::Read> {
    reader: BufReader<R>,
    line: String,
}

impl Reader<std::fs::File> {
    pub fn from_file<P: AsRef<Path>>(file: P) -> std::io::Result<Self> {
        let reader = std::fs::File::open(file).map(BufReader::new)?;
        let line = String::new();
        Ok(Self { reader, line })
    }
}

impl<R: io::Read> Reader<R> {
    pub fn new(reader: R) -> Self {
        let line = String::new();
        let reader = BufReader::new(reader);
        Self { reader, line }
    }
    pub fn read(&mut self, record: &mut Record) -> std::io::Result<usize> {
        record.clear();
        if self.line.is_empty() {
            self.reader.read_line(&mut self.line)?;
        }
        if !self.line.starts_with('>') {
            return Err(std::io::Error::from(std::io::ErrorKind::Other));
        } else {
            let mut header = self.line.split_whitespace();
            record.id = header.next().unwrap().trim_start_matches('>').to_string();
            record.desc = header.next().map(|e| e.to_string());
            loop {
                self.line.clear();
                self.reader.read_line(&mut self.line).unwrap();
                if self.line.starts_with('>') {
                    break;
                } else {
                    record.seq.push_str(&self.line);
                }
            }
            Ok(1)
        }
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
        match self.inner.read(&mut record) {
            Ok(_) if record.is_empty() => None,
            Ok(_) => Some(Ok(record)),
            Err(why) => Some(Err(why)),
        }
    }
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct Record {
    id: String,
    desc: Option<String>,
    seq: String,
}

impl Record {
    fn clear(&mut self) {
        self.id.clear();
        self.desc = None;
        self.seq.clear();
    }
    pub fn is_empty(&self) -> bool {
        self.id.is_empty() && self.seq.is_empty()
    }
    pub fn len(&self) -> usize {
        self.seq().len()
    }
    pub fn id(&self) -> &str {
        &self.id
    }
    pub fn seq(&self) -> &[u8] {
        self.seq.as_bytes()
    }
    pub fn desc(&self) -> &Option<String> {
        &self.desc
    }
}

impl std::fmt::Display for Record {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if let Some(ref desc) = self.desc {
            writeln!(f, ">{} {}\n{}", self.id, desc, self.seq)
        } else {
            writeln!(f, ">{}\n{}", self.id, self.seq)
        }
    }
}

/// Fastest method to open and parse fasta file.
pub fn parse_into_vec<P: AsRef<Path>>(file: P) -> std::io::Result<Vec<Record>> {
    let lines = std::fs::read_to_string(file)?;
    let mut result = Vec::with_capacity(bytecount::count(lines.as_bytes(), b'>'));
    let mut lines = lines.lines();
    let mut line = lines.next().unwrap();
    'outer: loop {
        let mut record = Record::default();
        let mut header = line[1..].splitn(2, ' ');
        record.id = header.next().unwrap().to_owned();
        record.desc = header.next().map(|e| e.to_owned());
        while let Some(next) = lines.next() {
            if next.starts_with('>') {
                line = next;
                break;
            } else {
                record.seq.push_str(next);
            }
        }
        if !record.seq.is_empty() {
            result.push(record);
        } else {
            return Ok(result);
        }
    }
}
