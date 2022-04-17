//! Multiple alignmnet format.
//! This is a quick implementation of MAF parser.
//! Should be refactored as soon as possible.
use std::fmt;
use std::fs;
use std::io;
use std::io::prelude::*;
use std::path::Path;
/// A MAF reader.
#[derive(Debug)]
pub struct Reader<R: io::Read> {
    reader: io::BufReader<R>,
    line: String,
}

impl Reader<fs::File> {
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(Self::new)
    }
}

impl<R: io::Read> Reader<R> {
    pub fn new(reader: R) -> Self {
        Reader {
            reader: io::BufReader::new(reader),
            line: String::new(),
        }
    }
    pub fn records(self) -> Records<R> {
        Records::new(self)
    }
    pub fn read(&mut self, record: &mut Record) -> io::Result<bool> {
        record.clear();
        if !self.line.trim().is_empty() {
            eprintln!("Error: the line is not empty:{}", self.line);
        }
        let mut is_first = true;
        loop {
            self.line.clear();
            let num_bytes = self.reader.read_line(&mut self.line)?;
            if num_bytes == 0 && is_first {
                // It reaches the EOF
                break Ok(false);
            } else if self.line.trim().is_empty() {
                // It reaches the boundary.
                break Ok(true);
            } else {
                is_first = false;
                record.add_line(&self.line);
            }
        }
    }
}
// impl<R:std::io::Read> Reader<R>{
//     pub fn show(&self){
//         for dur in self.time.iter(){
//             println!("TIME\t{}.{}",dur.as_secs(),dur.subsec_millis());
//         }
//     }
// }

#[derive(Debug, Clone, Default)]
pub struct Record {
    // Some field
    score: Option<f64>,
    pass: Option<u64>,
    header: Vec<(String, String)>,
    sequence: Vec<Seq>,
    sequence_index: usize,
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut res = String::new();
        res.push_str("a ");
        if let Some(score) = self.score {
            res.push_str(&format!("score={} ", score));
        };
        if let Some(pass) = self.score {
            res.push_str(&format!("pass={} ", pass));
        };
        for (key, val) in self.header.iter() {
            res.push_str(key);
            res.push('=');
            res.push_str(val);
            res.push(' ');
        }
        writeln!(f, "{}", res)?;
        for seq in &self.sequence {
            writeln!(f, "{}", seq)?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct Seq {
    name: String,
    start: u64,
    length: u64,
    strand: Strand,
    src_size: u64,
    text: Vec<u8>,
}

impl fmt::Display for Seq {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let text = String::from_utf8_lossy(&self.text).to_owned();
        write!(
            f,
            "s {} {} {} {} {} {}",
            self.name, self.start, self.length, self.strand, self.src_size, text
        )?;
        Ok(())
    }
}

impl Seq {
    fn clear(&mut self) {
        self.name.clear();
        self.text.clear();
    }
    fn update(&mut self, seq: Vec<&str>) -> Result<(), std::num::ParseIntError> {
        self.name.push_str(seq[1]);
        self.start = seq[2].parse()?;
        self.length = seq[3].parse()?;
        self.strand = if seq[4] == "+" {
            Strand::Forward
        } else if seq[4] == "-" {
            Strand::Reverse
        } else {
            unreachable!();
        };
        self.src_size = seq[5].parse()?;
        self.text.extend(seq[6].as_bytes());
        Ok(())
    }
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn start(&self) -> u64 {
        self.start
    }
    pub fn length(&self) -> u64 {
        self.length
    }
    pub fn strand(&self) -> Strand {
        self.strand
    }
    pub fn is_forward(&self) -> bool {
        match self.strand {
            Strand::Forward => true,
            Strand::Reverse => false,
        }
    }
    pub fn src_size(&self) -> u64 {
        self.src_size
    }
    pub fn text(&self) -> &[u8] {
        &self.text
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Strand::Forward => write!(f, "+")?,
            Strand::Reverse => write!(f, "-")?,
        };
        Ok(())
    }
}
impl Record {
    pub fn is_empty(&self) -> bool {
        self.score.is_none()
            && self.pass.is_none()
            && self.header.is_empty()
            && self.sequence.is_empty()
            && self.sequence_index == 0
    }
    pub fn clear(&mut self) {
        self.score = None;
        self.pass = None;
        self.header = Vec::new();
        self.sequence.iter_mut().for_each(|e| e.clear());
        self.sequence_index = 0;
    }
    pub fn score(&self) -> Option<f64> {
        self.score
    }
    pub fn pass(&self) -> Option<u64> {
        self.pass
    }
    pub fn other_header(&self) -> &[(String, String)] {
        &self.header
    }
    pub fn sequence(&self) -> &[Seq] {
        &self.sequence
    }
    /// Return the sequence with specified name, if any.
    /// If there are multiple sequence with that name,
    /// one of the maches would be returned.
    pub fn with_query_name(&self, name: &str) -> Option<&Seq> {
        self.sequence.iter().find(|&e| e.name() == name)
    }
    /// Return the first sequence mathces the pattern.
    /// Note that it accepet predicate, rather than just string.
    pub fn find_sequence<P>(&self, predicate: P) -> Option<&Seq>
    where
        P: FnMut(&&Seq) -> bool,
    {
        self.sequence.iter().find(predicate)
    }
    fn add_line(&mut self, line: &str) {
        // If the line is comment, ignore.
        if line.starts_with('a') {
            self.add_alignment(line);
        } else if line.starts_with('s') {
            self.add_sequence(line).unwrap();
        } else if !line.starts_with('#') {
            eprintln!("Currently we've not implemented the line such as {}", line);
        }
    }
    fn add_alignment(&mut self, line: &str) {
        let mut header: Vec<_> = vec![];
        for field in line.split_whitespace().skip(1) {
            let mut slots = field.split('=');
            let key = match slots.next() {
                Some(res) => res,
                None => continue,
            };
            let value = match slots.next() {
                Some(res) => res,
                None => continue,
            };
            if key == "score" {
                self.score = value.parse().ok();
            } else if key == "pass" {
                self.pass = value.parse().ok();
            } else {
                header.push((key.to_string(), value.to_string()));
            }
        }
        self.header = header;
    }
    fn add_sequence(&mut self, line: &str) -> Result<(), std::num::ParseIntError> {
        let seq: Vec<_> = line.split_whitespace().collect();
        if self.sequence_index < self.sequence().len() {
            self.sequence[self.sequence_index].update(seq)?;
        } else {
            let name = seq[1].to_string();
            let start: u64 = seq[2].parse()?;
            let length: u64 = seq[3].parse()?;
            let strand = if seq[4] == "+" {
                Strand::Forward
            } else if seq[4] == "-" {
                Strand::Reverse
            } else {
                unreachable!();
            };
            let src_size: u64 = seq[5].parse()?;
            let text: Vec<u8> = seq[6].as_bytes().to_vec();
            self.sequence.push(Seq {
                name,
                start,
                length,
                strand,
                src_size,
                text,
            });
        };
        self.sequence_index += 1;
        Ok(())
    }
}

#[derive(Debug)]
pub struct Records<R: io::Read> {
    reader: Reader<R>,
    has_error_occured: bool,
}

impl<R: io::Read> Records<R> {
    fn new(reader: Reader<R>) -> Self {
        Records {
            reader,
            has_error_occured: false,
        }
    }
}

impl<R: io::Read> Iterator for Records<R> {
    type Item = io::Result<Record>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.has_error_occured {
            None
        } else {
            let mut record = Record::default();
            match self.reader.read(&mut record) {
                Ok(_) if record.is_empty() => None,
                Ok(_) => Some(Ok(record)),
                Err(why) => {
                    self.has_error_occured = true;
                    Some(Err(why))
                }
            }
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_parser() {
        let file: Vec<_> = Reader::from_file("./testdata/test.maf")
            .unwrap()
            .records()
            .filter_map(|e| e.ok())
            .collect();
        assert_eq!(file.len(), 2);
    }

    #[test]
    fn test_file_content() {
        let file: Vec<_> = Reader::from_file("./testdata/test.maf")
            .unwrap()
            .records()
            .filter_map(|e| e.ok())
            .collect();
        assert!(!file[0].is_empty());
        assert_eq!(file[0].score(), Some(275.0));
        assert!(!file[0].other_header().is_empty());
        assert!(file[0].with_query_name("Chr11").is_some());
        assert!(file[0]
            .find_sequence(|e| e.name().starts_with("Ctg0"))
            .is_some());
        assert!(!file[1].is_empty());
        assert_eq!(file[1].score(), Some(251.0));
        assert!(!file[1].other_header().is_empty());
        assert!(file[1].with_query_name("Chr11").is_some());
        assert!(file[1].with_query_name("Ctg232").is_none());
    }
    #[test]
    fn test_file_content2() {
        let file: Vec<_> = Reader::from_file("./testdata/test.maf")
            .unwrap()
            .records()
            .filter_map(|e| e.ok())
            .collect();
        let seq1 = file[0].sequence();
        if seq1[0].name() == "Chr11" {
            assert_eq!(seq1[0].start(), 1122118);
            assert_eq!(seq1[0].length(), 385);
            assert_eq!(seq1[0].src_size(), 38115440);
        } else {
            debug_assert!(false, "{:?}", seq1[0]);
        }
        if seq1[1].name() == "Ctg0" {
            assert_eq!(seq1[1].start(), 26518);
            assert_eq!(seq1[1].length(), 385);
            assert_eq!(seq1[1].src_size(), 106115);
        } else {
            debug_assert!(false, "{:?}", seq1[0]);
        }
    }
    #[test]
    fn reuse_record() {
        // let mut answers: Vec<_> = Reader::from_file("./testdata/test.maf")
        //     .unwrap()
        //     .records()
        //     .filter_map(|e| e.ok())
        //     .collect();
        let mut res = vec![];
        let mut file = Reader::from_file("./testdata/test.maf").unwrap();
        let mut record = Record::default();
        let result = file.read(&mut record).unwrap();
        res.push(record.clone());
        eprintln!("{}", result);
        let result = file.read(&mut record).unwrap();
        res.push(record.clone());
        eprintln!("{}", result);
        let result = file.read(&mut record).unwrap();
        res.push(record.clone());
        eprintln!("{}", result);
        // debug_assert!(false,"{:?}",res);
    }
}
