//! Tiny library to read SAM file(read only).
use std::cmp::max;
use std::io::BufRead;

/// Coverage on a contig.
#[derive(Debug, Clone)]
pub struct Coverage {
    r_name: String,
    cov: Vec<(usize, u64)>,
}

fn combine_sam(mut acc: Vec<Coverage>, sam: &Record) -> Vec<Coverage> {
    match acc.binary_search_by(|probe| probe.r_name.cmp(&sam.r_name)) {
        Err(index) => acc.insert(index, sam.to_coverage()),
        Ok(index) => acc[index] = acc[index].merge(&sam.to_coverage()),
    }
    acc
}

impl Coverage {
    /// Name of the reference (or contig)
    pub fn r_name(&self) -> &str {
        &self.r_name
    }
    /// Get the coverage.
    pub fn cov(&self) -> &[(usize, u64)] {
        self.cov.as_slice()
    }
    /// Convert SAM records into the coverages.
    pub fn new(records: &[Record]) -> Vec<Coverage> {
        records.iter().fold(vec![], combine_sam)
    }
    fn merge(&self, cov: &Self) -> Self {
        let mut res = Vec::with_capacity(max(self.cov.len(), cov.cov.len()));
        if self.r_name != cov.r_name {
            panic!("merging {:?},{:?}", self, cov);
        } else {
            let mut selfiter = self.cov.iter().peekable();
            let mut coviter = cov.cov.iter().peekable();
            while let (Some((s_index, _)), Some((c_index, _))) = (selfiter.peek(), coviter.peek()) {
                match s_index.cmp(c_index) {
                    std::cmp::Ordering::Less => res.push(*selfiter.next().unwrap()),
                    std::cmp::Ordering::Greater => res.push(*coviter.next().unwrap()),
                    std::cmp::Ordering::Equal => {
                        let (s_index, s_depth) = selfiter.next().unwrap();
                        let (_c_index, c_depth) = coviter.next().unwrap();
                        res.push((*s_index, s_depth + c_depth))
                    }
                }
            }
            res.extend(selfiter);
            res.extend(coviter);
        }
        Coverage {
            r_name: self.r_name.clone(),
            cov: res,
        }
    }
}

/// SAM file.
#[derive(Debug, Clone)]
pub struct Sam {
    /// Headers(begin with '@')
    pub headers: Vec<Header>,
    /// SAM records.
    pub records: Vec<Record>,
}

impl Sam {
    /// Read SAM file from the BufferedReader.
    pub fn from_reader<R: BufRead>(rdr: R) -> Sam {
        let mut headers = vec![];
        let mut records = vec![];
        for line in rdr.lines().filter_map(|x| x.ok()) {
            if line.starts_with('@') {
                headers.push(Header::new(&line).unwrap());
            } else {
                records.push(line.parse::<Record>().unwrap());
            }
        }
        Self { headers, records }
    }
}

/// SAM header file.
#[derive(Debug, Clone)]
pub struct Header {
    /// Tag name (NN for `@NN`)
    pub tag: String,
    /// Attributes for this tag. Each attribute is separated by '\t',
    /// and in "AttributeName:AttributeValue" format.
    pub attrs: Vec<(String, String)>,
}

impl Header {
    fn new(line: &str) -> Option<Self> {
        let mut line = line.split('\t');
        let tag: String = line.next()?.trim_start_matches('@').to_string();
        let attrs: Vec<_> = line
            .filter_map(|attr| {
                let mut attr = attr.splitn(2, ':');
                let key = attr.next()?.to_string();
                let value = attr.next()?.to_string();
                Some((key, value))
            })
            .collect();
        Some(Self { tag, attrs })
    }
}

/// SAM Record. The files can be accessed via method calling, such as [`Record::q_name()`].
/// Since this struct implements [`std::str::FromStr`], it is possible to `let sam_record:Sam = line.parse().unwrap();` to parse the record.
#[derive(Debug, Clone)]
pub struct Record {
    q_name: String,
    flag: u32,
    r_name: String,
    pos: usize,
    mapq: usize,
    cigar: String,
    rnext: String,
    pnext: usize,
    tlen: usize,
    seq: String,
    qual: Vec<u8>,
    attr: Vec<String>,
}

use std::fmt;
impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.q_name,
            self.flag,
            self.r_name,
            self.pos,
            self.mapq,
            self.cigar_as_str(),
            self.rnext,
            self.pnext,
            self.tlen,
            self.seq,
            self.qual_as_str()
        )?;
        if !self.attr.is_empty() {
            write!(f, "\t{}", self.attr.join("\t"))?;
        }
        Ok(())
    }
}

/// Error raised when parsing of a SAM record fails
#[derive(Debug, Clone, Copy)]
pub struct ParseSamError {}

impl std::str::FromStr for Record {
    type Err = ParseSamError;
    fn from_str(input: &str) -> Result<Self, Self::Err> {
        let mut contents = input.split('\t');
        fn map<T>(x: Option<T>) -> Result<T, ParseSamError> {
            match x {
                Some(x) => Ok(x),
                None => Err(ParseSamError {}),
            }
        }
        let q_name = map(contents.next().map(|s| s.to_string()))?;
        let flag = map(contents.next().and_then(|x| x.parse().ok()))?;
        let r_name = map(contents.next())?.to_string();
        let pos = map(contents.next().and_then(|x| x.parse().ok()))?;
        let mapq = map(contents.next().and_then(|x| x.parse().ok()))?;
        let cigar = map(contents.next())?.to_string();
        let rnext = map(contents.next())?.to_string();
        let pnext = map(contents.next().and_then(|x| x.parse().ok()))?;
        let tlen = map(contents.next().and_then(|x| x.parse().ok()))?;
        let seq = map(contents.next())?.to_string();
        let qual = map(contents.next())?.bytes().map(|e| e - 33).collect();
        let attr = map(contents.next().map(|e| e.to_string()))?;
        let attr: Vec<_> = attr.split('\t').map(|x| x.to_string()).collect();
        Ok(Self {
            q_name,
            flag,
            r_name,
            pos,
            mapq,
            cigar,
            rnext,
            pnext,
            tlen,
            seq,
            qual,
            attr,
        })
    }
}

impl Record {
    pub fn to_coverage(&self) -> Coverage {
        let mut cov = vec![];
        let mut start = self.pos; // reference position
        for op in &self.cigar() {
            use self::Op::*;
            match *op {
                Align(b) | Match(b) => {
                    for i in 0..b {
                        cov.push((start + i, 1));
                    }
                    start += b;
                }
                Insertion(_) | SoftClip(_) | HardClip(_) | Padding(_) => {}
                Deletion(b) | Skipped(b) | Mismatch(b) => start += b,
            }
        }
        Coverage {
            r_name: self.r_name.to_string(),
            cov,
        }
    }
    pub fn q_name(&self) -> &str {
        &self.q_name
    }
    pub fn r_name(&self) -> &str {
        &self.r_name
    }
    pub fn ref_name(&self) -> &str {
        &self.r_name
    }
    pub fn mapq(&self) -> usize {
        self.mapq
    }
    pub fn is_primary(&self) -> bool {
        (self.flag & 0x900) == 0
    }
    pub fn is_template(&self) -> bool {
        (self.flag & 0b10000) != 0b10000
    }
    pub fn is_forward(&self) -> bool {
        (self.flag & 0b10000) != 0b10000
    }
    pub fn flag(&self) -> u32 {
        self.flag
    }
    pub fn seq(&self) -> &str {
        &self.seq
    }
    /// Mapped position (1-based)
    pub fn pos(&self) -> usize {
        self.pos
    }
    /// Return the mapping region with respect to the query (0-based).
    /// If wanted to get the range with respect to reference, use `Self::refr_aligned_region` instead.
    pub fn query_aligned_region(&self) -> (usize, usize) {
        use self::Op::*; // 0-BASED!!!!!
        let (head_clip, middle, _tail_clip, _) =
            self.cigar().iter().fold((0, 0, 0, true), |acc, x| match x {
                HardClip(b) | SoftClip(b) if acc.3 => (acc.0 + b, acc.1, acc.2, acc.3),
                HardClip(b) | SoftClip(b) if !acc.3 => (acc.0, acc.1, acc.2 + b, acc.3),
                Align(b) | Insertion(b) | Match(b) | Mismatch(b) => {
                    (acc.0, acc.1 + b, acc.2, false)
                }
                _ => acc,
            });
        (head_clip, head_clip + middle)
    }
    /// Return the mapping region with respect to the reference(0-based).
    /// If wanted to get the range with respect to the query, use `Self::query_aligned_region` instead.
    pub fn refr_aligned_region(&self) -> (usize, usize) {
        // Return the position of the genome(measured in template). 0-BASED!!!!
        let start = self.pos;
        if start == 0 {
            return (0, 0);
        };
        use self::Op::*;
        let len: usize = self
            .cigar()
            .iter()
            .map(|op| match *op {
                Align(b) | Match(b) | Deletion(b) | Skipped(b) | Mismatch(b) => b,
                Insertion(_) | SoftClip(_) | HardClip(_) | Padding(_) => 0,
            })
            .sum();
        (start - 1, start + len - 1)
    }
    /// Return the length of the query.
    pub fn query_length(&self) -> usize {
        self.cigar()
            .iter()
            .map(|e| match e {
                Op::HardClip(b)
                | Op::SoftClip(b)
                | Op::Align(b)
                | Op::Match(b)
                | Op::Mismatch(b)
                | Op::Insertion(b) => *b,
                _ => 0,
            })
            .sum()
    }
    /// Parse and return the Cigar string.
    /// This method takes `O(|L|)`-time, where `L` is the length of the Cigar string.
    pub fn cigar(&self) -> Vec<Op> {
        parse_cigar_string(&self.cigar)
    }
    fn cigar_as_str(&self) -> &str {
        &self.cigar
    }
    fn qual_as_str(&self) -> String {
        self.qual.iter().map(|e| (e + 33) as char).collect()
    }
    /// Attributes of this record. Usually, each element is formatted as "[TAG_NAME]:[TAG_TYPE]:[TAG_VALUE]".
    pub fn attr(&self) -> &[String] {
        self.attr.as_slice()
    }
}

/// Alignment operations. Insertions are insertions to the reference, and deletions are deletions from the reference.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum Op {
    Align(usize),     //M
    Insertion(usize), //I
    Deletion(usize),  //D
    Skipped(usize),   //N
    SoftClip(usize),  //S
    HardClip(usize),  //H
    Padding(usize),   //P
    Match(usize),     //=
    Mismatch(usize),  //X
}

impl Op {
    pub fn new(op: &str) -> Option<Op> {
        let (operation, num): (char, usize) = {
            let mut op = String::from(op);
            let operation = op.pop()?;
            (operation, op.parse().ok()?)
        };
        match operation {
            'M' => Some(Op::Align(num)),
            'I' => Some(Op::Insertion(num)),
            'D' => Some(Op::Deletion(num)),
            'N' => Some(Op::Skipped(num)),
            'S' => Some(Op::SoftClip(num)),
            'H' => Some(Op::HardClip(num)),
            'P' => Some(Op::Padding(num)),
            '=' => Some(Op::Match(num)),
            'X' => Some(Op::Mismatch(num)),
            _ => None,
        }
    }
    pub fn from(num: usize, op: u8) -> Option<Op> {
        match op {
            b'M' => Some(Op::Align(num)),
            b'I' => Some(Op::Insertion(num)),
            b'D' => Some(Op::Deletion(num)),
            b'N' => Some(Op::Skipped(num)),
            b'S' => Some(Op::SoftClip(num)),
            b'H' => Some(Op::HardClip(num)),
            b'P' => Some(Op::Padding(num)),
            b'=' => Some(Op::Match(num)),
            b'X' => Some(Op::Mismatch(num)),
            _ => None,
        }
    }
    pub fn as_str(&self) -> String {
        let (num, op) = match self {
            Op::Align(x) => (x, 'M'),
            Op::Insertion(x) => (x, 'I'),
            Op::Deletion(x) => (x, 'D'),
            Op::Skipped(x) => (x, 'N'),
            Op::SoftClip(x) => (x, 'S'),
            Op::HardClip(x) => (x, 'H'),
            Op::Padding(x) => (x, 'P'),
            Op::Match(x) => (x, '='),
            Op::Mismatch(x) => (x, 'X'),
        };
        format!("{}{}", num, op)
    }
}

/// Parse a given CIGAR string. If it is not a valid CIGAR, panic.
pub fn parse_cigar_string(cigar: &str) -> Vec<Op> {
    let mut ops = vec![];
    let mut num = 0;
    for x in cigar.bytes() {
        if x.is_ascii_digit() {
            num = 10 * num + (x - b'0') as usize;
        } else {
            if let Some(res) = Op::from(num, x) {
                ops.push(res);
            }
            num = 0;
        }
    }
    ops
}

#[test]
fn cigar_parse() {
    use super::sam::Op::*;
    let cigar = "101S33M2I66M";
    let processed = parse_cigar_string(cigar);
    eprintln!("{:?}", processed);
    assert_eq!(
        processed,
        vec![SoftClip(101), Align(33), Insertion(2), Align(66)]
    );
}
