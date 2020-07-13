use rayon::prelude::*;
// use regex::Regex;
use std::cmp::max;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::slice;

pub fn into_coverage(reads: &[Sam]) -> Vec<Coverage> {
    reads
        .par_iter()
        .fold(Vec::new, |acc, sam| combine_sam(acc, sam))
        .reduce(Vec::new, merge_two_coverages)
}

#[derive(Debug, Clone)]
pub struct Coverage {
    r_name: String,
    cov: Vec<(usize, u64)>,
}

fn merge_two_coverages(a: Vec<Coverage>, b: Vec<Coverage>) -> Vec<Coverage> {
    let mut res = Vec::with_capacity(a.len().max(b.len()));
    let mut aiter = a.into_iter().peekable();
    let mut biter = b.into_iter().peekable();
    while let (Some(arname), Some(brname)) = (
        aiter.peek().map(|e| e.r_name.clone()),
        biter.peek().map(|e| e.r_name.clone()),
    ) {
        match arname.cmp(&brname) {
            std::cmp::Ordering::Less => res.push(aiter.next().unwrap().clone()),
            std::cmp::Ordering::Greater => res.push(biter.next().unwrap().clone()),
            std::cmp::Ordering::Equal => {
                let acov = aiter.next().unwrap();
                let bcov = biter.next().unwrap();
                res.push(acov.merge(&bcov).clone());
            }
        }
    }
    res.extend(aiter);
    res.extend(biter);
    res
}

impl Coverage {
    pub fn r_name(&self) -> &str {
        &self.r_name
    }
    pub fn get_cov(&self) -> slice::Iter<(usize, u64)> {
        self.cov.iter()
    }

    fn merge(&self, cov: &Self) -> Self {
        let mut res = Vec::with_capacity(max(self.cov.len(), cov.cov.len()));
        if self.r_name != cov.r_name {
            panic!("merging {:?},{:?}", self, cov);
        } else {
            let mut selfiter = self.cov.iter().peekable();
            let mut coviter = cov.cov.iter().peekable();
            while let (Some((s_index, _)), Some((c_index, _))) = (selfiter.peek(), coviter.peek()) {
                match s_index.cmp(&c_index) {
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

fn combine_sam(mut acc: Vec<Coverage>, sam: &Sam) -> Vec<Coverage> {
    match acc.binary_search_by(|probe| probe.r_name.cmp(&sam.r_name)) {
        Err(index) => acc.insert(index, sam.to_coverage()),
        Ok(index) => acc[index] = acc[index].merge(&sam.to_coverage()),
    }
    acc
}

pub fn load_sam_file<P: AsRef<Path>>(file: P) -> io::Result<Vec<Sam>> {
    // annotation for retuened value:(b,cigar,pos) where,
    //b: true when the read mapped to template strand,
    // cigar: cigar string for the alignment,
    // pos: First position of matched base
    Ok(BufReader::new(File::open(file)?)
        .lines()
        .filter_map(|e| e.ok())
        .skip_while(|e| e.starts_with('@'))
        .filter_map(|e| Sam::new(&e))
        .collect())
}

#[derive(Debug, Clone)]
pub struct Sam {
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
    attr: Option<String>,
}
use std::fmt;
impl fmt::Display for Sam {
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
        match &self.attr {
            Some(ref attr) => write!(f, "\t{}", attr),
            None => write!(f, ""),
        }
    }
}

impl Sam {
    pub fn new(input: &str) -> Option<Self> {
        let mut contents = input.split('\t');
        let q_name = contents.next()?.to_string();
        let flag = contents.next()?.parse().ok()?;
        let r_name = contents.next()?.to_string();
        let pos = contents.next()?.parse().ok()?;
        let mapq = contents.next()?.parse().ok()?;
        let cigar = contents.next()?.to_string();
        let rnext = contents.next()?.to_string();
        let pnext = contents.next()?.parse().ok()?;
        let tlen = contents.next()?.parse().ok()?;
        let seq = contents.next()?.to_string();
        let qual = contents.next()?.bytes().map(|e| e - 33).collect();
        let attr = contents.next().map(|e| e.to_string());
        Some(Sam {
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
    pub fn pos(&self) -> usize {
        self.pos
    }
    /// Return the mapping region with respect to the query(0-based).
    /// If wanted to get the range with respect to reference, use `get_range` instead.
    pub fn mapped_region(&self) -> (usize, usize) {
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
    /// If wanted to get the range with respect to the query, use `mapped_region` instead.
    pub fn get_range(&self) -> (usize, usize) {
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
    #[inline]
    pub fn cigar(&self) -> Vec<Op> {
        parse_cigar_string(&self.cigar)
    }
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
    fn cigar_as_str(&self) -> &str {
        &self.cigar
    }
    fn qual_as_str(&self) -> String {
        self.qual.iter().map(|e| (e + 33) as char).collect()
    }
}

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
#[inline]
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

/// Reconstruct bam alignmnet and output pritty strings.
/// The seq1 should be the query, while seq2 is the reference.
/// The first position where the reference consumed should be `pos`.
/// The return value consists of three vectors of characters,
/// first for the query, second for operations, and the third for
/// reference. As to the user can output digit, the output value is vectors,
/// instead of `String`s.
/// Note that the sequences should be 'revcomp'ed if the alignment is revcomp.
pub fn recover_alignment(
    iter: &[Op],
    seq1: &[u8],
    seq2: &[u8],
    pos: usize,
) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let empty_string = |len| (0..len).map(|_| " ").collect::<String>();
    use Op::*;
    let (mut seq1_with_gap, mut seq2_with_gap, mut operations) = (vec![], vec![], vec![]);
    let (mut seq1idx, mut seq2idx) = (0, pos - 1);
    let seq1_header = match iter[0] {
        SoftClip(l) | HardClip(l) => {
            seq1idx += l;
            format!("[head {:05} base]", l)
        }
        _ => "[head 00000 base]".to_string(),
    };
    let seq2_header = format!("[head {:05} base]", pos);
    let ops_header = empty_string("[head 00000 base]".len());
    seq1_with_gap.extend(seq1_header.as_bytes());
    seq2_with_gap.extend(seq2_header.as_bytes());
    operations.extend(ops_header.as_bytes());
    for op in iter {
        match op {
            Align(l) => {
                let l = *l as usize;
                seq1_with_gap.extend_from_slice(&seq1[seq1idx..(seq1idx + l)]);
                seq2_with_gap.extend_from_slice(&seq2[seq2idx..(seq2idx + l)]);
                operations.extend(match_mismatch(
                    &seq1[seq1idx..(seq1idx + l)],
                    &seq2[seq2idx..(seq2idx + l)],
                ));
                seq1idx += l;
                seq2idx += l;
            }
            Deletion(l) => {
                let l = *l as usize;
                seq1_with_gap.extend(vec![b'-'; l]);
                seq2_with_gap.extend_from_slice(&seq2[seq2idx..(seq2idx + l)]);
                operations.extend(vec![b' '; l]);
                seq2idx += l;
            }
            Insertion(l) => {
                let l = *l as usize;
                seq1_with_gap.extend_from_slice(&seq1[seq1idx..(seq1idx + l)]);
                seq2_with_gap.extend(vec![b'-'; l]);
                operations.extend(vec![b' '; l]);
                seq1idx += l;
            }
            _ => {}
        }
    }
    let seq1_footer = match iter.last().unwrap() {
        SoftClip(l) | HardClip(l) => format!("[tail {:05} bese]", l),
        _ => "[tail 00000 base]".to_string(),
    };
    let seq2_footer = format!("[tail {:05} base]", seq2.len() - seq2idx);
    let ops_footer = empty_string("[tail 00000 base]".len());
    seq1_with_gap.extend(seq1_footer.as_bytes());
    seq2_with_gap.extend(seq2_footer.as_bytes());
    operations.extend(ops_footer.as_bytes());
    (seq1_with_gap, operations, seq2_with_gap)
}

fn match_mismatch(xs: &[u8], ys: &[u8]) -> Vec<u8> {
    xs.iter()
        .zip(ys.iter())
        .map(|(x, y)| if x == y { b'|' } else { b'X' })
        .collect()
}

#[test]
fn cigar_parse() {
    use super::sam::Op::*;
    let cigar = "101S33M2I66M";
    let processed = parse_cigar_string(&cigar);
    eprintln!("{:?}", processed);
    assert_eq!(
        processed,
        vec![SoftClip(101), Align(33), Insertion(2), Align(66)]
    );
}
