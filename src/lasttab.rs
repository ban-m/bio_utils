//! LastTAB is a struct to represent an alignment record produced by `last` program.

/// The direction of the alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Hash)]
pub enum Strand {
    /// The strand is forward.
    Forward,
    /// The strand is reverse. So is the alignment coordinate.
    /// For example, if the length of sequence is 90, the start position is 10,
    /// the length of the alignment is 30,
    /// and the direction is `Strand::Reverse`, then the start position
    /// with respect to forward strand is 90 - 10 - 30 + 1.
    Reverse,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
        }
    }
}
impl Strand {
    pub fn is_forward(self) -> bool {
        match self {
            Self::Forward => true,
            Self::Reverse => false,
        }
    }
}

/// This is the information of alignment for a single strand.
/// Usually, a TAB-formatted alignment is losslessly represented by two `AlignInfo`
/// , an alignment pattern, and scores.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignInfo {
    /// The name of the sequence.
    seqname: String,
    /// The start position of the alignment. 0-based.
    /// If the direction is reverse, then the
    /// start position is the position at the reverse complement.
    seqstart: usize,
    /// The number of based mathed. In other words,
    /// the alignment region is [seqstart..seqstart+matchlen).
    matchlen: usize,
    /// The direction of the alignment.
    direction: Strand,
    /// The length of the sequence.
    seqlen: usize,
}

impl std::fmt::Display for AlignInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.seqname, self.seqstart, self.matchlen, self.direction, self.seqlen
        )
    }
}

impl AlignInfo {
    fn from_splits(splits: &[&str]) -> Option<Self> {
        let seqname = splits[0].to_string();
        let seqstart: usize = splits[1].parse().ok()?;
        let matchlen = splits[2].parse().ok()?;
        let direction = if splits[3] == "+" {
            Strand::Forward
        } else {
            Strand::Reverse
        };
        let seqlen = splits[4].parse().ok()?;
        Some(Self {
            seqname,
            seqstart,
            matchlen,
            direction,
            seqlen,
        })
    }
    fn seqstart_from_forward(&self) -> usize {
        match self.direction {
            Strand::Forward => self.seqstart,
            Strand::Reverse => self.seqlen - self.matchlen - self.seqstart,
        }
    }
}

use crate::sam::Record;
use std::collections::HashMap;
pub fn try_from(value: &Record, length: &HashMap<String, usize>) -> Result<LastTAB, &'static str> {
    if value.pos() == 0 {
        return Err("Alignment Invalid");
    }
    let score: i64 = value
        .attr()
        .iter()
        .find(|x| x.starts_with("AS"))
        .ok_or("AS tag is not valid.")?
        .split(':')
        .nth(2)
        .ok_or("AS tag does not have sufficient fields.")?
        .parse::<i64>()
        .map_err(|_| "Parsing fail.")?;
    let score = score.max(0) as u64;
    let eg2 = 0.;
    let e = 0.;
    let mut alignment = vec![];
    let (mut head_clip, mut _tail_clip) = (0, 0);
    let cigar = value.cigar();
    for op in cigar.iter() {
        use crate::sam::Op::*;
        match &op {
            Align(l) | Match(l) | Mismatch(l) => alignment.push(Op::Match(*l)),
            Insertion(l) => alignment.push(Op::Seq1In(*l)),
            Deletion(l) => alignment.push(Op::Seq2In(*l)),
            SoftClip(l) | HardClip(l) if head_clip == 0 => head_clip = *l,
            SoftClip(l) | HardClip(l) => _tail_clip = *l,
            Skipped(_) => return Err("Skipped in Cigar."),
            Padding(_) => return Err("Padding in Cigar."),
        }
    }
    let matchlen_1 = cigar
        .iter()
        .map(|op| {
            use crate::sam::Op::*;
            match &op {
                Deletion(l) | Align(l) | Match(l) | Mismatch(l) => *l,
                _ => 0,
            }
        })
        .sum::<usize>();
    // Reference.
    let seq1_len = length.get(value.r_name()).ok_or("Invalid Reference Name")?;
    let seq1_information = AlignInfo {
        seqname: value.r_name().to_string(),
        seqstart: value.pos() - 1,
        matchlen: matchlen_1,
        direction: Strand::Forward,
        seqlen: *seq1_len,
    };
    let matchlen_2 = cigar
        .iter()
        .map(|op| {
            use crate::sam::Op::*;
            match &op {
                Insertion(l) | Align(l) | Match(l) | Mismatch(l) => *l,
                _ => 0,
            }
        })
        .sum::<usize>();
    let total = cigar
        .iter()
        .map(|op| {
            use crate::sam::Op::*;
            match &op {
                HardClip(l) | SoftClip(l) | Insertion(l) | Align(l) | Match(l) | Mismatch(l) => *l,
                _ => 0,
            }
        })
        .sum::<_>();
    let direction = if value.is_forward() {
        Strand::Forward
    } else {
        Strand::Reverse
    };
    let seq2_information = AlignInfo {
        seqname: value.q_name().to_string(),
        seqstart: head_clip,
        matchlen: matchlen_2,
        direction,
        seqlen: total,
    };
    let lt = LastTAB {
        score,
        alignment,
        eg2,
        e,
        seq1_information,
        seq2_information,
    };
    Ok(lt)
}

/// A struct to represent a last's TAB-format alignment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LastTAB {
    seq1_information: AlignInfo,
    seq2_information: AlignInfo,
    score: u64,
    alignment: Vec<Op>,
    eg2: f64,
    e: f64,
}

impl std::fmt::Display for LastTAB {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let alignment: Vec<_> = self.alignment.iter().map(|x| format!("{}", x)).collect();
        write!(
            f,
            "{}\t{}\t{}\t{}\tEG2={}\tE={}",
            self.score,
            self.seq1_information,
            self.seq2_information,
            alignment.join(","),
            self.eg2,
            self.e
        )
    }
}

impl PartialEq for LastTAB {
    fn eq(&self, other: &Self) -> bool {
        let seq1 = self.seq1_name() == other.seq1_name()
            && self.seq1_start_from_forward() == other.seq1_start_from_forward()
            && self.seq1_end_from_forward() == other.seq1_end_from_forward();
        let seq2 = self.seq2_name() == other.seq2_name()
            && self.seq2_start_from_forward() == other.seq2_start_from_forward()
            && self.seq2_end_from_forward() == other.seq2_end_from_forward();
        let temp = seq1 && seq2;
        let seq1 = self.seq1_name() == other.seq2_name()
            && self.seq1_start_from_forward() == other.seq2_start_from_forward()
            && self.seq1_end_from_forward() == other.seq2_end_from_forward();
        let seq2 = self.seq2_name() == other.seq1_name()
            && self.seq2_start_from_forward() == other.seq1_start_from_forward()
            && self.seq2_end_from_forward() == other.seq1_end_from_forward();
        let rev = seq1 && seq2;
        temp || rev
    }
}
impl Eq for LastTAB {}

impl LastTAB {
    pub fn from_line(line: &str) -> Option<Self> {
        let line: Vec<&str> = line.split('\t').collect();
        let score: u64 = line[0].parse().ok()?;
        let seq1_information = AlignInfo::from_splits(&line[1..=5])?;
        let seq2_information = AlignInfo::from_splits(&line[6..=10])?;
        let alignment = line[11].split(',').fold(vec![], |mut res, op| {
            Op::from_string(&mut res, op);
            res
        });
        let (mut eg2, mut e) = (2., 3.); // Dummy values
        if line.len() > 13 {
            if line[12].starts_with("E=") {
                e = match line[12][2..].parse() {
                    Ok(res) => res,
                    Err(why) => panic!("{},{}", why, &line[12][2..]),
                };
            } else if line[12].starts_with("EG2=") {
                eg2 = match line[12][4..].parse() {
                    Ok(res) => res,
                    Err(why) => panic!("{},{}", why, &line[12][4..]),
                };
            };
            if line[13].starts_with("E=") {
                e = match line[13][2..].parse() {
                    Ok(res) => res,
                    Err(why) => panic!("{},{}", why, &line[13][2..]),
                };
            } else if line[13].starts_with("EG2=") {
                eg2 = match line[13][4..].parse() {
                    Ok(res) => res,
                    Err(why) => panic!("{},{}", why, &line[13][4..]),
                };
            };
        }
        Some(Self {
            score,
            seq1_information,
            seq2_information,
            alignment,
            e,
            eg2,
        })
    }
    pub fn score(&self) -> u64 {
        self.score
    }
    pub fn seq1_name(&self) -> &str {
        &self.seq1_information.seqname
    }
    pub fn seq2_name(&self) -> &str {
        &self.seq2_information.seqname
    }
    pub fn seq1_start(&self) -> usize {
        self.seq1_information.seqstart
    }
    /// The location where the alignment start,
    /// counted from the start position of the sequence, regardless of the strand.
    /// Thus, if the strand is reversed, the actual alignment starts from seq[start+len] and
    /// end at seq[start], in rev cmp manner.
    pub fn seq1_start_from_forward(&self) -> usize {
        self.seq1_information.seqstart_from_forward()
    }
    pub fn seq2_start(&self) -> usize {
        self.seq2_information.seqstart
    }
    pub fn seq2_start_from_forward(&self) -> usize {
        self.seq2_information.seqstart_from_forward()
    }
    pub fn seq1_matchlen(&self) -> usize {
        self.seq1_information.matchlen
    }
    pub fn seq2_matchlen(&self) -> usize {
        self.seq2_information.matchlen
    }
    pub fn seq1_end_from_forward(&self) -> usize {
        self.seq1_start_from_forward() + self.seq1_matchlen()
    }
    pub fn seq2_end_from_forward(&self) -> usize {
        self.seq2_start_from_forward() + self.seq2_matchlen()
    }
    pub fn seq1_direction(&self) -> Strand {
        self.seq1_information.direction
    }
    pub fn seq2_direction(&self) -> Strand {
        self.seq2_information.direction
    }
    pub fn seq1_len(&self) -> usize {
        self.seq1_information.seqlen
    }
    pub fn seq2_len(&self) -> usize {
        self.seq2_information.seqlen
    }
    pub fn alignment(&self) -> &[Op] {
        &self.alignment
    }
    pub fn e_score(&self) -> f64 {
        self.e
    }
    pub fn eg2_score(&self) -> f64 {
        self.eg2
    }
    // Return alignment length. Not the length of the reference nor the query.
    pub fn alignment_length(&self) -> usize {
        self.alignment
            .iter()
            .map(|op| match &op {
                Op::Match(l) => l,
                Op::Seq1In(l) => l,
                Op::Seq2In(l) => l,
            })
            .sum()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Hash)]
pub enum Op {
    /// Match of `usize` length.
    Match(usize),
    /// Sequence insertion of `usize` length to the seq1.
    /// In other words, if we encounter Seq1In(l), the location of the
    /// sequence 2 would increase by l.
    Seq1In(usize),
    /// Sequence insertion with `usize` length to the seq2.
    /// In other words, if we encounter Seq1In(l), the location of the
    /// sequence 1 would increase by l.
    Seq2In(usize),
}
impl std::fmt::Display for Op {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Op::Match(l) => write!(f, "{}", l),
            Op::Seq1In(l) => write!(f, "0:{}", l),
            Op::Seq2In(l) => write!(f, "{}:0", l),
        }
    }
}

impl Op {
    fn from_string(res: &mut Vec<Op>, input: &str) {
        if input.contains(':') {
            let mut input = input.split(':');
            let seq1 = input.next().and_then(|e| e.parse().ok()).unwrap();
            let seq2 = input.next().and_then(|e| e.parse().ok()).unwrap();
            if seq1 != 0 {
                res.push(Op::Seq2In(seq1));
            }
            if seq2 != 0 {
                res.push(Op::Seq1In(seq2));
            }
        } else {
            res.push(Op::Match(input.parse().unwrap()))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const LAST_INPUT:&str ="1035\ttig00000001\t98045\t539\t+\t261026\tm54113_160913_184949/5570667/0_1125\t4\t527\t-\t1125\t10,1:0,5,1:0,3,1:0,5,0:1,3\tEG2=8.2e-86\tE=6.6e-95";
    #[test]
    fn last_parse_test() {
        eprintln!("Start");
        let aln = LastTAB::from_line(LAST_INPUT).unwrap();
        assert_eq!(aln.score(), 1035);
        assert_eq!(aln.seq1_name(), "tig00000001");
        assert_eq!(aln.seq1_start(), 98045);
        assert_eq!(aln.seq1_matchlen(), 539);
        assert_eq!(aln.seq1_direction(), Strand::Forward);
        assert_eq!(aln.seq1_len(), 261026);
        assert_eq!(aln.seq2_name(), "m54113_160913_184949/5570667/0_1125");
        assert_eq!(aln.seq2_start(), 4);
        assert_eq!(aln.seq2_matchlen(), 527);
        assert_eq!(aln.seq2_direction(), Strand::Reverse);
        assert_eq!(aln.seq2_len(), 1125);
        use Op::*;
        assert_eq!(
            aln.alignment(),
            vec![
                Match(10),
                Seq2In(1),
                Match(5),
                Seq2In(1),
                Match(3),
                Seq2In(1),
                Match(5),
                Seq1In(1),
                Match(3)
            ]
        );
        assert_eq!(aln.e_score(), 6.6e-95);
        assert_eq!(aln.eg2_score(), 8.2e-86);
        assert_eq!(aln.seq1_start_from_forward(), 98045);
        assert_eq!(aln.seq1_end_from_forward(), 98045 + 539);
        assert_eq!(aln.seq2_start_from_forward(), 1125 - 527 - 4);
        assert_eq!(aln.seq2_end_from_forward(), 1125 - 4);
    }
}
