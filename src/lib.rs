#[macro_use]
extern crate serde;
pub mod alignments;
pub mod fasta;
pub mod fastq;
pub mod lasttab;
pub mod maf;
pub mod paf;
pub mod sam;
#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}

#[inline]
pub fn revcmp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&e| match e {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => unreachable!(),
        })
        .collect()
}
