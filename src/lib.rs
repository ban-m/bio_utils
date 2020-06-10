extern crate bytecount;
extern crate rand;
extern crate rayon;
extern crate regex;
#[macro_use]
extern crate serde;
// pub mod bam;
pub mod fasta;
pub mod lasttab;
pub mod maf;
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
    seq.into_iter()
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
