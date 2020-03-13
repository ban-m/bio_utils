extern crate bio;
extern crate bytecount;
extern crate rand;
extern crate rayon;
extern crate regex;
#[macro_use]
extern crate serde;
// #[macro_use] extern crate lazy_static;
extern crate rust_htslib;
pub mod bam;
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
