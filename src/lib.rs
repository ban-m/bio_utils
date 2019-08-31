extern crate rand;
extern crate bio;
extern crate regex;
extern crate rayon;
// #[macro_use] extern crate lazy_static;
extern crate rust_htslib;
pub mod sam;
pub mod maf;
pub mod bam;
#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
