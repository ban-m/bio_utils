/// Reconstruct bam alignmnet and output pritty strings.
/// The seq1 should be the query, while seq2 is the reference.
/// The first position where the reference consumed should be `pos`.
/// The return value consists of three vectors of characters,
/// first for the query, second for operations, and the third for
/// reference. As to the user can output digit, the output value is vectors,
/// instead of `String`s.
/// Note that the sequences should be 'revcomp'ed if the alignment is revcomp.
pub fn recover_alignment<'a>(
    iter: &'a [rust_htslib::bam::record::Cigar],
    seq1: &[u8],
    seq2: &[u8],
    pos: usize,
) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let empty_string = |len| (0..len).map(|_| " ").collect::<String>();
    use rust_htslib::bam::record::Cigar::*;
    let (mut seq1_with_gap, mut seq2_with_gap, mut operations) = (vec![], vec![], vec![]);
    let (mut seq1idx, mut seq2idx) = (0, pos );
    let seq1_header = match iter[0] {
        SoftClip(l) | HardClip(l) => {
            seq1idx += l as usize;
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
            Match(l) => {
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
            Del(l) => {
                let l = *l as usize;
                seq1_with_gap.extend(vec![b'-'; l]);
                seq2_with_gap.extend_from_slice(&seq2[seq2idx..(seq2idx + l)]);
                operations.extend(vec![b' '; l]);
                seq2idx += l;
            }
            Ins(l) => {
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
