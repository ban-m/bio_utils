fn main() -> std::io::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let start = std::time::Instant::now();
    let records = bio_utils::fastq::parse_into_vec(&args[1])?;
    let end = std::time::Instant::now();
    let time = (end - start).as_millis();
    println!("{}\t{}", records.len(), time);
    Ok(())
}
