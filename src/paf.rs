#[derive(Debug, Clone)]
pub struct PAF {
    pub qname: String,
    pub qlen: usize,
    pub qstart: usize,
    pub qend: usize,
    // True to forward
    pub relstrand: bool,
    pub tname: String,
    pub tlen: usize,
    pub tstart: usize,
    pub tend: usize,
    pub matchnum: usize,
    pub blocklen: usize,
    pub mapq: u16,
    pub tags: std::collections::HashMap<String, String>,
}


impl PAF {
    pub fn new(line: &str) -> Option<Self> {
        let mut line = line.split("\t");
        let res = Self {
            qname: line.next()?.to_string(),
            qlen: line.next()?.parse().ok()?,
            qstart: line.next()?.parse().ok()?,
            qend: line.next()?.parse().ok()?,
            relstrand: line.next()? == "+",
            tname: line.next()?.to_string(),
            tlen: line.next()?.parse().ok()?,
            tstart: line.next()?.parse().ok()?,
            tend: line.next()?.parse().ok()?,
            matchnum: line.next()?.parse().ok()?,
            blocklen: line.next()?.parse().ok()?,
            mapq: line.next()?.parse().ok()?,
            tags: line
                .filter_map(|tag| {
                    let mut tag = tag.split(":");
                    let key = tag.next()?.to_string();
                    tag.next()?;
                    let value = tag.next()?.to_string();
                    Some((key, value))
                })
                .collect(),
        };
        Some(res)
    }
}
