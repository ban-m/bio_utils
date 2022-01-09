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
    pub tags: Vec<(String, String, String)>,
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
                    let mut tag = tag.split(':').map(|x| x.to_string());
                    let key = tag.next()?;
                    let type_str = tag.next()?;
                    let value = tag.next()?;
                    Some((key, type_str, value))
                })
                .collect(),
        };
        Some(res)
    }
    pub fn get_tag(&self, key: &str) -> Option<(&str, &str)> {
        self.tags
            .iter()
            .find(|(k, _, _)| key == k)
            .map(|(_, tag_type, value)| (tag_type.as_str(), value.as_str()))
    }
}

impl std::fmt::Display for PAF {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let &Self {
            ref qname,
            qlen,
            qstart,
            qend,
            relstrand,
            ref tname,
            tlen,
            tstart,
            tend,
            matchnum,
            blocklen,
            mapq,
            ref tags,
        }: &PAF = self;
        let relstrand = if relstrand { '+' } else { '-' };
        let tags: String = tags
            .iter()
            .fold(String::new(), |mut tag, (key, val_type, value)| {
                if !tag.is_empty() {
                    tag.push('\t');
                }
                tag += key;
                tag.push(':');
                tag += val_type;
                tag.push(':');
                tag += value;
                tag
            });
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            qname,
            qlen,
            qstart,
            qend,
            relstrand,
            tname,
            tlen,
            tstart,
            tend,
            matchnum,
            blocklen,
            mapq,
            tags,
        )
    }
}
