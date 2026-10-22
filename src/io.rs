//! Reading and writing reads.

use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;

use anyhow::{bail, Context, Result};
use needletail::errors::ParseErrorKind;
use needletail::parse_fastx_reader;

/// One read. `seq` is upper case and holds only A, C, G, T and N.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Record {
    pub name: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

impl Record {
    pub fn new(name: &[u8], seq: &[u8], qual: Option<&[u8]>) -> Self {
        Record {
            name: name.to_vec(),
            seq: normalise(seq),
            qual: qual.map(|q| q.to_vec()),
        }
    }
}

/// Upper-cases a sequence and turns anything that is not A, C, G or T into N.
pub fn normalise(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|b| match b.to_ascii_uppercase() {
            c @ (b'A' | b'C' | b'G' | b'T') => c,
            _ => b'N',
        })
        .collect()
}

/// Reads every record of a FASTA or FASTQ file, plain or gzipped, or of stdin
/// when `path` is "-". A file with no records gives an empty vector rather than
/// an error, so "nothing to correct" is not mistaken for a broken input.
pub fn read_all(path: &Path) -> Result<Vec<Record>> {
    let source: Box<dyn Read + Send> = if path == Path::new("-") {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).with_context(|| format!("cannot open {}", path.display()))?)
    };
    read_from(source).with_context(|| format!("cannot parse {}", path.display()))
}

pub fn read_from<R: Read + Send + 'static>(source: R) -> Result<Vec<Record>> {
    let mut reader = match parse_fastx_reader(source) {
        Ok(reader) => reader,
        Err(e) if e.kind == ParseErrorKind::EmptyFile => return Ok(Vec::new()),
        Err(e) => bail!(e),
    };
    let mut reads = Vec::new();
    while let Some(record) = reader.next() {
        let record = record?;
        // the name is the first word of the header, as minimap2 and samtools take it
        let id = record.id();
        let name = id.split(|b| b.is_ascii_whitespace()).next().unwrap_or(id);
        reads.push(Record::new(name, &record.seq(), record.qual()));
    }
    Ok(reads)
}

pub enum Format {
    Fastq,
    Fasta,
}

pub fn open_output(path: &Path) -> Result<Box<dyn Write>> {
    Ok(if path == Path::new("-") {
        Box::new(BufWriter::new(io::stdout().lock()))
    } else {
        let file =
            File::create(path).with_context(|| format!("cannot create {}", path.display()))?;
        Box::new(BufWriter::new(file))
    })
}

pub fn write_read(out: &mut dyn Write, read: &Record, format: &Format) -> io::Result<()> {
    match format {
        Format::Fastq => {
            out.write_all(b"@")?;
            out.write_all(&read.name)?;
            out.write_all(b"\n")?;
            out.write_all(&read.seq)?;
            out.write_all(b"\n+\n")?;
            match &read.qual {
                Some(q) => out.write_all(q)?,
                None => out.write_all(&vec![b'#'; read.seq.len()])?,
            }
            out.write_all(b"\n")
        }
        Format::Fasta => {
            out.write_all(b">")?;
            out.write_all(&read.name)?;
            out.write_all(b"\n")?;
            out.write_all(&read.seq)?;
            out.write_all(b"\n")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn empty_input_is_no_reads() {
        assert!(read_from(Cursor::new(Vec::<u8>::new())).unwrap().is_empty());
    }

    #[test]
    fn parses_fastq_and_fasta() {
        let fq = b"@r1 extra words\nacgtRY\n+\nIIIIII\n".to_vec();
        let reads = read_from(Cursor::new(fq)).unwrap();
        assert_eq!(reads.len(), 1);
        assert_eq!(reads[0].name, b"r1");
        assert_eq!(reads[0].seq, b"ACGTNN");
        assert_eq!(reads[0].qual.as_deref(), Some(&b"IIIIII"[..]));

        let fa = b">a\nACG\nTTT\n>b\nGG\n".to_vec();
        let reads = read_from(Cursor::new(fa)).unwrap();
        assert_eq!(reads.len(), 2);
        assert_eq!(reads[0].seq, b"ACGTTT");
        assert_eq!(reads[1].qual, None);
    }

    #[test]
    fn garbage_is_an_error() {
        assert!(read_from(Cursor::new(b"this is not fastx\n".to_vec())).is_err());
    }
}
