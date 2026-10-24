use std::io::Write;
use std::path::PathBuf;
use std::process::ExitCode;

use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};

use lrcorrect::io::{open_output, read_all, write_read, Format};
use lrcorrect::overlap::Preset;
use lrcorrect::{correct, subsample, Params};

#[derive(Clone, Copy, ValueEnum)]
enum ReadType {
    /// Oxford Nanopore (minimap2 ava-ont)
    Ont,
    /// PacBio CLR (minimap2 ava-pb)
    Pb,
}

/// Self-corrects a small set of long reads (plasmids, small genomes) without a
/// genome size, perl or java. Corrected reads are written as FASTQ (or FASTA).
/// No output reads, with exit status 0, means there was too little data to
/// correct anything; any other problem exits non-zero.
#[derive(Parser)]
#[command(version, about)]
struct Cli {
    /// Reads to correct: FASTQ or FASTA, optionally gzipped; "-" for stdin
    input: PathBuf,

    /// Output file; "-" for stdout
    #[arg(short, long, default_value = "-")]
    output: PathBuf,

    /// Threads
    #[arg(short, long, default_value_t = 1)]
    threads: usize,

    /// Read type
    #[arg(short, long, value_enum, default_value_t = ReadType::Ont)]
    preset: ReadType,

    /// Before correcting, randomly subsample the reads to about this many bases
    /// (0: use every read)
    #[arg(long, default_value_t = 0)]
    target_bases: u64,

    /// Seed for --target-bases
    #[arg(long, default_value_t = 11)]
    seed: u64,

    /// Rounds of correction
    #[arg(short, long, default_value_t = 2)]
    rounds: usize,

    /// Minimum number of other reads covering a base for it to be corrected
    #[arg(long, default_value_t = 3)]
    min_depth: u16,

    /// Minimum length of a corrected read (shorter pieces are dropped)
    #[arg(long, default_value_t = 500)]
    min_length: usize,

    /// Minimum overlap length between two reads
    #[arg(long, default_value_t = 500)]
    min_overlap: usize,

    /// Maximum number of overlapping reads used to correct each read
    /// [default: 60 for ont, 120 for pb]
    #[arg(long)]
    max_overlaps: Option<usize>,

    /// Largest unaligned stretch both reads may continue with beyond an
    /// overlap; above it, the match is taken to be a shared repeat
    #[arg(long, default_value_t = 1000)]
    max_overhang: usize,

    /// Minimum identity of an overlap alignment [default: 0.75 for ont, 0.65 for pb]
    #[arg(long)]
    min_identity: Option<f64>,

    /// Give every output base this quality character instead of one computed
    /// from the consensus vote
    #[arg(long)]
    qual_char: Option<char>,

    /// Cut reads where too few other reads cover them, instead of keeping those
    /// bases uncorrected. Unsupported read ends are trimmed either way
    #[arg(long)]
    split_unsupported: bool,

    /// Write FASTA instead of FASTQ
    #[arg(long)]
    fasta: bool,
}

fn run(cli: Cli) -> Result<()> {
    if cli.rounds == 0 {
        bail!("--rounds must be at least 1");
    }
    if cli.min_depth == 0 {
        bail!("--min-depth must be at least 1");
    }
    let fixed_qual = match cli.qual_char {
        Some(c) if c.is_ascii_graphic() => Some(c as u8),
        Some(c) => bail!("--qual-char must be a printable ASCII character, not {c:?}"),
        None => None,
    };
    let preset = match cli.preset {
        ReadType::Ont => Preset::Ont,
        ReadType::Pb => Preset::Pb,
    };
    let params = Params {
        preset,
        rounds: cli.rounds,
        min_depth: cli.min_depth,
        min_length: cli.min_length,
        min_overlap: cli.min_overlap,
        // noisier PacBio CLR reads need more of them for a stable consensus
        max_overlaps: cli.max_overlaps.unwrap_or(match preset {
            Preset::Ont => 60,
            Preset::Pb => 120,
        }),
        max_overhang: cli.max_overhang,
        min_identity: cli.min_identity.unwrap_or(match preset {
            Preset::Ont => 0.75,
            Preset::Pb => 0.65,
        }),
        fixed_qual,
        fill_holes: !cli.split_unsupported,
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads.max(1))
        .build_global()
        .context("cannot start the thread pool")?;

    let reads = read_all(&cli.input)?;
    let bases: usize = reads.iter().map(|r| r.seq.len()).sum();
    eprintln!(
        "lrcorrect: read {} reads ({} bp) from {}",
        reads.len(),
        bases,
        cli.input.display()
    );

    let lengths: Vec<usize> = reads.iter().map(|r| r.seq.len()).collect();
    let keep = subsample::choose(&lengths, cli.target_bases, cli.seed);
    let reads: Vec<_> = if keep.len() == reads.len() {
        reads
    } else {
        let mut keep = keep.into_iter().peekable();
        let picked: Vec<_> = reads
            .into_iter()
            .enumerate()
            .filter(|(i, _)| keep.next_if_eq(i).is_some())
            .map(|(_, r)| r)
            .collect();
        let kept: usize = picked.iter().map(|r| r.seq.len()).sum();
        eprintln!(
            "lrcorrect: subsampled to {} reads ({} bp)",
            picked.len(),
            kept
        );
        picked
    };

    let (corrected, stats) = if reads.len() < 2 {
        (Vec::new(), Vec::new())
    } else {
        correct(reads, &params)?
    };
    for (i, s) in stats.iter().enumerate() {
        eprintln!(
            "lrcorrect: round {}: {} reads ({} bp) in, {} with overlaps, {} corrected reads ({} bp) out; {:.1}s finding overlaps, {:.1}s correcting",
            i + 1,
            s.reads_in,
            s.bases_in,
            s.reads_with_overlaps,
            s.reads_out,
            s.bases_out,
            s.overlap_secs,
            s.consensus_secs
        );
    }

    let format = if cli.fasta {
        Format::Fasta
    } else {
        Format::Fastq
    };
    let mut out = open_output(&cli.output)?;
    for read in &corrected {
        write_read(&mut *out, read, &format)
            .with_context(|| format!("cannot write to {}", cli.output.display()))?;
    }
    out.flush()
        .with_context(|| format!("cannot write to {}", cli.output.display()))?;

    let out_bases: usize = corrected.iter().map(|r| r.seq.len()).sum();
    if corrected.is_empty() {
        eprintln!("lrcorrect: no reads could be corrected: too few reads overlap each other deeply enough");
    } else {
        eprintln!(
            "lrcorrect: wrote {} corrected reads ({} bp)",
            corrected.len(),
            out_bases
        );
    }
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("lrcorrect: error: {e:#}");
            ExitCode::FAILURE
        }
    }
}
