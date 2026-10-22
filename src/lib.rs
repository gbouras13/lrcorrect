//! lrcorrect: self-correction of small long-read sets.
//!
//! Reads are optionally subsampled to a target number of bases, then corrected
//! in one or more rounds. Each round finds overlaps between all reads with
//! minimap2, aligns every read's best partners to it, and calls a pileup
//! consensus. Bases too few reads support are cut away, so a corrected read can
//! come out shorter, in pieces, or not at all.

pub mod consensus;
pub mod io;
pub mod overlap;
pub mod subsample;

use rayon::prelude::*;

use crate::consensus::{Aligned, Pileup};
use crate::io::Record;
use crate::overlap::{find_overlaps, Preset, TargetAligner};

pub struct Params {
    pub preset: Preset,
    pub rounds: usize,
    pub min_depth: u16,
    pub min_length: usize,
    pub min_overlap: usize,
    pub max_overlaps: usize,
    pub min_identity: f64,
    pub fixed_qual: Option<u8>,
    pub fill_holes: bool,
}

#[derive(Debug, Default, Clone)]
pub struct RoundStats {
    pub reads_in: usize,
    pub bases_in: usize,
    pub reads_with_overlaps: usize,
    pub reads_out: usize,
    pub bases_out: usize,
}

/// Corrects `reads` and returns the corrected reads plus per-round counts. An
/// empty result is not an error: it means no read had enough support.
pub fn correct(
    reads: Vec<Record>,
    params: &Params,
) -> anyhow::Result<(Vec<Record>, Vec<RoundStats>)> {
    let mut current = reads;
    let mut stats = Vec::new();
    for round in 0..params.rounds {
        let last = round + 1 == params.rounds;
        let (next, s) = correct_round(&current, params, last)?;
        stats.push(s);
        current = next;
        if current.len() < 2 {
            break;
        }
    }
    Ok((current, stats))
}

fn correct_round(
    reads: &[Record],
    params: &Params,
    last: bool,
) -> anyhow::Result<(Vec<Record>, RoundStats)> {
    let seqs: Vec<Vec<u8>> = reads.iter().map(|r| r.seq.clone()).collect();
    let mut stats = RoundStats {
        reads_in: reads.len(),
        bases_in: seqs.iter().map(|s| s.len()).sum(),
        ..Default::default()
    };
    let partners = find_overlaps(
        &seqs,
        params.preset,
        params.min_overlap,
        params.max_overlaps,
    )?;
    stats.reads_with_overlaps = partners.iter().filter(|p| !p.is_empty()).count();

    // intermediate rounds keep qualities out of it: they are recomputed at the end
    let fixed_qual = if last { params.fixed_qual } else { Some(b'#') };
    let corrected: Vec<Vec<Record>> = (0..reads.len())
        .into_par_iter()
        .map(|t| -> anyhow::Result<Vec<Record>> {
            // a read needs at least min_depth partners to get any base called
            if partners[t].len() < params.min_depth as usize {
                return Ok(Vec::new());
            }
            let aligner = TargetAligner::new(&seqs[t], params.preset)?;
            let mut pileup = Pileup::new(&seqs[t]);
            for &p in &partners[t] {
                if let Some(a) = aligner.align(&seqs[p], params.min_identity)? {
                    pileup.add(&Aligned {
                        query: &a.query,
                        q_start: a.q_start,
                        t_start: a.t_start,
                        cigar: &a.cigar,
                    });
                }
            }
            let segments = pileup.consensus(
                params.min_depth,
                params.min_length,
                fixed_qual,
                params.fill_holes,
            );
            let many = segments.len() > 1;
            Ok(segments
                .into_iter()
                .enumerate()
                .map(|(k, s)| {
                    let mut name = reads[t].name.clone();
                    if many {
                        name.extend_from_slice(format!("_{}", k + 1).as_bytes());
                    }
                    Record {
                        name,
                        seq: s.seq,
                        qual: Some(s.qual),
                    }
                })
                .collect())
        })
        .collect::<anyhow::Result<_>>()?;
    let out: Vec<Record> = corrected.into_iter().flatten().collect();
    stats.reads_out = out.len();
    stats.bases_out = out.iter().map(|r| r.seq.len()).sum();
    Ok((out, stats))
}
