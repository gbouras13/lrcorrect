//! Overlaps and alignments with minimap2.
//!
//! Correction is two passes. First every read is mapped against an index of all
//! reads without base-level alignment, which is cheap, to find which reads
//! overlap. Then each read in turn becomes a one-sequence index and its chosen
//! partners are aligned to it with a CIGAR. That bounds memory to one read's
//! alignments at a time, and makes each read's correction depend only on its
//! own partners, so the output does not depend on the thread count.

use std::io::Write;

use anyhow::Context;
use minimap2::ffi::{MM_F_NO_DIAG, MM_F_NO_DUAL};
use minimap2::{Aligner, Built, Strand};
use rayon::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Preset {
    Ont,
    Pb,
}

impl Preset {
    fn builder(self) -> Aligner<minimap2::PresetSet> {
        match self {
            Preset::Ont => Aligner::builder().ava_ont(),
            Preset::Pb => Aligner::builder().ava_pb(),
        }
    }
}

fn read_id(i: usize) -> Vec<u8> {
    i.to_string().into_bytes()
}

/// For every read, its overlapping reads, best first. A pair found once
/// (minimap2's all-vs-all mode reports each pair a single time) is recorded for
/// both reads. Overlaps shorter than `min_overlap` bases are ignored, and each
/// read keeps at most `max_overlaps` partners, ranked by minimap2's count of
/// matching bases.
pub fn find_overlaps(
    seqs: &[Vec<u8>],
    preset: Preset,
    min_overlap: usize,
    max_overlaps: usize,
) -> anyhow::Result<Vec<Vec<usize>>> {
    let mut partners: Vec<Vec<(u32, usize)>> = vec![Vec::new(); seqs.len()];
    if seqs.len() < 2 {
        return Ok(vec![Vec::new(); seqs.len()]);
    }
    let ids: Vec<Vec<u8>> = (0..seqs.len()).map(read_id).collect();
    // minimap2-rs's in-memory multi-sequence index (with_seqs_and_ids) hands
    // minimap2 an array of CStrings where it expects char**, which crashes for
    // more than one sequence. Indexing a FASTA file goes through minimap2's own
    // reader instead.
    let dir = tempfile::tempdir().context("cannot create a temporary directory")?;
    let fasta = dir.path().join("reads.fa");
    {
        let mut out = std::io::BufWriter::new(
            std::fs::File::create(&fasta).context("cannot write the temporary index FASTA")?,
        );
        for (id, seq) in ids.iter().zip(seqs) {
            out.write_all(b">")?;
            out.write_all(id)?;
            out.write_all(b"\n")?;
            out.write_all(seq)?;
            out.write_all(b"\n")?;
        }
        out.flush()
            .context("cannot write the temporary index FASTA")?;
    }
    let aligner = preset
        .builder()
        .with_index_threads(rayon::current_num_threads())
        .with_index(&fasta, None)
        .map_err(|e| anyhow::anyhow!("minimap2 could not index the reads: {e}"))?;
    let hits: Vec<Vec<(usize, u32, usize)>> = (0..seqs.len())
        .into_par_iter()
        .map(|q| {
            let mut found = Vec::new();
            if seqs[q].is_empty() {
                return Ok(found);
            }
            let maps = aligner
                .map(&seqs[q], false, false, None, None, Some(&ids[q]))
                .map_err(|e| anyhow::anyhow!("minimap2 failed to map read {q}: {e}"))?;
            for m in maps {
                let t = m.target_id as usize;
                let len = (m.target_end - m.target_start).max(m.query_end - m.query_start);
                if t != q && len as usize >= min_overlap {
                    found.push((t, m.match_len.max(0) as u32, q));
                }
            }
            Ok(found)
        })
        .collect::<anyhow::Result<_>>()?;
    for (t, score, q) in hits.into_iter().flatten() {
        partners[t].push((score, q));
        partners[q].push((score, t));
    }
    Ok(partners
        .into_iter()
        .map(|mut p| {
            // best first; ties by read index, then one entry per partner
            p.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
            let mut seen = std::collections::HashSet::new();
            p.into_iter()
                .filter(|&(_, r)| seen.insert(r))
                .take(max_overlaps)
                .map(|(_, r)| r)
                .collect()
        })
        .collect())
}

/// An alignment of a partner read to a target read, ready for the pileup.
pub struct TargetAlignment {
    /// the partner in the target's orientation
    pub query: Vec<u8>,
    pub q_start: usize,
    pub t_start: usize,
    pub cigar: Vec<(u32, u8)>,
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

pub struct TargetAligner {
    aligner: Aligner<Built>,
}

impl TargetAligner {
    pub fn new(target: &[u8], preset: Preset) -> anyhow::Result<Self> {
        let mut aligner = preset
            .builder()
            .with_cigar()
            .with_seq_and_id(target, b"t")
            .map_err(|e| anyhow::anyhow!("minimap2 could not index a read: {e}"))?;
        // the partner is a different read, so neither all-vs-all filter applies
        aligner.mapopt.flag &= !((MM_F_NO_DIAG | MM_F_NO_DUAL) as i64);
        Ok(TargetAligner { aligner })
    }

    /// The longest alignment of `query` to the target, if it reaches
    /// `min_identity` (matching bases over alignment columns).
    pub fn align(
        &self,
        query: &[u8],
        min_identity: f64,
    ) -> anyhow::Result<Option<TargetAlignment>> {
        if query.is_empty() {
            return Ok(None);
        }
        let maps = self
            .aligner
            .map(query, false, false, None, None, Some(b"q"))
            .map_err(|e| anyhow::anyhow!("minimap2 failed to align a read: {e}"))?;
        let best = maps
            .into_iter()
            .filter(|m| {
                m.alignment
                    .as_ref()
                    .and_then(|a| a.cigar.as_ref())
                    .is_some()
            })
            .filter(|m| m.block_len > 0 && m.match_len as f64 / m.block_len as f64 >= min_identity)
            .max_by_key(|m| (m.block_len, -(m.target_start)));
        Ok(best.map(|m| {
            let cigar = m.alignment.unwrap().cigar.unwrap();
            let (query, q_start) = match m.strand {
                Strand::Forward => (query.to_vec(), m.query_start as usize),
                Strand::Reverse => (revcomp(query), query.len() - m.query_end as usize),
            };
            TargetAlignment {
                query,
                q_start,
                t_start: m.target_start as usize,
                cigar,
            }
        }))
    }
}
