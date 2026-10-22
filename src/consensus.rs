//! Pileup consensus of one read ("the target") from other reads aligned to it.
//!
//! Each aligned read votes for a base, or a deletion, at every target position
//! it covers, and for the bases it inserts between two target positions. Indels
//! are left-aligned first: without that, reads that all carry the same extra A
//! in a homopolymer place it at different offsets, and the vote splits.

pub const GAP: u8 = b'-';
const BASES: [u8; 4] = *b"ACGT";
const DEL: usize = 4;

/// One alignment of a query against the target, in target coordinates. `cigar`
/// holds (length, op) with minimap2's op codes; `query` is already in the
/// target's orientation, and `q_start` is where the alignment begins in it.
pub struct Aligned<'a> {
    pub query: &'a [u8],
    pub q_start: usize,
    pub t_start: usize,
    pub cigar: &'a [(u32, u8)],
}

const CIGAR_M: u8 = 0;
const CIGAR_I: u8 = 1;
const CIGAR_D: u8 = 2;
const CIGAR_EQ: u8 = 7;
const CIGAR_X: u8 = 8;

pub struct Pileup {
    target: Vec<u8>,
    /// votes per position for A, C, G, T and deletion
    counts: Vec<[u16; 5]>,
    /// number of aligned reads covering each position
    depth: Vec<u16>,
    /// (position the bases are inserted before, inserted bases)
    insertions: Vec<(u32, Vec<u8>)>,
}

/// A corrected stretch of the target.
pub struct Segment {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Phred quality, as a FASTQ character, of a call `agree` of `total` votes made.
/// The +1/+2 keeps a unanimous call finite; the cap matches common basecallers.
fn quality(agree: usize, total: usize) -> u8 {
    let err = (total - agree.min(total) + 1) as f64 / (total + 2) as f64;
    let q = (-10.0 * err.log10()).round().clamp(0.0, 40.0) as u8;
    q + 33
}

/// Shifts every gap run in `x` as far left as it can go without changing the
/// alignment score: the column just before the run must be a match whose base
/// also matches the last column of the run.
fn left_align(x: &mut [u8], y: &[u8]) {
    let n = x.len();
    let mut i = 0;
    while i < n {
        if x[i] != GAP {
            i += 1;
            continue;
        }
        let mut end = i;
        while end < n && x[end] == GAP {
            end += 1;
        }
        let (mut a, mut b) = (i, end);
        while a > 0 && x[a - 1] != GAP && x[a - 1] == y[a - 1] && x[a - 1] == y[b - 1] {
            x[b - 1] = x[a - 1];
            x[a - 1] = GAP;
            a -= 1;
            b -= 1;
        }
        i = end;
    }
}

impl Pileup {
    pub fn new(target: &[u8]) -> Self {
        Pileup {
            target: target.to_vec(),
            counts: vec![[0; 5]; target.len()],
            depth: vec![0; target.len()],
            insertions: Vec::new(),
        }
    }

    /// Adds one alignment's votes. Returns false, adding nothing, if the CIGAR
    /// runs off either sequence.
    pub fn add(&mut self, aln: &Aligned) -> bool {
        // expand to two gapped rows over the aligned region
        let mut t_row = Vec::new();
        let mut q_row = Vec::new();
        let (mut t, mut q) = (aln.t_start, aln.q_start);
        for &(len, op) in aln.cigar {
            let len = len as usize;
            match op {
                CIGAR_M | CIGAR_EQ | CIGAR_X => {
                    if t + len > self.target.len() || q + len > aln.query.len() {
                        return false;
                    }
                    t_row.extend_from_slice(&self.target[t..t + len]);
                    q_row.extend_from_slice(&aln.query[q..q + len]);
                    t += len;
                    q += len;
                }
                CIGAR_I => {
                    if q + len > aln.query.len() {
                        return false;
                    }
                    t_row.extend(std::iter::repeat(GAP).take(len));
                    q_row.extend_from_slice(&aln.query[q..q + len]);
                    q += len;
                }
                CIGAR_D => {
                    if t + len > self.target.len() {
                        return false;
                    }
                    t_row.extend_from_slice(&self.target[t..t + len]);
                    q_row.extend(std::iter::repeat(GAP).take(len));
                    t += len;
                }
                _ => {} // clipping and anything else carries no aligned bases
            }
        }
        left_align(&mut q_row, &t_row);
        left_align(&mut t_row, &q_row);

        let mut pos = aln.t_start;
        let mut inserted: Vec<u8> = Vec::new();
        for (&tb, &qb) in t_row.iter().zip(q_row.iter()) {
            if tb == GAP {
                inserted.push(qb);
                continue;
            }
            if !inserted.is_empty() {
                // bases before the first aligned column have no anchor; drop them
                if pos > aln.t_start {
                    self.insertions
                        .push((pos as u32, std::mem::take(&mut inserted)));
                }
                inserted.clear();
            }
            let slot = if qb == GAP { Some(DEL) } else { base_index(qb) };
            if let Some(s) = slot {
                self.counts[pos][s] = self.counts[pos][s].saturating_add(1);
            }
            self.depth[pos] = self.depth[pos].saturating_add(1);
            pos += 1;
        }
        true
    }

    /// Calls the consensus. Positions covered by fewer than `min_depth` reads
    /// are not corrected: the target is cut there, and only stretches of at
    /// least `min_len` corrected bases are returned. The target's own base
    /// breaks ties. With `fixed_qual`, every base gets that quality character.
    pub fn consensus(
        mut self,
        min_depth: u16,
        min_len: usize,
        fixed_qual: Option<u8>,
        fill_holes: bool,
    ) -> Vec<Segment> {
        self.insertions.sort_unstable();
        // with fill_holes, unsupported stretches between the first and last
        // supported positions keep the target's own bases instead of cutting it
        let first = self.depth.iter().position(|&d| d >= min_depth);
        let last = self.depth.iter().rposition(|&d| d >= min_depth);
        let mut segments = Vec::new();
        let mut seq = Vec::new();
        let mut qual = Vec::new();
        let mut ins_i = 0;
        let mut prev_in_segment = false;
        let push = |segments: &mut Vec<Segment>, seq: &mut Vec<u8>, qual: &mut Vec<u8>| {
            if seq.len() >= min_len {
                segments.push(Segment {
                    seq: std::mem::take(seq),
                    qual: std::mem::take(qual),
                });
            }
            seq.clear();
            qual.clear();
        };
        for pos in 0..self.target.len() {
            // the insertions anchored before this position
            let ins_start = ins_i;
            while ins_i < self.insertions.len() && self.insertions[ins_i].0 as usize == pos {
                ins_i += 1;
            }
            let depth = self.depth[pos];
            let inside = matches!((first, last), (Some(f), Some(l)) if f <= pos && pos <= l);
            if depth < min_depth && fill_holes && inside {
                seq.push(self.target[pos]);
                qual.push(fixed_qual.unwrap_or(b'!'));
                // an uncorrected base calls no insertion before the next one
                prev_in_segment = false;
                continue;
            }
            if depth < min_depth {
                if !seq.is_empty() {
                    push(&mut segments, &mut seq, &mut qual);
                }
                prev_in_segment = false;
                continue;
            }
            if prev_in_segment && ins_i > ins_start {
                let span = depth.min(self.depth[pos - 1]) as usize;
                let n_ins = ins_i - ins_start;
                if 2 * n_ins > span {
                    let (best, support) = most_common(&self.insertions[ins_start..ins_i]);
                    for &b in best {
                        seq.push(b);
                        qual.push(fixed_qual.unwrap_or_else(|| quality(support, span)));
                    }
                }
            }
            let c = &self.counts[pos];
            let own = base_index(self.target[pos]);
            let mut best = own.unwrap_or(DEL);
            for k in 0..5 {
                if c[k] > c[best] {
                    best = k;
                }
            }
            if best != DEL {
                let total: usize = c.iter().map(|&v| v as usize).sum();
                seq.push(BASES[best]);
                qual.push(fixed_qual.unwrap_or_else(|| quality(c[best] as usize, total)));
            }
            prev_in_segment = true;
        }
        push(&mut segments, &mut seq, &mut qual);
        segments
    }
}

/// The most frequent inserted string of a sorted run, with its count. Ties go
/// to the first in sort order, which keeps the call deterministic.
fn most_common(run: &[(u32, Vec<u8>)]) -> (&[u8], usize) {
    let mut best: (&[u8], usize) = (&run[0].1, 0);
    let mut i = 0;
    while i < run.len() {
        let mut j = i;
        while j < run.len() && run[j].1 == run[i].1 {
            j += 1;
        }
        if j - i > best.1 {
            best = (&run[i].1, j - i);
        }
        i = j;
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aligned<'a>(query: &'a [u8], cigar: &'a [(u32, u8)]) -> Aligned<'a> {
        Aligned {
            query,
            q_start: 0,
            t_start: 0,
            cigar,
        }
    }

    #[test]
    fn left_align_moves_homopolymer_gap() {
        let mut q = b"AAA-C".to_vec();
        let t = b"AAAAC".to_vec();
        left_align(&mut q, &t);
        assert_eq!(q, b"-AAAC");
    }

    #[test]
    fn majority_fixes_substitution_deletion_and_insertion() {
        //            target has an error at 3 (T for G), an extra C at 6, and is missing an A after 9
        let target = b"ACGTTACCGTTA";
        let truth = b"ACGGTACGTTAA";
        let mut p = Pileup::new(target);
        // truth aligned to target: ACG [G/T] TAC [C deleted] GTTA [+A]
        let cig = [(7u32, 0u8), (1, 2), (4, 0), (1, 1)];
        for _ in 0..5 {
            assert!(p.add(&aligned(truth, &cig)));
        }
        let seg = p.consensus(3, 1, None, false);
        assert_eq!(seg.len(), 1);
        // left-alignment anchors the trailing insertion before the last A
        assert_eq!(seg[0].seq, truth);
        assert_eq!(seg[0].qual.len(), seg[0].seq.len());
    }

    #[test]
    fn insertion_votes_agree_after_left_alignment() {
        let target = b"CCAAAGG";
        // two reads put the extra A at different offsets in the run
        let q = b"CCAAAAGG";
        let early = [(2u32, 0u8), (1, 1), (5, 0)];
        let late = [(5u32, 0u8), (1, 1), (2, 0)];
        let mut p = Pileup::new(target);
        p.add(&aligned(q, &early));
        p.add(&aligned(q, &late));
        p.add(&aligned(q, &late));
        let seg = p.consensus(3, 1, Some(b'#'), false);
        assert_eq!(seg[0].seq, b"CCAAAAGG");
        assert!(seg[0].qual.iter().all(|&c| c == b'#'));
    }

    #[test]
    fn low_depth_splits_and_trims() {
        let target = b"AAAACCCCGGGGTTTT";
        let mut p = Pileup::new(target);
        let left = [(6u32, 0u8)];
        let right = Aligned {
            query: &target[10..],
            q_start: 0,
            t_start: 10,
            cigar: &[(6u32, 0u8)],
        };
        for _ in 0..3 {
            p.add(&aligned(target, &left));
            p.add(&right);
        }
        let seg = p.consensus(3, 2, None, false);
        let seqs: Vec<&[u8]> = seg.iter().map(|s| &s.seq[..]).collect();
        assert_eq!(seqs, vec![&b"AAAACC"[..], &b"GGTTTT"[..]]);
        // and a minimum length drops both
        let mut p = Pileup::new(target);
        p.add(&aligned(target, &left));
        assert!(p.consensus(3, 1, None, false).is_empty());
    }

    #[test]
    fn fill_holes_keeps_unsupported_middle_but_trims_ends() {
        let target = b"AAAACCCCGGGGTTTT";
        let left = [(6u32, 0u8)];
        let right = Aligned {
            query: &target[10..],
            q_start: 0,
            t_start: 10,
            cigar: &[(6u32, 0u8)],
        };
        let mut p = Pileup::new(target);
        for _ in 0..3 {
            p.add(&aligned(target, &left));
            p.add(&right);
        }
        let seg = p.consensus(3, 2, None, true);
        assert_eq!(seg.len(), 1);
        assert_eq!(seg[0].seq, target);
        // the uncorrected middle is marked with the lowest quality
        assert_eq!(&seg[0].qual[6..10], b"!!!!");
    }

    #[test]
    fn cigar_past_the_end_is_rejected() {
        let mut p = Pileup::new(b"ACGT");
        assert!(!p.add(&aligned(b"ACGTACGT", &[(8u32, 0u8)])));
        assert!(p.consensus(1, 1, None, false).is_empty());
    }
}
