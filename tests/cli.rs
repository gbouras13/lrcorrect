//! End-to-end tests of the lrcorrect binary on simulated reads.

use std::fs;
use std::path::Path;
use std::process::{Command, Output};

use lrcorrect::io::read_all;
use lrcorrect::subsample::SplitMix64;

fn lrcorrect(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_lrcorrect"))
        .args(args)
        .output()
        .expect("cannot run lrcorrect")
}

fn random_genome(len: usize, rng: &mut SplitMix64) -> Vec<u8> {
    (0..len).map(|_| b"ACGT"[rng.below(4) as usize]).collect()
}

/// Reads from a circular genome with about `error_rate` errors per base,
/// split evenly between substitutions, insertions and deletions, half of them
/// reverse complemented.
fn simulate(
    genome: &[u8],
    n: usize,
    len: usize,
    error_rate: f64,
    rng: &mut SplitMix64,
) -> Vec<Vec<u8>> {
    let mut reads = Vec::new();
    for _ in 0..n {
        let start = rng.below(genome.len() as u64) as usize;
        let mut read = Vec::with_capacity(len + len / 10);
        for k in 0..len {
            let b = genome[(start + k) % genome.len()];
            let r = rng.below(1_000_000) as f64 / 1e6;
            if r < error_rate / 3.0 {
                read.push(b"ACGT"[rng.below(4) as usize]); // substitution (sometimes silent)
            } else if r < 2.0 * error_rate / 3.0 {
                read.push(b);
                read.push(b"ACGT"[rng.below(4) as usize]); // insertion
            } else if r < error_rate {
                // deletion
            } else {
                read.push(b);
            }
        }
        if rng.below(2) == 1 {
            read = lrcorrect::overlap::revcomp(&read);
        }
        reads.push(read);
    }
    reads
}

fn write_fastq(path: &Path, reads: &[Vec<u8>]) {
    let mut s = String::new();
    for (i, r) in reads.iter().enumerate() {
        s.push_str(&format!(
            "@read{i}\n{}\n+\n{}\n",
            String::from_utf8_lossy(r),
            "5".repeat(r.len())
        ));
    }
    fs::write(path, s).unwrap();
}

/// Fraction of a read's 15-mers found in the genome (either strand): a cheap,
/// alignment-free accuracy proxy.
fn kmer_accuracy(reads: &[Vec<u8>], genome: &[u8]) -> f64 {
    let k = 15;
    let mut circ = genome.to_vec();
    circ.extend_from_slice(&genome[..k]);
    let mut set = std::collections::HashSet::new();
    for s in [circ.clone(), lrcorrect::overlap::revcomp(&circ)] {
        for w in s.windows(k) {
            set.insert(w.to_vec());
        }
    }
    let (mut hit, mut total) = (0usize, 0usize);
    for r in reads {
        for w in r.windows(k) {
            total += 1;
            hit += set.contains(w) as usize;
        }
    }
    hit as f64 / total as f64
}

#[test]
fn corrects_and_is_thread_count_independent() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = SplitMix64::new(7);
    let genome = random_genome(20_000, &mut rng);
    let reads = simulate(&genome, 150, 4_000, 0.08, &mut rng);
    let input = dir.path().join("reads.fq");
    write_fastq(&input, &reads);

    let mut outputs = Vec::new();
    for threads in ["1", "4"] {
        let out = dir.path().join(format!("out{threads}.fq"));
        let res = lrcorrect(&[
            "-t",
            threads,
            "-o",
            out.to_str().unwrap(),
            input.to_str().unwrap(),
        ]);
        assert!(
            res.status.success(),
            "{}",
            String::from_utf8_lossy(&res.stderr)
        );
        outputs.push(fs::read(&out).unwrap());
    }
    assert_eq!(outputs[0], outputs[1], "output depends on the thread count");

    let corrected: Vec<Vec<u8>> = read_all(&dir.path().join("out1.fq"))
        .unwrap()
        .into_iter()
        .map(|r| r.seq)
        .collect();
    let before = kmer_accuracy(&reads, &genome);
    let after = kmer_accuracy(&corrected, &genome);
    let bases: usize = corrected.iter().map(|r| r.len()).sum();
    assert!(before < 0.5, "simulation too easy: {before}");
    assert!(
        after > 0.97,
        "15-mer accuracy only {after} after correction ({before} before)"
    );
    assert!(
        bases > 400_000,
        "only {bases} bp of the 600 kbp survived correction"
    );
}

#[test]
fn subsamples_to_target() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = SplitMix64::new(3);
    let genome = random_genome(10_000, &mut rng);
    let reads = simulate(&genome, 100, 3_000, 0.05, &mut rng);
    let input = dir.path().join("reads.fq");
    write_fastq(&input, &reads);
    let out = dir.path().join("out.fq");
    let res = lrcorrect(&[
        "--target-bases",
        "100000",
        "-r",
        "1",
        "-o",
        out.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert!(res.status.success());
    let stderr = String::from_utf8_lossy(&res.stderr);
    assert!(stderr.contains("subsampled to 34 reads"), "{stderr}");
}

#[test]
fn too_little_data_is_success_with_no_reads() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.fq");
    for content in ["", "@only\nACGTACGTAC\n+\nIIIIIIIIII\n"] {
        let input = dir.path().join("in.fq");
        fs::write(&input, content).unwrap();
        let res = lrcorrect(&["-o", out.to_str().unwrap(), input.to_str().unwrap()]);
        assert!(
            res.status.success(),
            "{}",
            String::from_utf8_lossy(&res.stderr)
        );
        assert_eq!(fs::read(&out).unwrap().len(), 0);
        assert!(String::from_utf8_lossy(&res.stderr).contains("no reads could be corrected"));
    }
}

#[test]
fn broken_input_is_an_error() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.fq");
    fs::write(&input, "not a fastq file\n").unwrap();
    let res = lrcorrect(&[input.to_str().unwrap()]);
    assert!(!res.status.success());
    assert!(String::from_utf8_lossy(&res.stderr).contains("lrcorrect: error:"));

    let res = lrcorrect(&["/no/such/file.fq"]);
    assert_eq!(res.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&res.stderr).contains("cannot open /no/such/file.fq"));
}
