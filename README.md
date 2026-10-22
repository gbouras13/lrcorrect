# lrcorrect

Fast, deterministic self-correction of small long-read sets: the reads left over
once a bacterial chromosome is removed, a plasmid, a phage, a small genome.

It was written to replace `canu -correct` in
[plassembler](https://github.com/gbouras13/plassembler), where corrected long
reads are given to Unicycler (`-s`) as if they were short reads. It is one
binary with minimap2 built in: no perl, no java, no shell scripts, no genome
size and no job scheduler.

```
lrcorrect -t 8 --target-bases 7500000 -o corrected.fastq reads.fastq.gz
```

## How it works

1. Optionally, subsample the reads at random to about `--target-bases` bases
   (seeded, so the same reads are picked every time).
2. Find which reads overlap with minimap2's all-vs-all preset (`ava-ont`, or
   `ava-pb` with `--preset pb`).
3. For every read, align its best `--max-overlaps` partners to it with minimap2,
   left-align the indels and call a majority-vote consensus. Bases covered by
   fewer than `--min-depth` other reads are cut out, so a read can come out
   shorter, in pieces, or not at all.
4. Repeat for `--rounds` rounds (default 2), each on the previous round's output.

Each read's correction depends only on the reads aligned to it, so the output is
byte-identical whatever the thread count.

## Output and exit status

Corrected reads are written as FASTQ, with qualities computed from the vote
(`--qual-char` gives every base one character instead), or as FASTA with
`--fasta`. Reads keep their names; a read cut into pieces gets `_1`, `_2`, ...

- Exit 0 with reads: corrected.
- Exit 0 with no reads: there was too little data, for example an empty input,
  a single read, or reads that do not overlap deeply enough. The reason is
  printed to stderr. This is the expected result for a sample with no plasmids.
- Exit 1: something is wrong, for example an unreadable or malformed input file
  or an output that cannot be written. The message on stderr says what.

## Options

```
  -o, --output <OUTPUT>              Output file; "-" for stdout [default: -]
  -t, --threads <THREADS>            Threads [default: 1]
  -p, --preset <PRESET>              Read type: ont or pb [default: ont]
      --target-bases <TARGET_BASES>  Subsample to about this many bases first (0: every read) [default: 0]
      --seed <SEED>                  Seed for --target-bases [default: 11]
  -r, --rounds <ROUNDS>              Rounds of correction [default: 2]
      --min-depth <MIN_DEPTH>        Minimum other reads covering a base to correct it [default: 3]
      --min-length <MIN_LENGTH>      Minimum length of a corrected read [default: 500]
      --min-overlap <MIN_OVERLAP>    Minimum overlap length [default: 500]
      --max-overlaps <MAX_OVERLAPS>  Maximum overlapping reads used per read [default: 60]
      --min-identity <MIN_IDENTITY>  Minimum overlap identity [default: 0.75 ont, 0.65 pb]
      --qual-char <QUAL_CHAR>        Fixed quality character for every base
      --fasta                        Write FASTA instead of FASTQ
```

## Scope

lrcorrect is meant for megabases, not gigabases: every read's overlaps are
found against an in-memory index of all reads. It is not haplotype-aware. PacBio
HiFi reads do not need it.

## Building

```
cargo build --release
cargo test --release
```

minimap2 is compiled from source by the `minimap2-sys` crate, so a C compiler
and zlib are needed.

## Licence

MIT
