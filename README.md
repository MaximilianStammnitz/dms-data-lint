# dms-data-lint

`dms_data_lint` checks that the two FASTQ files of a paired-end library are intact and consistent before anything else reads them. It decompresses each file completely, checking every gzip member's checksum and length, and checks every record of both files: the format of each record, equal numbers of records, mate names that match in the same order, and no read name twice. It is the raw-data gate of the nf-core pipeline [nf-core/deepmutscan](https://nf-co.re/deepmutscan) from v1.1.0 on.

The checks, their order and their messages are those of [`fq lint`](https://github.com/stjude-rust-labs/fq) (St. Jude Rust Labs) in its default settings, version 0.12.0, so that the two give the same verdict and can be compared case by case. `dms_data_lint` is a separate implementation in C and uses no code of fq. It decompresses with [ISA-L](https://github.com/intel/isa-l), reads the two files in two threads, and checks characters with vector instructions:

| 6 million read pairs, gzip, 2 CPUs | wall time | CPU time | peak memory |
| --- | --- | --- | --- |
| `fq lint` 0.12.0 | 15.2 s | 14.9 s | 27 MB |
| `dms_data_lint` 0.1.0 | 2.0 s | 3.5 s | 103 MB |

Memory grows by 8 bytes per read pair (about 0.8 GB for 94 million pairs): the duplicate-name check keeps a 64-bit hash of every name, so it needs no second pass over the first file.

## Installation

From Bioconda:

```bash
conda install -c bioconda -c conda-forge dms-data-lint
```

From source, with [ISA-L](https://github.com/intel/isa-l) installed (Debian/Ubuntu: `libisal-dev`; Homebrew: `isa-l`; conda-forge: `isa-l`):

```bash
make -f build.mk ISAL_PREFIX=/usr      # prefix of the ISA-L installation (include/, lib/)
bash test/run_test.sh                  # checks the build against the bundled test data
```

## Usage

```bash
dms_data_lint sample_R1.fastq.gz sample_R2.fastq.gz
dms_data_lint --max-base-qual sample.max_base_qual.tsv sample_R1.fastq.gz sample_R2.fastq.gz
```

| option | meaning |
| --- | --- |
| `R1 R2` | the two files of a read pair; gzip-compressed when the name ends in `.gz`, else plain text |
| `--max-base-qual FILE` | for a pair that passes: the highest base quality (Phred, quality characters offset by 33) among the first 10,000 records of R1 and of R2, tab-separated. nf-core/deepmutscan uses it to warn when an instrument's qualities stay below the counting threshold. |
| `--version` | the version of `dms_data_lint` and of ISA-L |
| `--help` | usage |

Exit status 0: the pair passes (one summary line on standard error). 1: it does not, with one message on standard error. 2: wrong usage.

## Checks

For each pair of records, in this order, R1 before R2 for each check (the codes are fq's):

| code | check |
| --- | --- |
| S003 | the name line starts with `@` |
| S004 | name, sequence, plus line and quality are not empty |
| S002 | the sequence consists of `ACGTNacgtn` |
| S001 | the plus line starts with `+` |
| S005 | sequence and quality have the same length |
| S006 | quality characters are printable ASCII (`!` to `~`) |
| P001 | R1 and R2 have the same name |

Then, once all records passed: S007, no name occurs twice in R1. Names are compared after truncation at their last `/` or space, so `@read/1` and `@read/2`, or `@read 1:N:0` and `@read 2:N:0`, are mates. Lines end at a line feed; one carriage return before it is dropped. Records whose last lines are missing count as records with empty lines. The first problem stops the check.

Messages:

```
sample_R1.fastq.gz:30:11: [S002] AlphabetValidator: invalid sequence character: 'R'
Error: r2-src unexpectedly ended
Error: I/O error in sample_R1.fastq.gz: corrupt gzip stream does not have a matching checksum
```

The first gives file, line and column. The second says that R2 has fewer records than R1 (`r1-src` the reverse). The third is a file that cannot be read or decompressed, naming the file (fq's message does not).

## Tests

- `test/run_test.sh`: 5,000 read pairs of the nf-core/deepmutscan test data (GID1A, public), and 20 small cases in `test/cases/` (made by `test/make_cases.py`) with the expected exit status and message of fq lint 0.12.0 on the same files.
- `test/compare_with_fq/` (needs fq): `equivalence.py` runs both tools on 122 constructed pairs (every check in R1 and in R2, the order of checks within and across records, name truncation, CRLF, missing final newline, blank lines, truncated records, duplicates, multi-member and damaged gzip, a record longer than the read buffer); `fuzz.py` on randomly damaged pairs. With fq 0.12.0: all 122 identical; 1,000 of 1,000 randomly edited plain files identical; 1,000 of 1,000 damaged gzip files rejected by both. For about four in ten of the latter the first message differs: fq's gzip decoder carries on through some kinds of damage (for example, a back-reference beyond the start of the data becomes zeros) and reports a record check on the garbage, while ISA-L stops at the damage and `dms_data_lint` reports a read error.

## Citation

If you use `dms_data_lint`, please cite nf-core/deepmutscan:

> Wehnert B, et al. nf-core/deepmutscan (preprint in preparation). https://nf-co.re/deepmutscan

and fq, whose checks it implements:

> St. Jude Rust Labs. fq: FASTQ file manipulation and validation. https://github.com/stjude-rust-labs/fq

## License

MIT, see [LICENSE](LICENSE).
