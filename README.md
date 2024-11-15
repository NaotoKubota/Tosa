# Tosa (v0.1.0)

Rust implementation of junction read counting from RNA-seq BAM files.

## Usage

```bash
Counts junction reads from RNA-seq data

Usage: Tosa [OPTIONS] <mode> <bam_file> <output_prefix>

Arguments:
  <mode>           Mode of operation: 'bulk' or 'single' [possible values: bulk, single]
  <bam_file>       Path to the BAM file
  <output_prefix>  Output prefix for the output files

Options:
  -a, --anchor-length <anchor_length>
          Minimum anchor length for both sides of junctions [default: 8]
  -m, --min-intron-length <min_intron_length>
          Minimum intron length for junctions [default: 70]
  -M, --max-intron-length <max_intron_length>
          Maximum intron length for junctions [default: 500000]
  -v, --verbose
          Enable verbose output to print all arguments
  -h, --help
          Print help
  -V, --version
          Print version
```

## Build

```bash
cargo build --release
```

## Example

```bash
# Count junction reads from bulk RNA-seq BAM file
./target/release/Tosa bulk -a 8 -m 70 -M 500000 example.bam example.tsv
# Count junction reads from single-cell RNA-seq BAM file
./target/release/Tosa single -a 8 -m 70 -M 500000 example.bam example.tsv
```
