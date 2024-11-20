[![GitHub License](https://img.shields.io/github/license/NaotoKubota/Tosa)](https://github.com/NaotoKubota/Tosa/blob/main/LICENSE)
[![Rust](https://github.com/NaotoKubota/Tosa/actions/workflows/rust.yaml/badge.svg)](https://github.com/NaotoKubota/Tosa/actions/workflows/rust.yaml)

# Tosa (v0.1.0)

Rust implementation of junction read counting from RNA-seq BAM files.

## Usage

```bash
Extract junction reads from RNA-seq/scRNA-seq bam files

Usage: tosa [OPTIONS] <mode> <bam_file> <output_prefix>

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
  -l, --max-loci <max_loci>
          Maximum number of loci the read maps to [default: 1]
  -c, --cell-barcodes <cell_barcode_file>
          Optional file specifying cell barcodes of interest
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
./target/release/tosa bulk -a 8 -m 70 -M 500000 example.bam output_example
# Count junction reads from single-cell RNA-seq BAM file
./target/release/tosa single -a 8 -m 70 -M 500000 example.bam output_example
```
