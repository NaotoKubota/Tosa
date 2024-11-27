[![GitHub License](https://img.shields.io/github/license/NaotoKubota/Tosa)](https://github.com/NaotoKubota/Tosa/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/888699325.svg)](https://doi.org/10.5281/zenodo.14202094)
[![GitHub Release](https://img.shields.io/github/v/release/NaotoKubota/Tosa?style=flat)](https://github.com/NaotoKubota/Tosa/releases)
[![GitHub Release Date](https://img.shields.io/github/release-date/NaotoKubota/Tosa)](https://github.com/NaotoKubota/Tosa/releases)
[![Rust](https://github.com/NaotoKubota/Tosa/actions/workflows/rust.yaml/badge.svg)](https://github.com/NaotoKubota/Tosa/actions/workflows/rust.yaml)
[![Create Release and Build Docker Image](https://github.com/NaotoKubota/Tosa/actions/workflows/release-docker-build-push.yaml/badge.svg)](https://github.com/NaotoKubota/Tosa/actions/workflows/release-docker-build-push.yaml)
[![Docker Pulls](https://img.shields.io/docker/pulls/naotokubota/tosa)](https://hub.docker.com/r/naotokubota/tosa)
[![Docker Image Size](https://img.shields.io/docker/image-size/naotokubota/tosa)](https://hub.docker.com/r/naotokubota/tosa)

# Tosa (v0.3.0)

Rust implementation of junction read counting from RNA-seq BAM files.

## Usage

```bash
Extract junction reads from RNA-seq/scRNA-seq bam files

Usage: tosa [OPTIONS] <mode> <bam_file> <output_dir>

Arguments:
  <mode>        Mode of operation: 'bulk' or 'single' [possible values: bulk, single]
  <bam_file>    Path to the BAM file
  <output_dir>  Output directory for the output files

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
./target/release/tosa bulk example.bam output_example
# Count junction reads from single-cell RNA-seq BAM file
./target/release/tosa single example.bam output_example
```
