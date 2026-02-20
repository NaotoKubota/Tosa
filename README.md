# Tosa

[![GitHub Release](https://img.shields.io/github/v/release/NaotoKubota/Tosa?style=flat)](https://github.com/NaotoKubota/Tosa/releases)
[![GitHub Release Date](https://img.shields.io/github/release-date/NaotoKubota/Tosa)](https://github.com/NaotoKubota/Tosa/releases)
[![Create Release](https://github.com/NaotoKubota/Tosa/actions/workflows/release.yml/badge.svg)](https://github.com/NaotoKubota/Tosa/actions/workflows/release.yml)
[![Rust](https://github.com/NaotoKubota/Tosa/actions/workflows/rust.yaml/badge.svg)](https://github.com/NaotoKubota/Tosa/actions/workflows/rust.yaml)
[![Publish to crates.io](https://github.com/NaotoKubota/Tosa/actions/workflows/publish.yml/badge.svg)](https://github.com/NaotoKubota/Tosa/actions/workflows/publish.yml)
[![crates.io](https://img.shields.io/crates/v/tosa)](https://crates.io/crates/tosa)
[![GitHub License](https://img.shields.io/github/license/NaotoKubota/Tosa)](https://github.com/NaotoKubota/Tosa/blob/main/LICENSE)
[![Docker](https://img.shields.io/docker/v/naotokubota/tosa?color=blue&label=Docker)](https://hub.docker.com/r/naotokubota/tosa)
[![Docker Pulls](https://img.shields.io/docker/pulls/naotokubota/tosa)](https://hub.docker.com/r/naotokubota/tosa)
[![Docker Image Size](https://img.shields.io/docker/image-size/naotokubota/tosa)](https://hub.docker.com/r/naotokubota/tosa)

Fast junction and exon-intron boundary read counting from RNA-seq/scRNA-seq BAM files.

## Features

- **Junction read counting** from spliced alignments (CIGAR `N` operations)
- **Exon-intron boundary read counting** using GTF annotation (1-base boundary coordinates)
- **Strand specificity** support: unstranded, XS tag, RF (first-strand), FR (second-strand)
- **Bulk and single-cell** modes (10x Genomics-style cell barcodes)
- Paired-end read deduplication (same junction/boundary counted once per read pair)

## Usage

```
Extract junction and boundary reads from RNA-seq/scRNA-seq BAM files

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
  -s, --strand <strand>
          Strand specificity of RNA library: RF (first-strand), FR (second-strand),
          XS (use XS tags). Omit for unstranded [possible values: RF, FR, XS]
  -g, --gtf <gtf_file>
          GTF annotation file for exon-intron boundary read counting
  -v, --verbose
          Enable verbose output to print all arguments
  -h, --help
          Print help
  -V, --version
          Print version
```

## Installation

```bash
# Build and install to ~/.cargo/bin (make sure ~/.cargo/bin is in your PATH)
cargo install --path .
```

After installation you can run `tosa` directly from anywhere:

```bash
tosa --version
```

> [!TIP]
> If you only want to build without installing, use `cargo build --release`.
> The binary will be at `./target/release/tosa`.

## Quick start

Example data is included in the `examples/` directory so you can try Tosa right away:

```bash
# Count junction reads (bulk mode)
tosa bulk examples/example.bam output_example

# Count junction + boundary reads with GTF annotation
tosa bulk -g examples/annotation.gtf examples/example.bam output_example

# View results
gunzip -c output_example_junction.tsv.gz | head
gunzip -c output_example_boundary.tsv.gz | head
```

## More examples

```bash
# Count junction reads with strand specificity (first-strand library)
tosa bulk -s RF input.bam output_prefix

# Count junction + boundary reads with GTF annotation and strand
tosa bulk -s RF -g annotation.gtf input.bam output_prefix

# Count junction reads from single-cell RNA-seq BAM file
tosa single -c barcodes.tsv input.bam output_prefix
```

## Output

### Junction output (bulk mode)

- `{prefix}_junction.tsv.gz`: Tab-separated file with columns `Junction`, `Strand`, `Count`

### Junction output (single-cell mode)

- `{prefix}_matrix.mtx.gz`: MatrixMarket sparse matrix
- `{prefix}_barcodes.tsv.gz`: Cell barcodes
- `{prefix}_features.tsv.gz`: Junction coordinates with strand
- `{prefix}_junction_barcodes.tsv.gz`: Detailed TSV with `Feature`, `Strand`, `Barcode`, `Count`

### Boundary output (when `--gtf` is provided)

- Bulk: `{prefix}_boundary.tsv.gz` with columns `Boundary`, `Type`, `Strand`, `Count`
- Single: `{prefix}_boundary_matrix.mtx.gz`, `{prefix}_boundary_barcodes.tsv.gz`, `{prefix}_boundary_features.tsv.gz`

Boundary coordinates use 1-base intervals at intron ends. For example, intron `chr2:6545675-6547042` produces:

- 5' boundary: `chr2:6545675-6545676`
- 3' boundary: `chr2:6547041-6547042`

## License

MIT License

## Contributing

Thank you for wanting to improve Tosa! If you have any bugs or questions, feel free to [open an issue](https://github.com/NaotoKubota/Tosa/issues) or pull request.

## Authors

- Naoto Kubota ([0000-0003-0612-2300](https://orcid.org/0000-0003-0612-2300))
