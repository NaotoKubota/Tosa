# Tosa

Rust implementation of junction read counting from deduplicated scRNA-seq BAM files.

## Usage

```bash
Counts junction reads from scRNA-seq data

Usage: Tosa [OPTIONS] <bam_file> <output_file>

Arguments:
  <bam_file>     Path to the BAM file
  <output_file>  Path to the output file

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
./target/release/Tosa -a 8 -m 70 -M 500000 example.bam example.tsv
```
