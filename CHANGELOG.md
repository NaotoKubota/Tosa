# Change log

All notable changes to this Tosa project will be documented in this file.

<<<<<<< HEAD
## [v0.4.0] - 2026-XX-XX

### Added

- Add `--strand` option for strand specificity: RF (first-strand), FR (second-strand), XS (use XS tags), or omit for unstranded.
- Add `--gtf` option for exon-intron boundary read counting using GTF annotation.
- Add boundary read counting with 1-base boundary coordinates at intron ends (5' and 3' splice sites).
- Add Strand column to junction and boundary output files.
- Add unit tests for junction counting, boundary counting, GTF parsing, and strand determination.
- Add integration tests using test BAM file.
- Add version synchronization tests (version-sync).
- Add `build.rs` for automatic VERSION file generation from Cargo.toml.

### Changed

- Major refactoring: split monolithic main.rs into modular files (cli.rs, types.rs, bam_reader.rs, junction.rs, boundary.rs, gtf.rs, output.rs).
- Add library crate (lib.rs) for better testability and potential crate reuse.
- Junction output now includes Strand column: `Junction\tStrand\tCount`.
- Single-cell features.tsv.gz now includes strand information.

### Fixed

- Fixed `Cigar::Ins` incorrectly advancing reference position (insertions don't consume reference bases).
- Fixed off-by-one error in single-cell mode barcode indexing for MatrixMarket output.

=======
>>>>>>> 33a432e75c1c3be5acd1d209fd1f1af23b49563d
## [v0.3.0] - 2024-11-27

### Added

### Changed

- Output directory is now specified by `--output-dir` instead of `--output-prefix`

### Fixed

- Fix barcode index mapping in sparse matrix data processing

## [v0.2.0] - 2024-11-21

### Added

- Add `--cell-barcodes` option to count reads by cell barcodes of interest.
- Add `--max-loci` option to count reads that map to a maximum number of loci.

### Changed

- Count paired reads that span the same junction only once, not twice.
- Refactor anchor length calculation in junction processing

### Fixed

- Fixed a bug that reads having softclip are not counted

## [v0.1.0] - 2024-11-16

### Added

- Initial release

### Changed

### Fixed
