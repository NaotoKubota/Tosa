# Change log

All notable changes to this Tosa project will be documented in this file.

## [v0.2.1] - 2024-??-??

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
