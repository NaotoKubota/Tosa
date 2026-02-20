//! # Tosa
//!
//! Fast junction and exon-intron boundary read counting from RNA-seq/scRNA-seq BAM files.
//!
//! Tosa processes mapped BAM files to extract:
//! - **Junction read counts**: Reads spanning splice junctions (identified by `N` CIGAR operations)
//! - **Boundary read counts**: Reads overlapping exon-intron boundaries (when GTF annotation is provided)
//!
//! Supports both bulk RNA-seq and single-cell RNA-seq (10x Genomics-style barcodes),
//! with configurable strand specificity (unstranded, XS, RF, FR).

pub mod types;
pub mod cli;
pub mod data_loader;
pub mod bam_reader;
pub mod junction;
pub mod boundary;
pub mod gtf;
pub mod output;
