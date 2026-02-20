//! Common type definitions for Tosa.

use std::fmt;

/// Strand specificity mode for RNA-seq library preparation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StrandMode {
    /// No strand specificity. Uses XS tag if available, otherwise Unknown.
    Unstranded,
    /// Use XS tags provided by aligner.
    XS,
    /// First-strand (RF): read1 reverse = +, read1 forward = -.
    RF,
    /// Second-strand (FR): read1 forward = +, read1 reverse = -.
    FR,
}

impl StrandMode {
    /// Parse strand mode from a string argument.
    pub fn from_str_opt(s: Option<&String>) -> Self {
        match s.map(|v| v.as_str()) {
            Some("RF") => StrandMode::RF,
            Some("FR") => StrandMode::FR,
            Some("XS") => StrandMode::XS,
            _ => StrandMode::Unstranded,
        }
    }
}

impl fmt::Display for StrandMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            StrandMode::Unstranded => write!(f, "unstranded"),
            StrandMode::XS => write!(f, "XS"),
            StrandMode::RF => write!(f, "RF"),
            StrandMode::FR => write!(f, "FR"),
        }
    }
}

/// Strand of a read or feature.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
            Strand::Unknown => write!(f, "."),
        }
    }
}

/// Type of exon-intron boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum BoundaryType {
    /// 5' splice site (exon → intron boundary).
    FivePrime,
    /// 3' splice site (intron → exon boundary).
    ThreePrime,
}

impl fmt::Display for BoundaryType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BoundaryType::FivePrime => write!(f, "5p"),
            BoundaryType::ThreePrime => write!(f, "3p"),
        }
    }
}

/// Configuration for a Tosa run, parsed from CLI arguments.
#[derive(Debug, Clone)]
pub struct RunConfig {
    /// Mode of operation: "bulk" or "single".
    pub mode: String,
    /// Path to the BAM file.
    pub bam_file: String,
    /// Output prefix for output files.
    pub output_prefix: String,
    /// Minimum anchor length for both sides of junctions.
    pub min_anchor_length: i64,
    /// Minimum intron length for junctions.
    pub min_intron_length: i64,
    /// Maximum intron length for junctions.
    pub max_intron_length: i64,
    /// Maximum number of loci the read maps to (NH tag).
    pub max_loci: u32,
    /// Optional path to cell barcode file.
    pub cell_barcode_file: Option<String>,
    /// Strand specificity mode.
    pub strand_mode: StrandMode,
    /// Optional path to GTF annotation file.
    pub gtf_file: Option<String>,
    /// Enable verbose (debug) logging.
    pub verbose: bool,
}
