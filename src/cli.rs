//! CLI argument definition for Tosa.

use clap::{Arg, Command};
use crate::types::{RunConfig, StrandMode};

/// Build the CLI command definition.
pub fn build_cli() -> Command {
    Command::new("tosa")
        .version(env!("CARGO_PKG_VERSION"))
        .author("NaotoKubota")
        .about("Extract junction and boundary reads from RNA-seq/scRNA-seq BAM files")
        .arg(Arg::new("mode")
            .required(true)
            .value_parser(["bulk", "single"])
            .help("Mode of operation: 'bulk' or 'single'"))
        .arg(Arg::new("bam_file")
            .required(true)
            .help("Path to the BAM file"))
        .arg(Arg::new("output_prefix")
            .required(true)
            .help("Output prefix for the output files"))
        .arg(Arg::new("anchor_length")
            .short('a')
            .long("anchor-length")
            .default_value("8")
            .value_parser(clap::value_parser!(i64))
            .help("Minimum anchor length for both sides of junctions"))
        .arg(Arg::new("min_intron_length")
            .short('m')
            .long("min-intron-length")
            .default_value("70")
            .value_parser(clap::value_parser!(i64))
            .help("Minimum intron length for junctions"))
        .arg(Arg::new("max_intron_length")
            .short('M')
            .long("max-intron-length")
            .default_value("500000")
            .value_parser(clap::value_parser!(i64))
            .help("Maximum intron length for junctions"))
        .arg(Arg::new("max_loci")
            .short('l')
            .long("max-loci")
            .default_value("1")
            .value_parser(clap::value_parser!(u32))
            .help("Maximum number of loci the read maps to"))
        .arg(Arg::new("cell_barcode_file")
            .short('c')
            .long("cell-barcodes")
            .value_parser(clap::value_parser!(String))
            .help("Optional file specifying cell barcodes of interest"))
        .arg(Arg::new("strand")
            .short('s')
            .long("strand")
            .value_parser(["RF", "FR", "XS"])
            .help("Strand specificity of RNA library: RF (first-strand), FR (second-strand), XS (use XS tags). Omit for unstranded"))
        .arg(Arg::new("gtf_file")
            .short('g')
            .long("gtf")
            .value_parser(clap::value_parser!(String))
            .help("GTF annotation file for exon-intron boundary read counting"))
        .arg(Arg::new("verbose")
            .short('v')
            .long("verbose")
            .action(clap::ArgAction::SetTrue)
            .help("Enable verbose output to print all arguments"))
}

/// Parse CLI matches into a RunConfig.
pub fn parse_config(matches: &clap::ArgMatches) -> RunConfig {
    RunConfig {
        mode: matches.get_one::<String>("mode").unwrap().clone(),
        bam_file: matches.get_one::<String>("bam_file").unwrap().clone(),
        output_prefix: matches.get_one::<String>("output_prefix").unwrap().clone(),
        min_anchor_length: *matches.get_one::<i64>("anchor_length").unwrap(),
        min_intron_length: *matches.get_one::<i64>("min_intron_length").unwrap(),
        max_intron_length: *matches.get_one::<i64>("max_intron_length").unwrap(),
        max_loci: *matches.get_one::<u32>("max_loci").unwrap(),
        cell_barcode_file: matches.get_one::<String>("cell_barcode_file").cloned(),
        strand_mode: StrandMode::from_str_opt(matches.get_one::<String>("strand")),
        gtf_file: matches.get_one::<String>("gtf_file").cloned(),
        verbose: matches.get_flag("verbose"),
    }
}
