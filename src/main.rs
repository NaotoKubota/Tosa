use std::collections::HashSet;
use log::{info, LevelFilter};

use tosa::cli;
use tosa::data_loader;
use tosa::bam_reader;
use tosa::gtf;
use tosa::output;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse CLI arguments
    let matches = cli::build_cli().get_matches();
    let config = cli::parse_config(&matches);

    // Initialize logger
    if config.verbose {
        env_logger::Builder::from_default_env()
            .filter(None, LevelFilter::Debug)
            .init();
    } else {
        env_logger::Builder::from_default_env()
            .filter(None, LevelFilter::Info)
            .init();
    }

    // Log configuration
    info!("Running tosa v{}", env!("CARGO_PKG_VERSION"));
    info!("Mode: {}", config.mode);
    info!("BAM file: {}", config.bam_file);
    info!("Output prefix: {}", config.output_prefix);
    info!("Minimum anchor length: {}", config.min_anchor_length);
    info!("Minimum intron length: {}", config.min_intron_length);
    info!("Maximum intron length: {}", config.max_intron_length);
    info!("Maximum loci (NH): {}", config.max_loci);
    info!("Strand mode: {}", config.strand_mode);

    // Load cell barcodes of interest (single mode only)
    let cell_barcodes_of_interest = if config.mode == "single" {
        let barcodes = data_loader::load_cell_barcodes(config.cell_barcode_file.as_ref())?;
        info!(
            "Cell barcodes of interest: {}",
            if barcodes.is_empty() {
                "None (processing all reads)".to_string()
            } else {
                format!("{} barcodes", barcodes.len())
            }
        );
        barcodes
    } else {
        HashSet::new()
    };

    // Parse GTF for boundary counting (if provided)
    let boundary_index = if let Some(ref gtf_path) = config.gtf_file {
        info!("GTF file: {}", gtf_path);
        Some(gtf::parse_gtf(gtf_path)?)
    } else {
        None
    };

    // Process BAM file
    let result = bam_reader::process_bam_records(
        &config,
        &cell_barcodes_of_interest,
        boundary_index.as_ref(),
    )?;

    // Write output files
    info!("Writing output files");
    if config.mode == "single" {
        output::write_junction_single(
            &config.output_prefix,
            &result.junction_counts,
            &result.cell_barcodes,
            &result.junction_strands,
        )?;

        if boundary_index.is_some() {
            output::write_boundary_single(
                &config.output_prefix,
                &result.boundary_counts,
                &result.cell_barcodes,
                &result.boundary_types,
                &result.boundary_strands,
            )?;
        }
    } else {
        output::write_junction_bulk(
            &config.output_prefix,
            &result.junction_totals,
            &result.junction_strands,
        )?;

        if boundary_index.is_some() {
            output::write_boundary_bulk(
                &config.output_prefix,
                &result.boundary_totals,
                &result.boundary_types,
                &result.boundary_strands,
            )?;
        }
    }

    info!("Finished processing");
    Ok(())
}
