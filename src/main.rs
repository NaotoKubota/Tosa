use std::collections::HashSet;
use log::{info, LevelFilter};

use tosa::cli;
use tosa::data_loader;
use tosa::bam_reader;
use tosa::gtf;
use tosa::output;

fn main() -> Result<(), Box<dyn std::error::Error>> {
<<<<<<< HEAD
    // Parse CLI arguments
    let matches = cli::build_cli().get_matches();
    let config = cli::parse_config(&matches);

    // Initialize logger
    if config.verbose {
=======
    // Set up command-line arguments using clap
    let matches = Command::new("tosa")
        .version("0.3.0")
        .author("NaotoKubota")
        .about("Extract junction reads from RNA-seq/scRNA-seq bam files")
        .arg(Arg::new("mode")
            .required(true)
            .value_parser(["bulk", "single"])
            .help("Mode of operation: 'bulk' or 'single'"))
        .arg(Arg::new("bam_file")
            .required(true)
            .help("Path to the BAM file"))
        .arg(Arg::new("output_dir")
            .required(true)
            .help("Output directory for the output files"))
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
        .arg(Arg::new("verbose")
            .short('v')
            .long("verbose")
            .action(clap::ArgAction::SetTrue)
            .help("Enable verbose output to print all arguments"))
        .get_matches();

    // Parse arguments
    let mode = matches.get_one::<String>("mode").unwrap();
    let bam_file = matches.get_one::<String>("bam_file").unwrap();
    let output_dir = matches.get_one::<String>("output_dir").unwrap();
    let cell_barcode_file = matches.get_one::<String>("cell_barcode_file");
    let min_anchor_length = *matches.get_one::<i64>("anchor_length").unwrap();
    let min_intron_length = *matches.get_one::<i64>("min_intron_length").unwrap();
    let max_intron_length = *matches.get_one::<i64>("max_intron_length").unwrap();
    let max_loci = *matches.get_one::<u32>("max_loci").unwrap();
    let verbose = matches.get_flag("verbose");

    // Initialize the logger with the appropriate level
    if verbose {
>>>>>>> 33a432e75c1c3be5acd1d209fd1f1af23b49563d
        env_logger::Builder::from_default_env()
            .filter(None, LevelFilter::Debug)
            .init();
    } else {
        env_logger::Builder::from_default_env()
            .filter(None, LevelFilter::Info)
            .init();
    }

<<<<<<< HEAD
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
=======
    // Log all arguments if verbose is enabled
    info!("Running tosa");
    info!("Mode: {}", mode);
    info!("BAM file: {}", bam_file);
    info!("Output prefix: {}", output_dir);
    info!("Minimum anchor length: {}", min_anchor_length);
    info!("Minimum intron length: {}",min_intron_length);
    info!("Maximum intron length: {}", max_intron_length);
    info!("Maximum loci (NH): {}", max_loci);
    // Load cell barcodes of interest
    let cell_barcodes_of_interest = if mode == "single" {
        let barcodes = data_loader::load_cell_barcodes(cell_barcode_file)?;
>>>>>>> 33a432e75c1c3be5acd1d209fd1f1af23b49563d
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
<<<<<<< HEAD
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
=======
    if mode == "single" {
        // Prepare output files with compression
        let mut matrix_file = GzEncoder::new(File::create(format!("{}/matrix.mtx.gz", output_dir))?, Compression::default());
        let mut barcodes_file = GzEncoder::new(File::create(format!("{}/barcodes.tsv.gz", output_dir))?, Compression::default());
        let mut features_file = GzEncoder::new(File::create(format!("{}/features.tsv.gz", output_dir))?, Compression::default());
        let mut output_tsv = GzEncoder::new(File::create(format!("{}/junction_barcodes.tsv.gz", output_dir))?, Compression::default());

        // Write barcodes.tsv.gz
        debug!("Writing barcodes.tsv.gz");
        let barcode_list: Vec<_> = cell_barcodes.iter().sorted().collect();
        for barcode in &barcode_list {
            writeln!(barcodes_file, "{}", barcode)?;
        }

        // Write features.tsv.gz
        debug!("Writing features.tsv.gz");
        let feature_list: Vec<_> = junction_counts.keys().sorted().collect();
        for feature in &feature_list {
            writeln!(features_file, "{}", feature)?;
        }

        // Buffers to accumulate lines for matrix.mtx.gz and output.tsv.gz
        let mut matrix_buffer: Vec<String> = Vec::new();
        let mut tsv_buffer: Vec<String> = Vec::new();

        // Add the header lines to the matrix buffer
        matrix_buffer.push("%%MatrixMarket matrix coordinate integer general".to_string());
        matrix_buffer.push("%".to_string());
        matrix_buffer.push(format!(
            "{} {} {}",
            feature_list.len(),
            barcode_list.len(),
            junction_counts.values().map(|c| c.len()).sum::<usize>()
        ));

        // Add sparse matrix data and TSV data to the buffers
        debug!("Writing matrix.mtx.gz and junction_barcodes.tsv.gz");
        let barcode_map: HashMap<_, _> = barcode_list.iter().enumerate().map(|(i, b)| (b.as_str(), i)).collect();
        tsv_buffer.push("Feature\tBarcode\tCount".to_string());
        for (i, feature) in feature_list.iter().enumerate() {
            if let Some(cell_counts) = junction_counts.get(*feature) {
                for (barcode, count) in cell_counts {
                    if let Some(&j) = barcode_map.get(barcode.as_str()) {
                        matrix_buffer.push(format!("{} {} {}", i + 1, j + 1, count));
                        tsv_buffer.push(format!("{}\t{}\t{}", feature, barcode, count));
                    }
                }
            }
        }

        // Write the accumulated lines to the compressed output files
        for line in matrix_buffer {
            writeln!(matrix_file, "{}", line)?;
        }
        for line in tsv_buffer {
            writeln!(output_tsv, "{}", line)?;
        }

    } else if mode == "bulk" {
        let mut output_file = GzEncoder::new(File::create(format!("{}/junction.tsv.gz", output_dir))?, Compression::default());
        debug!("Writing junction.tsv.gz");
        writeln!(output_file, "Junction\tCount")?;
        for (junction, count) in junction_totals.iter().sorted() {
            writeln!(output_file, "{}\t{}", junction, count)?;
>>>>>>> 33a432e75c1c3be5acd1d209fd1f1af23b49563d
        }
    }

    info!("Finished processing");
    Ok(())
}
