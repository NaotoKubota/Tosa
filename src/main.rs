use clap::{Arg, Command};
use rust_htslib::bam::{self, Read};
use rust_htslib::bam::record::{Aux, Cigar};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{Write, BufWriter};
use log::{info, LevelFilter};
use env_logger;
use itertools::Itertools;
use flate2::write::GzEncoder;
use flate2::Compression;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Set up command-line arguments using clap
    let matches = Command::new("Tosa")
        .version("0.1.0")
        .author("NaotoKubota")
        .about("Counts junction reads from scRNA-seq data")
        .arg(Arg::new("bam_file")
            .required(true)
            .help("Path to the BAM file"))
        .arg(Arg::new("output_prefix")
            .required(true)
            .help("Output prefix for the .mtx, barcodes.tsv, features.tsv, and output.tsv files"))
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
        .arg(Arg::new("verbose")
            .short('v')
            .long("verbose")
            .action(clap::ArgAction::SetTrue)
            .help("Enable verbose output to print all arguments"))
        .get_matches();

    // Parse arguments
    let bam_file = matches.get_one::<String>("bam_file").unwrap();
    let output_prefix = matches.get_one::<String>("output_prefix").unwrap();
    let min_anchor_length = *matches.get_one::<i64>("anchor_length").unwrap();
    let min_intron_length = *matches.get_one::<i64>("min_intron_length").unwrap();
    let max_intron_length = *matches.get_one::<i64>("max_intron_length").unwrap();
    let verbose = matches.get_flag("verbose");

    // Initialize the logger with the appropriate level
    if verbose {
        env_logger::Builder::from_default_env()
            .filter(None, LevelFilter::Info)
            .init();
    } else {
        env_logger::Builder::from_default_env()
            .filter(None, LevelFilter::Warn)
            .init();
    }

    // Log all arguments if verbose is enabled
    info!("Running Tosa with the following arguments:");
    info!("BAM file: {}", bam_file);
    info!("Output prefix: {}", output_prefix);
    info!("Minimum anchor length: {}", min_anchor_length);
    info!("Minimum intron length: {}", min_intron_length);
    info!("Maximum intron length: {}", max_intron_length);

    // Count total reads in the BAM file in a preliminary pass
    let total_reads = {
        let mut bam_reader = bam::Reader::from_path(bam_file)?;
        bam_reader.records().count()
    };
    info!("Total number of reads: {}", total_reads);

    // Open the BAM file again for processing
    let mut bam_reader = bam::Reader::from_path(bam_file)?;

    // Get reference names (chromosome names)
    let header = bam_reader.header().to_owned();
    let reference_names: Vec<String> = header
        .target_names()
        .iter()
        .map(|name| String::from_utf8_lossy(name).to_string())
        .collect();

    // HashMaps to store counts by junction and cell barcode
    let mut junction_counts: HashMap<String, HashMap<String, u32>> = HashMap::new();
    let mut cell_barcodes: HashSet<String> = HashSet::new();
    let mut junctions: HashSet<String> = HashSet::new();

    // Counter for tracking the number of reads processed
    let mut read_count = 0;
    let mut last_percentage = 0;

    // Iterate over each read in the BAM file
    for result in bam_reader.records() {
        let record = result?;
        read_count += 1;

        // Calculate and log progress at each 1% increment
        let progress_percentage = (read_count * 100) / total_reads;
        if progress_percentage > last_percentage {
            info!("Progress: {}% ({} / {})", progress_percentage, read_count, total_reads);
            last_percentage = progress_percentage;
        }

        // Extract reference name (chromosome) and start position
        let ref_name = reference_names[record.tid() as usize].clone();
        let mut current_pos = record.pos(); // Start of the alignment

        // Extract Cell Barcode (CB) from tags
        let cell_barcode = match record.aux(b"CB") {
            Ok(Aux::String(cb_str)) => Some(cb_str.to_string()),
            _ => None,
        };

        // If a cell barcode is present, process each CIGAR element for junctions
        if let Some(cb_str) = cell_barcode {
            cell_barcodes.insert(cb_str.clone());
            let cigar_vec = record.cigar();  // Create a longer-lived binding for the cigar data
            let cigars: Vec<_> = cigar_vec.iter().collect();
            for i in 0..cigars.len() {
                if let Cigar::RefSkip(len) = cigars[i] {
                    // Check intron length constraints
                    let intron_length = *len as i64;
                    if intron_length < min_intron_length || intron_length > max_intron_length {
                        // Skip junctions outside the specified intron length range
                        current_pos += intron_length;
                        continue;
                    }

                    // Check for minimum anchor length on both sides
                    let has_left_anchor = if i > 0 {
                        match cigars[i - 1] {
                            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => *l as i64 >= min_anchor_length,
                            _ => false,
                        }
                    } else {
                        false
                    };

                    let has_right_anchor = if i < cigars.len() - 1 {
                        match cigars[i + 1] {
                            Cigar::Match(r) | Cigar::Equal(r) | Cigar::Diff(r) => *r as i64 >= min_anchor_length,
                            _ => false
                        }
                    } else {
                        false
                    };

                    if has_left_anchor && has_right_anchor {
                        // Found a valid junction; calculate the coordinates
                        let start = current_pos + 1; // BAM is 0-based, BED is 1-based
                        let end = start + intron_length;
                        let junction_coords = format!("{}:{}-{}", ref_name, start, end);

                        // Increment count for this junction and cell barcode
                        junctions.insert(junction_coords.clone());
                        let junction_entry = junction_counts.entry(junction_coords).or_insert_with(HashMap::new);
                        *junction_entry.entry(cb_str.clone()).or_insert(0) += 1;
                    }

                    // Update current position to the end of this junction
                    current_pos += intron_length;
                } else {
                    // Update current position for non-junction elements
                    match cigars[i] {
                        Cigar::Match(l) | Cigar::Ins(l) | Cigar::SoftClip(l) | Cigar::Del(l) => {
                            current_pos += *l as i64;
                        }
                        _ => {}
                    }
                }
            }
        }
    }

    // Final log statement for total reads processed
    info!("Total reads processed: {}", read_count);

    // Prepare output files
    let mut matrix_file = BufWriter::new(File::create(format!("{}_matrix.mtx", output_prefix))?);
    let mut barcodes_file = GzEncoder::new(File::create(format!("{}_barcodes.tsv.gz", output_prefix))?, Compression::default());
    let mut features_file = GzEncoder::new(File::create(format!("{}_features.tsv.gz", output_prefix))?, Compression::default());
    let mut flat_tsv = GzEncoder::new(File::create(format!("{}_flat.tsv.gz", output_prefix))?, Compression::default());

    // Write barcodes.tsv
    let barcode_list: Vec<_> = cell_barcodes.iter().sorted().collect();
    for barcode in &barcode_list {
        writeln!(barcodes_file, "{}", barcode)?;
    }

    // Write features.tsv
    let feature_list: Vec<_> = junctions.iter().sorted().collect();
    for feature in &feature_list {
        writeln!(features_file, "{}", feature)?;
    }

    // Write matrix.mtx header
    writeln!(matrix_file, "%%MatrixMarket matrix coordinate integer general")?;
    writeln!(matrix_file, "%")?;
    writeln!(matrix_file, "{} {} {}", feature_list.len(), barcode_list.len(), junction_counts.values().map(|c| c.len()).sum::<usize>())?;

    // Write the sparse matrix data and output.tsv
    for (i, feature) in feature_list.iter().enumerate() {
        if let Some(cell_counts) = junction_counts.get(*feature) {
            for (barcode, count) in cell_counts {
                if let Some(j) = barcode_list.iter().position(|b| *b == barcode) {
                    writeln!(matrix_file, "{} {} {}", i + 1, j + 1, count)?;
                    writeln!(flat_tsv, "{}\t{}\t{}", feature, barcode, count)?;
                }
            }
        }
    }

    Ok(())
}
