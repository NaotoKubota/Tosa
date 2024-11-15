use clap::{Arg, Command};
use rust_htslib::bam::{self, Read};
use rust_htslib::bam::record::{Aux, Cigar};
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use log::{info, LevelFilter};
use env_logger;
use itertools::Itertools;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Set up command-line arguments using clap
    let matches = Command::new("Tosa")
        .version("1.0")
        .author("Your Name")
        .about("Counts junction reads with a minimum anchor length requirement")
        .arg(Arg::new("bam_file")
            .required(true)
            .help("Path to the BAM file"))
        .arg(Arg::new("output_file")
            .required(true)
            .help("Path to the output file"))
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
    let output_file = matches.get_one::<String>("output_file").unwrap();
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
    info!("Output file: {}", output_file);
    info!("Minimum anchor length: {}", min_anchor_length);
    info!("Minimum intron length: {}", min_intron_length);
    info!("Maximum intron length: {}", max_intron_length);

    // Open the BAM file
    let mut bam_reader = bam::Reader::from_path(bam_file)?;

    // Get reference names (chromosome names)
    let header = bam_reader.header().to_owned();
    let reference_names: Vec<String> = header
        .target_names()
        .iter()
        .map(|name| String::from_utf8_lossy(name).to_string())
        .collect();

    // HashMap to store counts by junction and cell barcode
    let mut junction_counts: HashMap<String, HashMap<String, u32>> = HashMap::new();

    // Counter for tracking the number of reads processed
    let mut read_count = 0;

    // Iterate over each read in the BAM file
    for result in bam_reader.records() {
        let record = result?;
        read_count += 1;

        // Log every 10,000 reads processed
        if read_count % 10000 == 0 {
            info!("Processed {} reads", read_count);
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

    // Open the output file for writing
    let mut output = File::create(output_file)?;

    // Collect junction entries, parse coordinates, and sort by chromosomal position
    let mut sorted_junctions: Vec<_> = junction_counts.iter().collect();
    sorted_junctions.sort_by_key(|(junction, _)| {
        let parts: Vec<&str> = junction.split([':', '-']).collect();
        let chrom = parts[0].to_string();
        let start: i64 = parts[1].parse().unwrap_or(0);
        let end: i64 = parts[2].parse().unwrap_or(0);
        (chrom, start, end)
    });

    // Write the sorted read counts per junction, with breakdown per cell barcode sorted by read count
    writeln!(output, "Junction\tCount\tCell Barcode")?;
    for (junction, cell_counts) in sorted_junctions {
        // Calculate total count for this junction
        let total_count: u32 = cell_counts.values().sum();

        // Sort cell barcodes by read count in descending order and format as comma-separated list
        let cell_barcode_info: String = cell_counts
            .iter()
            .sorted_by(|a, b| b.1.cmp(a.1))
            .map(|(cell_barcode, count)| format!("{}:{}", cell_barcode, count))
            .collect::<Vec<String>>()
            .join(",");

        writeln!(output, "{}\t{}\t{}", junction, total_count, cell_barcode_info)?;
    }

    Ok(())
}
