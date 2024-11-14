use rust_htslib::bam::{self, Read};
use rust_htslib::bam::record::{Aux, Cigar};
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::Write;
use log::{info, LevelFilter};
use env_logger;
use itertools::Itertools;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Initialize the logger
    env_logger::Builder::from_default_env()
        .filter(None, LevelFilter::Info)
        .init();

    // Get command-line arguments
    let args: Vec<String> = env::args().collect();

    // Check if both the BAM file path and output file path are provided
    if args.len() < 3 {
        eprintln!("Usage: {} <bam_file> <output_file>", args[0]);
        std::process::exit(1);
    }

    // Use the BAM file and output file paths from the command line
    let bam_file = &args[1];
    let output_file = &args[2];

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

        // Log every 10000 reads processed
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
            for cigar in record.cigar().iter() {
                match cigar {
                    Cigar::RefSkip(len) => {
                        // Found a junction; calculate the coordinates
                        let start = current_pos + 1; // BAM is 0-based, BED is 1-based
                        let end = start + (*len as i64);
                        let junction_coords = format!("{}:{}-{}", ref_name, start, end);

                        // Increment count for this junction and cell barcode
                        let junction_entry = junction_counts.entry(junction_coords).or_insert_with(HashMap::new);
                        *junction_entry.entry(cb_str.clone()).or_insert(0) += 1;

                        // Update current position to the end of this junction
                        current_pos = end - 1;
                    }
                    Cigar::Match(len) | Cigar::Ins(len) | Cigar::SoftClip(len) | Cigar::Del(len) => {
                        current_pos += *len as i64;
                    }
                    _ => {}
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
