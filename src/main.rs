use clap::{Arg, Command};
use rust_htslib::bam::{self, Read};
use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::record::{Aux, Cigar};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{Write};
use log::{info, debug, LevelFilter};
use env_logger;
use itertools::Itertools;
use flate2::write::GzEncoder;
use flate2::Compression;

mod data_loader;
mod junction;

fn main() -> Result<(), Box<dyn std::error::Error>> {
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
        env_logger::Builder::from_default_env()
            .filter(None, LevelFilter::Debug)
            .init();
    } else {
        env_logger::Builder::from_default_env()
            .filter(None, LevelFilter::Info)
            .init();
    }

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

    // Count total mapped reads in the BAM file
    let mut bam_index_reader = IndexedReader::from_path(bam_file)?;
    let stats = bam_index_reader.index_stats()?;
    debug!("stats: {:?}", stats);
    // Sum the mapped reads from all targets
    let total_mapped_reads: u64 = stats.iter().map(|(_, _, mapped, _)| mapped).sum();
    info!("Total number of reads: {}", total_mapped_reads);

    // Open the BAM file again for processing
    let mut bam_reader = bam::Reader::from_path(bam_file)?;

    // Get reference names (chromosome names)
    let header = bam_reader.header().to_owned();
    let reference_names: Vec<String> = header
        .target_names()
        .iter()
        .map(|name| String::from_utf8_lossy(name).to_string())
        .collect();

    // HashMaps to store counts by junction and optionally by cell barcode
    let mut junction_counts: HashMap<String, HashMap<String, u32>> = HashMap::new();
    let mut junction_totals: HashMap<String, u32> = HashMap::new();
    let mut cell_barcodes: HashSet<String> = HashSet::new();

    // Counter for tracking the number of reads processed
    let mut read_count = 0;
    let mut last_percentage = 0;

    // HashSet to store supported junctions and HashMap to store buffered reads
    let mut supported_junctions: HashSet<String> = HashSet::new();
    let mut buffered_reads: HashMap<String, Vec<(Option<String>, i64)>> = HashMap::new();

    // HashMap to store processed reads by junction
    let mut processed_reads: HashMap<String, HashSet<String>> = HashMap::new();

    // Iterate over each read in the BAM file
    for result in bam_reader.records() {
        let record = result?;
        read_count += 1;

        // Calculate and log progress at each 1% increment
        let progress_percentage = (read_count * 100) / total_mapped_reads;
        if progress_percentage > last_percentage {
            info!("Progress: {}% ({} / {})", progress_percentage, read_count, total_mapped_reads);
            last_percentage = progress_percentage;
        }

        // Skip read if NH tag exceeds max_loci
        if let Ok(Aux::U8(nh)) = record.aux(b"NH") {
            if nh > max_loci as u8 {
                // debug!("Skipping read {} with NH > max_loci ({})", std::str::from_utf8(record.qname()).unwrap(), nh);
                continue; // Skip this read
            }
        } else if let Ok(Aux::I32(nh)) = record.aux(b"NH") {
            if nh > max_loci as i32 {
                // debug!("Skipping read {} with NH > max_loci ({})", std::str::from_utf8(record.qname()).unwrap(), nh);
                continue; // Skip this read
            }
        }

        // Extract reference name (chromosome) and start position
        let ref_name = reference_names[record.tid() as usize].clone();
        let mut current_pos = record.pos(); // Start of the alignment

        // Extract Cell Barcode (CB) from tags if in single mode
        let cell_barcode = if mode == "single" {
            match record.aux(b"CB") {
                Ok(Aux::String(cb_str)) => Some(cb_str.to_string()),
                _ => None,
            }
        } else {
            None
        };

        // Skip read if its barcode is not in the list of interest
        if let Some(cb) = &cell_barcode {
            if cell_barcode_file.is_some() && !cell_barcodes_of_interest.is_empty() && !cell_barcodes_of_interest.contains(cb) {
            continue;
            }
        }

        // If a cell barcode is present (for single mode), or always process for bulk mode
        if mode == "bulk" || cell_barcode.is_some() {
            if let Some(cb_str) = &cell_barcode {
                cell_barcodes.insert(cb_str.clone());
            }

            let cigar_vec = record.cigar(); // Create a longer-lived binding for the cigar data
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

                    // Calculate left anchor length by accumulating lengths before the RefSkip
                    let mut left_anchor_length = 0;
                    let mut j = i; // Start from the current CIGAR index
                    while j > 0 {
                        j -= 1; // Move to the previous CIGAR element
                        match cigars[j] {
                            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                                left_anchor_length += *l as i64;
                                if left_anchor_length >= min_anchor_length {
                                    break; // Stop if the threshold is met
                                }
                            }
                            Cigar::RefSkip(_) => continue, // Skip RefSkip and keep checking alignment elements
                            _ => break, // Stop accumulating for other operations
                        }
                    }
                    let has_left_anchor = left_anchor_length >= min_anchor_length;

                    // Calculate right anchor length by accumulating lengths after the RefSkip
                    let mut right_anchor_length = 0;
                    let mut k = i + 1; // Start from the next CIGAR index
                    while k < cigars.len() {
                        match cigars[k] {
                            Cigar::Match(r) | Cigar::Equal(r) | Cigar::Diff(r) => {
                                right_anchor_length += *r as i64;
                                if right_anchor_length >= min_anchor_length {
                                    break; // Stop if the threshold is met
                                }
                            }
                            Cigar::RefSkip(_) => { k += 1; continue; } // Skip RefSkip and keep checking alignment elements
                            _ => break, // Stop accumulating for other operations
                        }
                        k += 1; // Move to the next CIGAR element
                    }
                    let has_right_anchor = right_anchor_length >= min_anchor_length;

                    let start = current_pos;
                    let end = start + intron_length + 1;
                    let junction_coords = format!("{}:{}-{}", ref_name, start, end);

                    if has_left_anchor && has_right_anchor {
                        // Mark as supported and process buffered reads
                        supported_junctions.insert(junction_coords.clone());
                        if let Some(buffered) = buffered_reads.remove(&junction_coords) {
                            for (buffered_cb, _buffered_pos) in buffered {
                                junction::process_junction(
                                    &junction_coords,
                                    buffered_cb.as_ref(),
                                    &mut junction_counts,
                                    &mut junction_totals,
                                    &mut processed_reads, // Pass the processed reads map
                                    std::str::from_utf8(record.qname()).unwrap(), // Pass read name
                                    &mode,
                                );
                            }
                        }
                    }

                    // Process or buffer the current read
                    if supported_junctions.contains(&junction_coords) {
                        junction::process_junction(
                            &junction_coords,
                            cell_barcode.as_ref(),
                            &mut junction_counts,
                            &mut junction_totals,
                            &mut processed_reads, // Pass the processed reads map
                            std::str::from_utf8(record.qname()).unwrap(), // Pass read name
                            &mode,
                        );
                    } else {
                        buffered_reads
                            .entry(junction_coords.clone())
                            .or_insert_with(Vec::new)
                            .push((cell_barcode.clone(), current_pos));
                    }
                    current_pos += intron_length;
                } else if let Cigar::SoftClip(_len) = cigars[i] {
                    continue;
                } else {
                    current_pos += match cigars[i] {
                        Cigar::Match(l) | Cigar::Ins(l) | Cigar::Del(l) => *l as i64,
                        _ => 0,
                    };
                }
            }
        }
    }

    // Write results based on mode
    info!("Writing output files");
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
        }
    }

    info!("Finished processing");
    Ok(())
}
