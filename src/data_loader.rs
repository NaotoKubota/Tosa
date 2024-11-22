// Modules for data loading
use std::collections::{HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use rust_htslib::bam::IndexedReader;

// Function to count total mapped reads in the BAM file
pub fn count_total_mapped_reads(bam_file: &str) -> Result<u64, Box<dyn std::error::Error>> {
	let mut bam_index_reader = IndexedReader::from_path(bam_file)?;
	let stats = bam_index_reader.index_stats()?;
	// Sum the mapped reads from all targets
	let total_mapped_reads: u64 = stats.iter().map(|(_, _, mapped, _)| mapped).sum();
	Ok(total_mapped_reads)
}

// Function to load the cell barcodes
pub fn load_cell_barcodes(file_path: Option<&String>) -> Result<HashSet<String>, Box<dyn std::error::Error>> {
    let mut barcodes = HashSet::new();
    if let Some(path) = file_path {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let barcode = line?.trim().to_string();
            barcodes.insert(barcode);
        }
    }
    Ok(barcodes)
}

// Function to load the intron list
// pub fn load_introns(file_path: Option<&String>) -> Result<HashSet<(String, i64, i64)>, Box<dyn std::error::Error>> {
//     let mut introns = HashSet::new();
//     if let Some(path) = file_path {
//         let file = File::open(path)?;
//         let reader = BufReader::new(file);
//         for line in reader.lines() {
//             let line = line?;
//             let parts: Vec<&str> = line.split(':').collect();
//             if parts.len() == 2 {
//                 let chrom = parts[0].to_string();
//                 let range: Vec<&str> = parts[1].split('-').collect();
//                 if range.len() == 2 {
//                     let start = range[0].parse::<i64>()?;
//                     let end = range[1].parse::<i64>()?;
//                     introns.insert((chrom, start, end));
//                 }
//             }
//         }
//     }
//     Ok(introns)
// }
