//! BAM file processing: read iteration, CIGAR parsing, junction extraction, and boundary counting.

use rust_htslib::bam::{self, Read};
use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::record::{Aux, Cigar};
use std::collections::{HashMap, HashSet};
use log::{info, debug};

use crate::types::{RunConfig, StrandMode, Strand};
use crate::junction;
use crate::boundary::{self, BoundaryIndex};

/// Results from processing a BAM file.
pub struct ProcessingResult {
    /// Per-cell junction counts: junction_key -> (barcode -> count). Used in single mode.
    pub junction_counts: HashMap<String, HashMap<String, u32>>,
    /// Total junction counts: junction_key -> count. Used in bulk mode.
    pub junction_totals: HashMap<String, u32>,
    /// Strand assigned to each junction key.
    pub junction_strands: HashMap<String, Strand>,
    /// All observed cell barcodes.
    pub cell_barcodes: HashSet<String>,
    /// Per-cell boundary counts: boundary_key -> (barcode -> count). Used in single mode.
    pub boundary_counts: HashMap<String, HashMap<String, u32>>,
    /// Total boundary counts: boundary_key -> count. Used in bulk mode.
    pub boundary_totals: HashMap<String, u32>,
    /// Boundary type for each boundary key.
    pub boundary_types: HashMap<String, crate::types::BoundaryType>,
    /// Strand assigned to each boundary key.
    pub boundary_strands: HashMap<String, Strand>,
}

/// Count total mapped reads using the BAM index.
pub fn count_total_reads(bam_file: &str) -> Result<u64, Box<dyn std::error::Error>> {
    let mut bam_index_reader = IndexedReader::from_path(bam_file)?;
    let stats = bam_index_reader.index_stats()?;
    debug!("stats: {:?}", stats);
    let total_mapped_reads: u64 = stats.iter().map(|(_, _, mapped, _)| mapped).sum();
    Ok(total_mapped_reads)
}

/// Determine strand of a read based on the strand mode.
pub fn determine_strand(record: &bam::Record, strand_mode: &StrandMode) -> Strand {
    match strand_mode {
        StrandMode::Unstranded => {
            Strand::Unknown
        }
        StrandMode::XS => {
            match record.aux(b"XS") {
                Ok(Aux::Char(c)) => match c {
                    b'+' => Strand::Plus,
                    b'-' => Strand::Minus,
                    _ => Strand::Unknown,
                },
                _ => Strand::Unknown,
            }
        }
        StrandMode::RF => {
            // First-strand: read1+reverse => +, read1+forward => -
            //               read2+reverse => -, read2+forward => +
            let is_reverse = record.is_reverse();
            if record.is_first_in_template() {
                if is_reverse { Strand::Plus } else { Strand::Minus }
            } else if record.is_last_in_template() {
                if is_reverse { Strand::Minus } else { Strand::Plus }
            } else {
                // Single-end read in RF mode: reverse => +, forward => -
                if is_reverse { Strand::Plus } else { Strand::Minus }
            }
        }
        StrandMode::FR => {
            // Second-strand: read1+forward => +, read1+reverse => -
            //                read2+forward => -, read2+reverse => +
            let is_reverse = record.is_reverse();
            if record.is_first_in_template() {
                if is_reverse { Strand::Minus } else { Strand::Plus }
            } else if record.is_last_in_template() {
                if is_reverse { Strand::Plus } else { Strand::Minus }
            } else {
                // Single-end read in FR mode: forward => +, reverse => -
                if is_reverse { Strand::Minus } else { Strand::Plus }
            }
        }
    }
}

/// Extract aligned segments from a CIGAR string.
/// Returns a list of (start, end) reference coordinate intervals for M/=/X operations.
pub fn extract_aligned_segments(record: &bam::Record) -> Vec<(i64, i64)> {
    let mut segments = Vec::new();
    let mut current_pos = record.pos();

    let cigar_view = record.cigar();
    for op in cigar_view.iter() {
        match op {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                let len = *l as i64;
                segments.push((current_pos, current_pos + len));
                current_pos += len;
            }
            Cigar::Del(l) | Cigar::RefSkip(l) => {
                current_pos += *l as i64;
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {
                // These do not consume reference bases
            }
        }
    }
    segments
}

/// Process all records in a BAM file, extracting junction and boundary counts.
pub fn process_bam_records(
    config: &RunConfig,
    cell_barcodes_of_interest: &HashSet<String>,
    boundary_index: Option<&BoundaryIndex>,
) -> Result<ProcessingResult, Box<dyn std::error::Error>> {
    let total_mapped_reads = count_total_reads(&config.bam_file)?;
    info!("Total number of reads: {}", total_mapped_reads);

    let mut bam_reader = bam::Reader::from_path(&config.bam_file)?;

    // Get reference names (chromosome names)
    let header = bam_reader.header().to_owned();
    let reference_names: Vec<String> = header
        .target_names()
        .iter()
        .map(|name| String::from_utf8_lossy(name).to_string())
        .collect();

    // Junction state
    let mut junction_counts: HashMap<String, HashMap<String, u32>> = HashMap::new();
    let mut junction_totals: HashMap<String, u32> = HashMap::new();
    let mut junction_strands: HashMap<String, Strand> = HashMap::new();
    let mut cell_barcodes: HashSet<String> = HashSet::new();
    let mut supported_junctions: HashSet<String> = HashSet::new();
    let mut buffered_reads: HashMap<String, Vec<(Option<String>, i64, Strand)>> = HashMap::new();
    let mut processed_reads: HashMap<String, HashSet<String>> = HashMap::new();

    // Boundary state
    let mut boundary_counts: HashMap<String, HashMap<String, u32>> = HashMap::new();
    let mut boundary_totals: HashMap<String, u32> = HashMap::new();
    let mut boundary_types: HashMap<String, crate::types::BoundaryType> = HashMap::new();
    let mut boundary_strands: HashMap<String, Strand> = HashMap::new();
    let mut processed_boundary_reads: HashMap<String, HashSet<String>> = HashMap::new();

    // Progress tracking
    let mut read_count: u64 = 0;
    let mut last_percentage: u64 = 0;

    for result in bam_reader.records() {
        let record = result?;
        read_count += 1;

        // Progress logging
        if total_mapped_reads > 0 {
            let progress_percentage = (read_count * 100) / total_mapped_reads;
            if progress_percentage > last_percentage {
                info!("Progress: {}% ({} / {})", progress_percentage, read_count, total_mapped_reads);
                last_percentage = progress_percentage;
            }
        }

        // Skip read if NH tag exceeds max_loci
        if let Ok(Aux::U8(nh)) = record.aux(b"NH") {
            if nh > config.max_loci as u8 {
                continue;
            }
        } else if let Ok(Aux::I32(nh)) = record.aux(b"NH") {
            if nh > config.max_loci as i32 {
                continue;
            }
        }

        // Extract reference name (chromosome) and start position
        let tid = record.tid();
        if tid < 0 {
            continue; // Unmapped read
        }
        let ref_name = reference_names[tid as usize].clone();
        let mut current_pos = record.pos();

        // Extract Cell Barcode (CB) from tags if in single mode
        let cell_barcode = if config.mode == "single" {
            match record.aux(b"CB") {
                Ok(Aux::String(cb_str)) => Some(cb_str.to_string()),
                _ => None,
            }
        } else {
            None
        };

        // Skip read if its barcode is not in the list of interest
        if let Some(cb) = &cell_barcode {
            if config.cell_barcode_file.is_some()
                && !cell_barcodes_of_interest.is_empty()
                && !cell_barcodes_of_interest.contains(cb)
            {
                continue;
            }
        }

        // If a cell barcode is present (for single mode), or always process for bulk mode
        if config.mode == "bulk" || cell_barcode.is_some() {
            if let Some(cb_str) = &cell_barcode {
                cell_barcodes.insert(cb_str.clone());
            }

            // Determine strand for this read
            let strand = determine_strand(&record, &config.strand_mode);
            let read_name = std::str::from_utf8(record.qname()).unwrap_or("unknown");

            // --- Junction extraction from CIGAR ---
            let cigar_view = record.cigar();
            let cigars: Vec<_> = cigar_view.iter().collect();
            for i in 0..cigars.len() {
                if let Cigar::RefSkip(len) = cigars[i] {
                    let intron_length = *len as i64;
                    if intron_length < config.min_intron_length || intron_length > config.max_intron_length {
                        current_pos += intron_length;
                        continue;
                    }

                    // Calculate left anchor length
                    let mut left_anchor_length: i64 = 0;
                    let mut j = i;
                    while j > 0 {
                        j -= 1;
                        match cigars[j] {
                            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                                left_anchor_length += *l as i64;
                                if left_anchor_length >= config.min_anchor_length {
                                    break;
                                }
                            }
                            Cigar::RefSkip(_) => continue,
                            _ => break,
                        }
                    }
                    let has_left_anchor = left_anchor_length >= config.min_anchor_length;

                    // Calculate right anchor length
                    let mut right_anchor_length: i64 = 0;
                    let mut k = i + 1;
                    while k < cigars.len() {
                        match cigars[k] {
                            Cigar::Match(r) | Cigar::Equal(r) | Cigar::Diff(r) => {
                                right_anchor_length += *r as i64;
                                if right_anchor_length >= config.min_anchor_length {
                                    break;
                                }
                            }
                            Cigar::RefSkip(_) => { k += 1; continue; }
                            _ => break,
                        }
                        k += 1;
                    }
                    let has_right_anchor = right_anchor_length >= config.min_anchor_length;

                    // Use 1-based inclusive coordinates for junction IDs
                    let start = current_pos + 1;          // 1-based intron start
                    let end = current_pos + intron_length; // 1-based intron end
                    let junction_coords = format!("{}:{}-{}", ref_name, start, end);

                    if has_left_anchor && has_right_anchor {
                        // Mark as supported and process buffered reads
                        supported_junctions.insert(junction_coords.clone());
                        junction_strands.entry(junction_coords.clone()).or_insert(strand);

                        if let Some(buffered) = buffered_reads.remove(&junction_coords) {
                            for (buffered_cb, _buffered_pos, buffered_strand) in buffered {
                                junction::process_junction(
                                    &junction_coords,
                                    buffered_cb.as_ref(),
                                    buffered_strand,
                                    &mut junction_counts,
                                    &mut junction_totals,
                                    &mut junction_strands,
                                    &mut processed_reads,
                                    read_name,
                                    &config.mode,
                                );
                            }
                        }
                    }

                    // Process or buffer the current read
                    if supported_junctions.contains(&junction_coords) {
                        junction::process_junction(
                            &junction_coords,
                            cell_barcode.as_ref(),
                            strand,
                            &mut junction_counts,
                            &mut junction_totals,
                            &mut junction_strands,
                            &mut processed_reads,
                            read_name,
                            &config.mode,
                        );
                    } else {
                        buffered_reads
                            .entry(junction_coords.clone())
                            .or_default()
                            .push((cell_barcode.clone(), current_pos, strand));
                    }
                    current_pos += intron_length;
                } else if let Cigar::SoftClip(_) = cigars[i] {
                    // SoftClip does not consume reference bases
                    continue;
                } else {
                    // FIX: Ins does NOT consume reference bases
                    current_pos += match cigars[i] {
                        Cigar::Match(l) | Cigar::Del(l) => *l as i64,
                        _ => 0,
                    };
                }
            }

            // --- Boundary counting ---
            if let Some(bi) = boundary_index {
                let segments = extract_aligned_segments(&record);
                boundary::count_boundaries(
                    &ref_name,
                    &segments,
                    bi,
                    cell_barcode.as_ref(),
                    strand,
                    &mut boundary_counts,
                    &mut boundary_totals,
                    &mut boundary_types,
                    &mut boundary_strands,
                    &mut processed_boundary_reads,
                    read_name,
                    &config.mode,
                );
            }
        }
    }

    Ok(ProcessingResult {
        junction_counts,
        junction_totals,
        junction_strands,
        cell_barcodes,
        boundary_counts,
        boundary_totals,
        boundary_types,
        boundary_strands,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_determine_strand_unstranded_no_xs() {
        // Without a real BAM record with XS tag, test the logic path
        // For unit testing determine_strand, we'd need mock records.
        // This is a placeholder — full tests use generated BAM files.
        let mode = StrandMode::Unstranded;
        assert_eq!(format!("{}", mode), "unstranded");
    }

    #[test]
    fn test_strand_mode_display() {
        assert_eq!(format!("{}", StrandMode::RF), "RF");
        assert_eq!(format!("{}", StrandMode::FR), "FR");
        assert_eq!(format!("{}", StrandMode::XS), "XS");
        assert_eq!(format!("{}", StrandMode::Unstranded), "unstranded");
    }
}
