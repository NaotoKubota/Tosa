//! Junction read counting logic.

use std::collections::{HashMap, HashSet};
use crate::types::Strand;

/// Process a junction read, incrementing counts and tracking duplicates.
///
/// Junctions are keyed by `(junction_coords, strand)` internally.
/// Duplicate reads (same read name for same junction) are counted only once.
#[allow(clippy::too_many_arguments)]
pub fn process_junction(
    junction_coords: &str,
    cell_barcode: Option<&String>,
    strand: Strand,
    junction_counts: &mut HashMap<String, HashMap<String, u32>>,
    junction_totals: &mut HashMap<String, u32>,
    junction_strands: &mut HashMap<String, Strand>,
    processed_reads: &mut HashMap<String, HashSet<String>>,
    read_name: &str,
    mode: &str,
) {
    // Composite key includes strand to distinguish same-position different-strand junctions
    let key = format!("{}:{}", junction_coords, strand);

    // Check if the read was already processed for this junction
    if let Some(reads) = processed_reads.get_mut(&key) {
        if reads.contains(read_name) {
            return; // Skip counting
        }
        reads.insert(read_name.to_string());
    } else {
        let mut reads_set = HashSet::new();
        reads_set.insert(read_name.to_string());
        processed_reads.insert(key.clone(), reads_set);
    }

    // Record strand for this junction key
    junction_strands.entry(key.clone()).or_insert(strand);

    // Count the read for the junction
    if mode == "single" {
        if let Some(cb_str) = cell_barcode {
            let junction_entry = junction_counts
                .entry(key)
                .or_default();
            *junction_entry.entry(cb_str.clone()).or_insert(0) += 1;
        }
    } else {
        *junction_totals
            .entry(key)
            .or_insert(0) += 1;
    }
}

/// Parse a composite junction key "chr:start-end:strand" back into parts.
/// Returns (junction_coords, strand_str).
pub fn parse_junction_key(key: &str) -> (&str, &str) {
    // Key format: "chr:start-end:strand" where strand is +, -, or .
    // We need to split from the last ':'
    if let Some(last_colon) = key.rfind(':') {
        let coords = &key[..last_colon];
        let strand = &key[last_colon + 1..];
        (coords, strand)
    } else {
        (key, ".")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_process_junction_bulk() {
        let mut junction_counts = HashMap::new();
        let mut junction_totals = HashMap::new();
        let mut junction_strands = HashMap::new();
        let mut processed_reads = HashMap::new();

        process_junction(
            "chr1:100-200",
            None,
            Strand::Plus,
            &mut junction_counts,
            &mut junction_totals,
            &mut junction_strands,
            &mut processed_reads,
            "read1",
            "bulk",
        );

        assert_eq!(junction_totals.get("chr1:100-200:+"), Some(&1));
        assert_eq!(*junction_strands.get("chr1:100-200:+").unwrap(), Strand::Plus);
    }

    #[test]
    fn test_process_junction_dedup() {
        let mut junction_counts = HashMap::new();
        let mut junction_totals = HashMap::new();
        let mut junction_strands = HashMap::new();
        let mut processed_reads = HashMap::new();

        // Process same read twice for same junction
        for _ in 0..2 {
            process_junction(
                "chr1:100-200",
                None,
                Strand::Plus,
                &mut junction_counts,
                &mut junction_totals,
                &mut junction_strands,
                &mut processed_reads,
                "read1",
                "bulk",
            );
        }

        // Should only be counted once
        assert_eq!(junction_totals.get("chr1:100-200:+"), Some(&1));
    }

    #[test]
    fn test_process_junction_different_strands() {
        let mut junction_counts = HashMap::new();
        let mut junction_totals = HashMap::new();
        let mut junction_strands = HashMap::new();
        let mut processed_reads = HashMap::new();

        process_junction(
            "chr1:100-200",
            None,
            Strand::Plus,
            &mut junction_counts,
            &mut junction_totals,
            &mut junction_strands,
            &mut processed_reads,
            "read1",
            "bulk",
        );

        process_junction(
            "chr1:100-200",
            None,
            Strand::Minus,
            &mut junction_counts,
            &mut junction_totals,
            &mut junction_strands,
            &mut processed_reads,
            "read2",
            "bulk",
        );

        // Same coords, different strand → separate counts
        assert_eq!(junction_totals.get("chr1:100-200:+"), Some(&1));
        assert_eq!(junction_totals.get("chr1:100-200:-"), Some(&1));
    }

    #[test]
    fn test_process_junction_single_mode() {
        let mut junction_counts = HashMap::new();
        let mut junction_totals = HashMap::new();
        let mut junction_strands = HashMap::new();
        let mut processed_reads = HashMap::new();

        let barcode = "ACGT-1".to_string();
        process_junction(
            "chr1:100-200",
            Some(&barcode),
            Strand::Unknown,
            &mut junction_counts,
            &mut junction_totals,
            &mut junction_strands,
            &mut processed_reads,
            "read1",
            "single",
        );

        let key = "chr1:100-200:.";
        assert_eq!(
            *junction_counts.get(key).unwrap().get("ACGT-1").unwrap(),
            1
        );
    }

    #[test]
    fn test_parse_junction_key() {
        let (coords, strand) = parse_junction_key("chr1:100-200:+");
        assert_eq!(coords, "chr1:100-200");
        assert_eq!(strand, "+");

        let (coords, strand) = parse_junction_key("chrX:5000-6000:.");
        assert_eq!(coords, "chrX:5000-6000");
        assert_eq!(strand, ".");
    }
}
