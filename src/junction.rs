//! Junction read counting logic.

use std::collections::{HashMap, HashSet};
#[allow(unused_imports)] // Used by tests via `use super::*`
use crate::types::{JunctionKey, Mode, Strand, hash_barcode_umi};

/// Process a junction read, incrementing counts and tracking duplicates.
///
/// Duplicate reads (same read-name hash for same junction) are counted only once.
/// In single-cell mode, when a UMI is available, reads with the same
/// (barcode, UMI) pair for the same junction are deduplicated so that
/// each unique molecule is counted only once.
pub fn process_junction(
    key: JunctionKey,
    cell_barcode: Option<&String>,
    umi: Option<&String>,
    junction_counts: &mut HashMap<JunctionKey, HashMap<String, u32>>,
    junction_totals: &mut HashMap<JunctionKey, u32>,
    processed_reads: &mut HashMap<JunctionKey, HashSet<u64>>,
    processed_umis: &mut HashMap<JunctionKey, HashSet<u64>>,
    read_name_hash: u64,
    mode: Mode,
) {
    // Check if the read was already processed for this junction (via u64 hash)
    // This deduplicates paired-end reads (same QNAME) mapping to the same junction.
    let reads = processed_reads.entry(key).or_default();
    if !reads.insert(read_name_hash) {
        return; // Already counted
    }

    // Count the read for the junction
    match mode {
        Mode::Single => {
            if let Some(cb_str) = cell_barcode {
                // UMI deduplication: if a UMI is present, ensure that
                // (barcode, UMI) is unique per junction before counting.
                if let Some(umi_str) = umi {
                    let umi_hash = hash_barcode_umi(cb_str, umi_str);
                    let umis = processed_umis.entry(key).or_default();
                    if !umis.insert(umi_hash) {
                        return; // Same barcode+UMI already counted for this junction
                    }
                }
                *junction_counts
                    .entry(key)
                    .or_default()
                    .entry(cb_str.clone())
                    .or_insert(0) += 1;
            }
        }
        Mode::Bulk => {
            *junction_totals.entry(key).or_insert(0) += 1;
        }
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
        let mut processed_reads = HashMap::new();
        let mut processed_umis = HashMap::new();

        let key = JunctionKey { tid: 0, start: 100, end: 200, strand: Strand::Plus };

        process_junction(
            key,
            None,
            None,
            &mut junction_counts,
            &mut junction_totals,
            &mut processed_reads,
            &mut processed_umis,
            12345, // read_name_hash
            Mode::Bulk,
        );

        assert_eq!(junction_totals.get(&key), Some(&1));
    }

    #[test]
    fn test_process_junction_dedup() {
        let mut junction_counts = HashMap::new();
        let mut junction_totals = HashMap::new();
        let mut processed_reads = HashMap::new();
        let mut processed_umis = HashMap::new();

        let key = JunctionKey { tid: 0, start: 100, end: 200, strand: Strand::Plus };

        // Process same read hash twice for same junction
        for _ in 0..2 {
            process_junction(
                key,
                None,
                None,
                &mut junction_counts,
                &mut junction_totals,
                &mut processed_reads,
                &mut processed_umis,
                12345,
                Mode::Bulk,
            );
        }

        // Should only be counted once
        assert_eq!(junction_totals.get(&key), Some(&1));
    }

    #[test]
    fn test_process_junction_different_strands() {
        let mut junction_counts = HashMap::new();
        let mut junction_totals = HashMap::new();
        let mut processed_reads = HashMap::new();
        let mut processed_umis = HashMap::new();

        let key_plus = JunctionKey { tid: 0, start: 100, end: 200, strand: Strand::Plus };
        let key_minus = JunctionKey { tid: 0, start: 100, end: 200, strand: Strand::Minus };

        process_junction(
            key_plus,
            None,
            None,
            &mut junction_counts,
            &mut junction_totals,
            &mut processed_reads,
            &mut processed_umis,
            11111,
            Mode::Bulk,
        );

        process_junction(
            key_minus,
            None,
            None,
            &mut junction_counts,
            &mut junction_totals,
            &mut processed_reads,
            &mut processed_umis,
            22222,
            Mode::Bulk,
        );

        // Same coords, different strand → separate counts
        assert_eq!(junction_totals.get(&key_plus), Some(&1));
        assert_eq!(junction_totals.get(&key_minus), Some(&1));
    }

    #[test]
    fn test_process_junction_single_mode() {
        let mut junction_counts = HashMap::new();
        let mut junction_totals = HashMap::new();
        let mut processed_reads = HashMap::new();
        let mut processed_umis = HashMap::new();

        let key = JunctionKey { tid: 0, start: 100, end: 200, strand: Strand::Unknown };
        let barcode = "ACGT-1".to_string();

        process_junction(
            key,
            Some(&barcode),
            None, // no UMI
            &mut junction_counts,
            &mut junction_totals,
            &mut processed_reads,
            &mut processed_umis,
            12345,
            Mode::Single,
        );

        assert_eq!(
            *junction_counts.get(&key).unwrap().get("ACGT-1").unwrap(),
            1
        );
    }

    #[test]
    fn test_process_junction_umi_dedup() {
        // Same barcode + same UMI + same junction → counted once
        let mut junction_counts = HashMap::new();
        let mut junction_totals = HashMap::new();
        let mut processed_reads = HashMap::new();
        let mut processed_umis = HashMap::new();

        let key = JunctionKey { tid: 0, start: 100, end: 200, strand: Strand::Plus };
        let barcode = "AAAA-1".to_string();
        let umi = "ACGT".to_string();

        // Two different reads with same barcode+UMI
        for hash in [11111u64, 22222u64] {
            process_junction(
                key,
                Some(&barcode),
                Some(&umi),
                &mut junction_counts,
                &mut junction_totals,
                &mut processed_reads,
                &mut processed_umis,
                hash,
                Mode::Single,
            );
        }

        assert_eq!(
            *junction_counts.get(&key).unwrap().get("AAAA-1").unwrap(),
            1,
            "Same barcode+UMI should be counted only once"
        );
    }

    #[test]
    fn test_process_junction_different_umis() {
        // Same barcode, different UMIs → each counted
        let mut junction_counts = HashMap::new();
        let mut junction_totals = HashMap::new();
        let mut processed_reads = HashMap::new();
        let mut processed_umis = HashMap::new();

        let key = JunctionKey { tid: 0, start: 100, end: 200, strand: Strand::Plus };
        let barcode = "AAAA-1".to_string();
        let umi1 = "ACGT".to_string();
        let umi2 = "TGCA".to_string();

        process_junction(
            key, Some(&barcode), Some(&umi1),
            &mut junction_counts, &mut junction_totals,
            &mut processed_reads, &mut processed_umis,
            11111, Mode::Single,
        );
        process_junction(
            key, Some(&barcode), Some(&umi2),
            &mut junction_counts, &mut junction_totals,
            &mut processed_reads, &mut processed_umis,
            22222, Mode::Single,
        );

        assert_eq!(
            *junction_counts.get(&key).unwrap().get("AAAA-1").unwrap(),
            2,
            "Different UMIs should be counted separately"
        );
    }

    #[test]
    fn test_process_junction_no_umi_fallback() {
        // No UMI → falls back to read-name dedup only; different read hashes both count
        let mut junction_counts = HashMap::new();
        let mut junction_totals = HashMap::new();
        let mut processed_reads = HashMap::new();
        let mut processed_umis = HashMap::new();

        let key = JunctionKey { tid: 0, start: 100, end: 200, strand: Strand::Plus };
        let barcode = "AAAA-1".to_string();

        for hash in [11111u64, 22222u64] {
            process_junction(
                key, Some(&barcode), None,
                &mut junction_counts, &mut junction_totals,
                &mut processed_reads, &mut processed_umis,
                hash, Mode::Single,
            );
        }

        assert_eq!(
            *junction_counts.get(&key).unwrap().get("AAAA-1").unwrap(),
            2,
            "Without UMI, each unique read should be counted"
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

    #[test]
    fn test_parse_junction_key_no_colon() {
        // Key with no colon should return the whole key as coords and "." as strand
        let (coords, strand) = parse_junction_key("no_colon_key");
        assert_eq!(coords, "no_colon_key");
        assert_eq!(strand, ".");
    }
}
