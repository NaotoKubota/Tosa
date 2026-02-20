//! Exon-intron boundary read counting.
//!
//! Boundary coordinates are represented as 1-base intervals at each end of an intron.
//! For example, intron `chr2:6545675-6547042` produces:
//! - 5' boundary: `chr2:6545675-6545676`
//! - 3' boundary: `chr2:6547041-6547042`

use std::collections::{HashMap, HashSet, BTreeMap};
use crate::types::{BoundaryType, Strand};

/// A single boundary entry derived from GTF annotation.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct BoundaryEntry {
    /// 1-base boundary coordinate string, e.g. "chr2:6545675-6545676"
    pub boundary_id: String,
    /// Start position of the 1-base interval.
    pub start: i64,
    /// End position of the 1-base interval (start + 1).
    pub end: i64,
    /// Type of boundary (5' or 3').
    pub boundary_type: BoundaryType,
    /// Strand from GTF annotation.
    pub strand: Strand,
}

/// Index of all exon-intron boundaries, organized by chromosome.
///
/// For each chromosome, boundaries are stored in a BTreeMap keyed by start position
/// for efficient range lookups.
#[derive(Debug, Clone, Default)]
pub struct BoundaryIndex {
    /// chrom -> (start_pos -> Vec<BoundaryEntry>)
    pub boundaries: HashMap<String, BTreeMap<i64, Vec<BoundaryEntry>>>,
}

impl BoundaryIndex {
    /// Create an empty BoundaryIndex.
    pub fn new() -> Self {
        BoundaryIndex {
            boundaries: HashMap::new(),
        }
    }

    /// Add a boundary entry to the index.
    pub fn add(&mut self, chrom: &str, entry: BoundaryEntry) {
        self.boundaries
            .entry(chrom.to_string())
            .or_default()
            .entry(entry.start)
            .or_default()
            .push(entry);
    }

    /// Find all boundaries on a given chromosome whose 1-base interval is completely
    /// contained within [seg_start, seg_end).
    pub fn find_overlapping(&self, chrom: &str, seg_start: i64, seg_end: i64) -> Vec<&BoundaryEntry> {
        let mut results = Vec::new();
        if let Some(chrom_boundaries) = self.boundaries.get(chrom) {
            // Use BTreeMap range to efficiently find candidates
            for (_pos, entries) in chrom_boundaries.range(seg_start..seg_end) {
                for entry in entries {
                    // Boundary [entry.start, entry.end) must be fully within [seg_start, seg_end)
                    if entry.start >= seg_start && entry.end <= seg_end {
                        results.push(entry);
                    }
                }
            }
        }
        results
    }
}

/// Derive 5' and 3' boundary entries from an intron coordinate.
///
/// Given intron `chr:intron_start-intron_end`:
/// - 5' boundary: `chr:intron_start-(intron_start+1)`
/// - 3' boundary: `chr:(intron_end-1)-intron_end`
pub fn intron_to_boundaries(chrom: &str, intron_start: i64, intron_end: i64, strand: Strand) -> (BoundaryEntry, BoundaryEntry) {
    let five_prime = BoundaryEntry {
        boundary_id: format!("{}:{}-{}", chrom, intron_start, intron_start + 1),
        start: intron_start,
        end: intron_start + 1,
        boundary_type: BoundaryType::FivePrime,
        strand,
    };
    let three_prime = BoundaryEntry {
        boundary_id: format!("{}:{}-{}", chrom, intron_end - 1, intron_end),
        start: intron_end - 1,
        end: intron_end,
        boundary_type: BoundaryType::ThreePrime,
        strand,
    };
    (five_prime, three_prime)
}

/// Count boundary overlaps for a read's aligned segments.
///
/// For each aligned segment, finds boundaries in the index that are fully contained,
/// and increments counts accordingly.
#[allow(clippy::too_many_arguments)]
pub fn count_boundaries(
    chrom: &str,
    aligned_segments: &[(i64, i64)],
    boundary_index: &BoundaryIndex,
    cell_barcode: Option<&String>,
    strand: Strand,
    boundary_counts: &mut HashMap<String, HashMap<String, u32>>,
    boundary_totals: &mut HashMap<String, u32>,
    boundary_types: &mut HashMap<String, BoundaryType>,
    boundary_strands: &mut HashMap<String, Strand>,
    processed_boundary_reads: &mut HashMap<String, HashSet<String>>,
    read_name: &str,
    mode: &str,
) {
    for (seg_start, seg_end) in aligned_segments {
        let overlapping = boundary_index.find_overlapping(chrom, *seg_start, *seg_end);
        for entry in overlapping {
            let key = &entry.boundary_id;

            // Check dedup: same read should not count same boundary twice
            if let Some(reads) = processed_boundary_reads.get_mut(key) {
                if reads.contains(read_name) {
                    continue;
                }
                reads.insert(read_name.to_string());
            } else {
                let mut reads_set = HashSet::new();
                reads_set.insert(read_name.to_string());
                processed_boundary_reads.insert(key.clone(), reads_set);
            }

            // Record type and strand
            boundary_types.entry(key.clone()).or_insert(entry.boundary_type);
            boundary_strands.entry(key.clone()).or_insert(strand);

            if mode == "single" {
                if let Some(cb_str) = cell_barcode {
                    let boundary_entry = boundary_counts
                        .entry(key.clone())
                        .or_default();
                    *boundary_entry.entry(cb_str.clone()).or_insert(0) += 1;
                }
            } else {
                *boundary_totals.entry(key.clone()).or_insert(0) += 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_intron_to_boundaries() {
        let (five_p, three_p) = intron_to_boundaries("chr2", 6545675, 6547042, Strand::Plus);

        assert_eq!(five_p.boundary_id, "chr2:6545675-6545676");
        assert_eq!(five_p.start, 6545675);
        assert_eq!(five_p.end, 6545676);
        assert_eq!(five_p.boundary_type, BoundaryType::FivePrime);

        assert_eq!(three_p.boundary_id, "chr2:6547041-6547042");
        assert_eq!(three_p.start, 6547041);
        assert_eq!(three_p.end, 6547042);
        assert_eq!(three_p.boundary_type, BoundaryType::ThreePrime);
    }

    #[test]
    fn test_intron_to_boundaries_example2() {
        let (five_p, three_p) = intron_to_boundaries("chr1", 1000, 2000, Strand::Minus);

        assert_eq!(five_p.boundary_id, "chr1:1000-1001");
        assert_eq!(three_p.boundary_id, "chr1:1999-2000");
        assert_eq!(five_p.strand, Strand::Minus);
    }

    #[test]
    fn test_boundary_index_find_overlapping() {
        let mut index = BoundaryIndex::new();

        let (five_p, three_p) = intron_to_boundaries("chr2", 6545675, 6547042, Strand::Plus);
        index.add("chr2", five_p);
        index.add("chr2", three_p);

        // Aligned segment that spans the 5' boundary
        let overlapping = index.find_overlapping("chr2", 6545600, 6545700);
        assert_eq!(overlapping.len(), 1);
        assert_eq!(overlapping[0].boundary_id, "chr2:6545675-6545676");
        assert_eq!(overlapping[0].boundary_type, BoundaryType::FivePrime);

        // Aligned segment that does NOT span the 5' boundary
        let overlapping = index.find_overlapping("chr2", 6545677, 6545800);
        assert_eq!(overlapping.len(), 0);

        // Aligned segment that spans the 3' boundary
        let overlapping = index.find_overlapping("chr2", 6547000, 6547100);
        assert_eq!(overlapping.len(), 1);
        assert_eq!(overlapping[0].boundary_id, "chr2:6547041-6547042");

        // Different chromosome
        let overlapping = index.find_overlapping("chr1", 6545600, 6545700);
        assert_eq!(overlapping.len(), 0);
    }

    #[test]
    fn test_count_boundaries_bulk() {
        let mut index = BoundaryIndex::new();
        let (five_p, three_p) = intron_to_boundaries("chr1", 1000, 2000, Strand::Plus);
        index.add("chr1", five_p);
        index.add("chr1", three_p);

        let mut boundary_counts = HashMap::new();
        let mut boundary_totals = HashMap::new();
        let mut boundary_types = HashMap::new();
        let mut boundary_strands = HashMap::new();
        let mut processed = HashMap::new();

        // Segment that spans the 5' boundary at 1000-1001
        let segments = vec![(900, 1100)];
        count_boundaries(
            "chr1",
            &segments,
            &index,
            None,
            Strand::Plus,
            &mut boundary_counts,
            &mut boundary_totals,
            &mut boundary_types,
            &mut boundary_strands,
            &mut processed,
            "read1",
            "bulk",
        );

        assert_eq!(boundary_totals.get("chr1:1000-1001"), Some(&1));
        assert_eq!(boundary_totals.get("chr1:1999-2000"), None); // Not overlapping
    }

    #[test]
    fn test_count_boundaries_dedup() {
        let mut index = BoundaryIndex::new();
        let (five_p, _three_p) = intron_to_boundaries("chr1", 1000, 2000, Strand::Plus);
        index.add("chr1", five_p);

        let mut boundary_counts = HashMap::new();
        let mut boundary_totals = HashMap::new();
        let mut boundary_types = HashMap::new();
        let mut boundary_strands = HashMap::new();
        let mut processed = HashMap::new();

        let segments = vec![(900, 1100)];

        // Same read twice
        for _ in 0..2 {
            count_boundaries(
                "chr1",
                &segments,
                &index,
                None,
                Strand::Plus,
                &mut boundary_counts,
                &mut boundary_totals,
                &mut boundary_types,
                &mut boundary_strands,
                &mut processed,
                "read1",
                "bulk",
            );
        }

        assert_eq!(boundary_totals.get("chr1:1000-1001"), Some(&1));
    }
}
