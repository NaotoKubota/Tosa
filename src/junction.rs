// Modules for handling junctions
use std::collections::{HashMap, HashSet};

pub fn process_junction(
    junction_coords: &str,
    cell_barcode: Option<&String>,
    junction_counts: &mut HashMap<String, HashMap<String, u32>>,
    junction_totals: &mut HashMap<String, u32>,
    processed_reads: &mut HashMap<String, HashSet<String>>,
    read_name: &str,                                        // Read name for tracking
    mode: &str,
) {
    // Check if the read was already processed for this junction
    if let Some(reads) = processed_reads.get_mut(junction_coords) {
        if reads.contains(read_name) {
            return; // Skip counting
        }
        reads.insert(read_name.to_string());
    } else {
        let mut reads_set = HashSet::new();
        reads_set.insert(read_name.to_string());
        processed_reads.insert(junction_coords.to_string(), reads_set);
    }

    // Count the read for the junction
    if mode == "single" {
        if let Some(cb_str) = cell_barcode {
            let junction_entry = junction_counts
                .entry(junction_coords.to_string())
                .or_insert_with(HashMap::new);
            *junction_entry.entry(cb_str.clone()).or_insert(0) += 1;
        }
    } else {
        *junction_totals
            .entry(junction_coords.to_string())
            .or_insert(0) += 1;
    }
}
