// Modules for handling exon-intron boundaries

pub fn count_exon_intron_boundaries(
    cell_barcode: Option<&String>,
    introns: &HashSet<(String, i64, i64)>,
    chrom: &str,
    start: i64,
    end: i64,
    left_counts: &mut HashMap<String, HashMap<String, u32>>,
    right_counts: &mut HashMap<String, HashMap<String, u32>>,
    left_totals: &mut HashMap<String, u32>,
    right_totals: &mut HashMap<String, u32>,
    processed_boundary_reads: &mut HashMap<String, HashSet<String>>,
    read_name: &str,
    mode: &str,
) {
    // Check if the read was already processed for this boundary
    if let Some(reads) = processed_boundary_reads.get_mut(chrom) {
        if reads.contains(read_name) {
            return; // Skip counting
        }
        reads.insert(read_name.to_string());
    } else {
        let mut reads_set = HashSet::new();
        reads_set.insert(read_name.to_string());
        processed_boundary_reads.insert(chrom.to_string(), reads_set);
    }

    // Count the read for the boundary
    if mode == "single" {
        if let Some(cb_str) = cell_barcode {
            for (intron_chrom, intron_start, intron_end) in introns {
                let key = format!("{}:{}-{}", intron_chrom, intron_start, intron_end);
                if chrom == intron_chrom {
                    if start <= *intron_start && end >= *intron_start {
                        let boundary_entry = left_counts
                            .entry(key.to_string())
                            .or_insert_with(HashMap::new);
                        *boundary_entry.entry(cb_str.clone()).or_insert(0) += 1;
                    }
                    if start <= *intron_end && end >= *intron_end {
                        let boundary_entry = right_counts
                            .entry(key.to_string())
                            .or_insert_with(HashMap::new);
                        *boundary_entry.entry(cb_str.clone()).or_insert(0) += 1;
                    }
                }
            }
        }
    } else {
        for (intron_chrom, intron_start, intron_end) in introns {
            let key = format!("{}:{}-{}", intron_chrom, intron_start, intron_end);
            if chrom == intron_chrom {
                if start <= *intron_start && end >= *intron_start {
                    *left_totals.entry(key.to_string()).or_insert(0) += 1;
                }
                if start <= *intron_end && end >= *intron_end {
                    *right_totals.entry(key.to_string()).or_insert(0) += 1;
                }
            }
        }
    }
}