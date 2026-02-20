//! Generates a small synthetic BAM file for trying out Tosa.
//!
//! Usage:
//!   cargo run --example generate_example_bam [output_dir]
//!
//! Creates `example.bam` and `example.bam.bai` in the specified directory
//! (defaults to `examples/`).
//!
//! The BAM contains reads mapped to chr1 with:
//! - 13 unique junction reads spanning intron 1 (chr1:1201-1499)
//!   - 10 with XS:A:+ and 3 with XS:A:-
//!   - 1 duplicate pair (same QNAME) for dedup testing
//! - 5 junction reads spanning intron 2 (chr1:1701-1999), all XS:A:+
//! - 6 non-spliced reads overlapping exon-intron boundaries
//! - 1 multi-mapped read (NH:i:2, filtered by default)
//! - All reads carry a CB (cell barcode) tag for single-cell mode testing:
//!   AAAA-1 (junction 1 + strand), BBBB-1 (junction 1 - strand),
//!   CCCC-1 (junction 2 + boundary reads)
//!
//! Matching GTF annotation is in `examples/annotation.gtf`.

use rust_htslib::bam::{self, record::Aux, record::Cigar, record::CigarString, Record, Writer};
use rust_htslib::bam::header::{Header, HeaderRecord};

fn main() {
    let out_dir = std::env::args().nth(1).unwrap_or_else(|| "examples".to_string());
    let bam_path = format!("{}/example.bam", out_dir);

    // --- BAM header: one reference sequence chr1 (10 kb) ---
    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", "10000");
    header.push_record(&sq);

    let mut writer = Writer::from_path(&bam_path, &header, bam::Format::Bam)
        .expect("Failed to create BAM writer");

    // --- Gene model (GTF 1-based, + strand) ---
    //   Exon 1: 1000–1200   Exon 2: 1500–1700   Exon 3: 2000–2300
    //   Intron 1: 1201–1499 (0-based 1200..1499, len 299)
    //   Intron 2: 1701–1999 (0-based 1700..1999, len 299)
    //
    // Junction keys (1-based): chr1:1201-1499, chr1:1701-1999
    // Boundaries (0-based pairs):
    //   Intron 1 → 5': chr1:1200-1201  3': chr1:1498-1499
    //   Intron 2 → 5': chr1:1700-1701  3': chr1:1998-1999

    struct ReadSpec {
        qname: &'static str,
        pos: i64,
        cigar: Vec<Cigar>,
        xs: Option<u8>,
        nh: u8,
        cb: &'static str,
    }

    // All reads are 100 bp, sorted by position for coordinate-sorted BAM.
    // + and - strand reads are interleaved by position to maintain sort order.
    // Cell barcodes: AAAA-1 = junction1 + strand, BBBB-1 = junction1 - strand, CCCC-1 = junction2 + boundaries
    let reads: Vec<ReadSpec> = vec![
        // ── Junction 1 reads (intron at 0-based 1200..1499) ─────────
        ReadSpec { qname: "read01", pos: 1140, cigar: vec![Cigar::Match(60), Cigar::RefSkip(299), Cigar::Match(40)], xs: Some(b'+'), nh: 1, cb: "AAAA-1" },
        ReadSpec { qname: "read02", pos: 1145, cigar: vec![Cigar::Match(55), Cigar::RefSkip(299), Cigar::Match(45)], xs: Some(b'+'), nh: 1, cb: "AAAA-1" },
        ReadSpec { qname: "read03", pos: 1148, cigar: vec![Cigar::Match(52), Cigar::RefSkip(299), Cigar::Match(48)], xs: Some(b'+'), nh: 1, cb: "AAAA-1" },
        ReadSpec { qname: "read10", pos: 1148, cigar: vec![Cigar::Match(52), Cigar::RefSkip(299), Cigar::Match(48)], xs: Some(b'-'), nh: 1, cb: "BBBB-1" },
        ReadSpec { qname: "read04", pos: 1150, cigar: vec![Cigar::Match(50), Cigar::RefSkip(299), Cigar::Match(50)], xs: Some(b'+'), nh: 1, cb: "AAAA-1" },
        ReadSpec { qname: "read05", pos: 1150, cigar: vec![Cigar::Match(50), Cigar::RefSkip(299), Cigar::Match(50)], xs: Some(b'+'), nh: 1, cb: "AAAA-1" },
        // Duplicate pair: same QNAME → counted only once
        ReadSpec { qname: "read_dup", pos: 1150, cigar: vec![Cigar::Match(50), Cigar::RefSkip(299), Cigar::Match(50)], xs: Some(b'+'), nh: 1, cb: "AAAA-1" },
        ReadSpec { qname: "read_dup", pos: 1150, cigar: vec![Cigar::Match(50), Cigar::RefSkip(299), Cigar::Match(50)], xs: Some(b'+'), nh: 1, cb: "AAAA-1" },
        ReadSpec { qname: "read06", pos: 1150, cigar: vec![Cigar::Match(50), Cigar::RefSkip(299), Cigar::Match(50)], xs: Some(b'+'), nh: 1, cb: "AAAA-1" },
        ReadSpec { qname: "read11", pos: 1150, cigar: vec![Cigar::Match(50), Cigar::RefSkip(299), Cigar::Match(50)], xs: Some(b'-'), nh: 1, cb: "BBBB-1" },
        ReadSpec { qname: "read07", pos: 1152, cigar: vec![Cigar::Match(48), Cigar::RefSkip(299), Cigar::Match(52)], xs: Some(b'+'), nh: 1, cb: "AAAA-1" },
        ReadSpec { qname: "read12", pos: 1152, cigar: vec![Cigar::Match(48), Cigar::RefSkip(299), Cigar::Match(52)], xs: Some(b'-'), nh: 1, cb: "BBBB-1" },
        ReadSpec { qname: "read08", pos: 1155, cigar: vec![Cigar::Match(45), Cigar::RefSkip(299), Cigar::Match(55)], xs: Some(b'+'), nh: 1, cb: "AAAA-1" },
        ReadSpec { qname: "read09", pos: 1160, cigar: vec![Cigar::Match(40), Cigar::RefSkip(299), Cigar::Match(60)], xs: Some(b'+'), nh: 1, cb: "AAAA-1" },

        // ── Boundary reads (non-spliced, overlap exon-intron boundaries) ─
        ReadSpec { qname: "read13", pos: 1190, cigar: vec![Cigar::Match(100)], xs: Some(b'+'), nh: 1, cb: "CCCC-1" },  // 5' boundary intron 1
        ReadSpec { qname: "read14", pos: 1195, cigar: vec![Cigar::Match(100)], xs: Some(b'+'), nh: 1, cb: "CCCC-1" },  // 5' boundary intron 1
        ReadSpec { qname: "read15", pos: 1450, cigar: vec![Cigar::Match(100)], xs: Some(b'+'), nh: 1, cb: "CCCC-1" },  // 3' boundary intron 1
        ReadSpec { qname: "read16", pos: 1455, cigar: vec![Cigar::Match(100)], xs: Some(b'+'), nh: 1, cb: "CCCC-1" },  // 3' boundary intron 1

        // ── Junction 2 reads (intron at 0-based 1700..1999) ─────────
        ReadSpec { qname: "read17", pos: 1645, cigar: vec![Cigar::Match(55), Cigar::RefSkip(299), Cigar::Match(45)], xs: Some(b'+'), nh: 1, cb: "CCCC-1" },
        ReadSpec { qname: "read18", pos: 1648, cigar: vec![Cigar::Match(52), Cigar::RefSkip(299), Cigar::Match(48)], xs: Some(b'+'), nh: 1, cb: "CCCC-1" },
        ReadSpec { qname: "read19", pos: 1650, cigar: vec![Cigar::Match(50), Cigar::RefSkip(299), Cigar::Match(50)], xs: Some(b'+'), nh: 1, cb: "CCCC-1" },
        ReadSpec { qname: "read20", pos: 1652, cigar: vec![Cigar::Match(48), Cigar::RefSkip(299), Cigar::Match(52)], xs: Some(b'+'), nh: 1, cb: "CCCC-1" },
        ReadSpec { qname: "read21", pos: 1655, cigar: vec![Cigar::Match(45), Cigar::RefSkip(299), Cigar::Match(55)], xs: Some(b'+'), nh: 1, cb: "CCCC-1" },

        // ── More boundary reads ─────────────────────────────────────
        ReadSpec { qname: "read22", pos: 1690, cigar: vec![Cigar::Match(100)], xs: Some(b'+'), nh: 1, cb: "CCCC-1" },  // 5' boundary intron 2
        ReadSpec { qname: "read23", pos: 1950, cigar: vec![Cigar::Match(100)], xs: Some(b'+'), nh: 1, cb: "CCCC-1" },  // 3' boundary intron 2

        // ── Multi-mapped read (filtered by default NH ≤ 1) ──────────
        ReadSpec { qname: "read_multi", pos: 5000, cigar: vec![Cigar::Match(100)], xs: Some(b'+'), nh: 2, cb: "AAAA-1" },
    ];

    for rd in &reads {
        let mut rec = Record::new();
        let cigar_str = CigarString(rd.cigar.clone());

        // Query length = sum of query-consuming CIGAR ops (M, I, S, =, X)
        let qlen: u32 = rd.cigar.iter().map(|c| match c {
            Cigar::Match(l) | Cigar::Ins(l) | Cigar::SoftClip(l)
            | Cigar::Equal(l) | Cigar::Diff(l) => *l,
            _ => 0,
        }).sum();

        let seq = vec![b'A'; qlen as usize];
        let qual = vec![40u8; qlen as usize];

        rec.set(rd.qname.as_bytes(), Some(&cigar_str), &seq, &qual);
        rec.set_flags(0);  // mapped, forward strand, unpaired
        rec.set_tid(0);   // chr1
        rec.set_pos(rd.pos);
        rec.set_mapq(255);
        rec.push_aux(b"NH", Aux::U8(rd.nh)).unwrap();
        if let Some(xs) = rd.xs {
            rec.push_aux(b"XS", Aux::Char(xs)).unwrap();
        }
        rec.push_aux(b"CB", Aux::String(rd.cb)).unwrap();

        writer.write(&rec).unwrap();
    }

    drop(writer);

    // Build BAI index
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

    println!("Created {bam_path} and {bam_path}.bai");
    println!("  26 records on chr1 (25 unique names + 1 duplicate pair)");
    println!("  Junction 1: chr1:1201-1499 (13 unique reads: 10 +strand, 3 -strand)");
    println!("  Junction 2: chr1:1701-1999 (5 reads, all +strand)");
    println!("  Boundary reads: 6 (2 per intron-end × 2 introns + 2 extra)");
    println!("  Multi-mapped (NH=2, filtered): 1");
}
