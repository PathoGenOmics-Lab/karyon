//! Does "+N more" survive show_counts(false), and where does it sit?
use karyon::{Figure, MsaSequence, Region, SnpTrack};

fn main() {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut aln = vec![MsaSequence::new("reference", "ACGTACGTACGTACGTACGT")];
    for i in 0..60usize {
        let seq: String = (0..20)
            .map(|c| if (c + i) % 3 == 0 { bases[(i + c) % 4] as char } else { b'A' as char })
            .collect();
        aln.push(MsaSequence::new(format!("sample{i}"), seq));
    }
    let region = Region::new("sites", 0, 20).unwrap();
    let off = Figure::new(region)
        .push(SnpTrack::from_alignment(0, &aln).show_counts(false))
        .to_svg();
    for line in off.lines() {
        if line.contains("more") {
            println!("{}", line.trim());
        }
    }
    println!("still has 'more': {}", off.contains("more"));
}
