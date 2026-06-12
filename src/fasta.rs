/// Simple FASTA parser and translation helpers.

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InputSequenceType {
    Protein,
    Dna,
}

#[derive(Clone, Debug)]
pub struct TranslatedFrame {
    pub label: String,
    pub sequence: String,
    pub nt_offset: usize,
    pub is_reverse: bool,
    /// Positions (indices into `sequence`) where a stop codon was read through
    /// as a placeholder `X` rather than terminating the frame.
    pub stop_positions: Vec<usize>,
}

/// Placeholder residue emitted for an internal stop codon. The HMM scores any
/// non-standard symbol neutrally (averaged emissions), so a read-through stop
/// no longer fragments the ORF; the stop is surfaced via `stop_positions`.
const STOP_PLACEHOLDER: char = 'X';

pub fn parse_fasta(input: &str) -> Vec<(String, String)> {
    let mut sequences = Vec::new();
    let mut current_name = String::new();
    let mut current_seq = String::new();

    for line in input.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if !current_name.is_empty() && !current_seq.is_empty() {
                sequences.push((current_name.clone(), current_seq.clone()));
                current_seq.clear();
            }
            current_name = header.trim().to_string();
        } else {
            append_sequence_line(&mut current_seq, line);
        }
    }

    if !current_name.is_empty() && !current_seq.is_empty() {
        sequences.push((current_name, current_seq));
    }

    sequences
}

fn append_sequence_line(current_seq: &mut String, line: &str) {
    current_seq.extend(
        line.chars()
            .filter(|ch| !matches!(ch, '-' | '.'))
            .map(|ch| ch.to_ascii_uppercase()),
    );
}

pub fn validate_protein_sequence(seq: &str) -> bool {
    // `X` (unknown) and `*` (stop read-through) are scored neutrally by the HMM,
    // so they are permitted in both pasted protein input and translated frames.
    const AMINO_ACIDS: &[u8] = b"ACDEFGHIKLMNPQRSTVWYX*";
    seq.len() < 10000
        && seq
            .bytes()
            .all(|b| AMINO_ACIDS.contains(&b.to_ascii_uppercase()))
}

pub fn validate_dna_sequence(seq: &str) -> bool {
    const DNA: &[u8] = b"ACGTUNRYKMSWBDHV";
    seq.len() < 100000
        && seq
            .bytes()
            .all(|b| DNA.contains(&b.to_ascii_uppercase()))
}

pub fn translate_six_frames(dna: &str) -> Vec<TranslatedFrame> {
    let dna = dna.as_bytes();
    let revcomp = reverse_complement(dna);
    let mut translated = Vec::new();

    for offset in 0..3 {
        translated.extend(translate_frame(dna, offset, false));
    }
    for offset in 0..3 {
        translated.extend(translate_frame(&revcomp, offset, true));
    }

    translated
}

fn translate_frame(sequence: &[u8], offset: usize, is_reverse: bool) -> Vec<TranslatedFrame> {
    let frame_no = offset + 1;
    let label = if is_reverse {
        format!("-{frame_no}")
    } else {
        format!("+{frame_no}")
    };

    let mut aa = String::new();
    let mut stop_positions = Vec::new();

    for codon_start in (offset..=sequence.len().saturating_sub(3)).step_by(3) {
        let codon = &sequence[codon_start..codon_start + 3];
        match translate_codon(codon) {
            Codon::Residue(residue) => aa.push(residue),
            Codon::Stop => {
                // Read through the stop as a neutral placeholder, keeping one
                // continuous peptide, and remember where it occurred so the
                // caller can flag it (e.g. a base-calling error vs. read-through).
                stop_positions.push(aa.chars().count());
                aa.push(STOP_PLACEHOLDER);
            }
            // Ambiguous codon (e.g. contains N): unknown residue, also scored
            // neutrally, but not a stop — so it is not flagged for review.
            Codon::Unknown => aa.push(STOP_PLACEHOLDER),
        }
    }

    if aa.is_empty() {
        return Vec::new();
    }

    vec![TranslatedFrame {
        label,
        sequence: aa,
        nt_offset: offset,
        is_reverse,
        stop_positions,
    }]
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' | b'U' => b'A',
            b'R' => b'Y',
            b'Y' => b'R',
            b'K' => b'M',
            b'M' => b'K',
            b'S' => b'S',
            b'W' => b'W',
            b'B' => b'V',
            b'D' => b'H',
            b'H' => b'D',
            b'V' => b'B',
            _ => b'N',
        })
        .collect()
}

enum Codon {
    Residue(char),
    Stop,
    Unknown,
}

fn translate_codon(codon: &[u8]) -> Codon {
    let mut normalized = [b'N'; 3];
    for (i, base) in codon.iter().enumerate().take(3) {
        normalized[i] = match base.to_ascii_uppercase() {
            b'U' => b'T',
            other => other,
        };
    }

    match &normalized {
        b"TTT" | b"TTC" => Codon::Residue('F'),
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => Codon::Residue('L'),
        b"ATT" | b"ATC" | b"ATA" => Codon::Residue('I'),
        b"ATG" => Codon::Residue('M'),
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => Codon::Residue('V'),
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => Codon::Residue('S'),
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => Codon::Residue('P'),
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => Codon::Residue('T'),
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => Codon::Residue('A'),
        b"TAT" | b"TAC" => Codon::Residue('Y'),
        b"TAA" | b"TAG" | b"TGA" => Codon::Stop,
        b"CAT" | b"CAC" => Codon::Residue('H'),
        b"CAA" | b"CAG" => Codon::Residue('Q'),
        b"AAT" | b"AAC" => Codon::Residue('N'),
        b"AAA" | b"AAG" => Codon::Residue('K'),
        b"GAT" | b"GAC" => Codon::Residue('D'),
        b"GAA" | b"GAG" => Codon::Residue('E'),
        b"TGT" | b"TGC" => Codon::Residue('C'),
        b"TGG" => Codon::Residue('W'),
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => Codon::Residue('R'),
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => Codon::Residue('G'),
        _ => Codon::Unknown,
    }
}

#[cfg(test)]
mod tests {
    use super::{parse_fasta, translate_frame};

    #[test]
    fn internal_stop_is_read_through_as_x_in_one_frame() {
        // ATG GAA TGA AAA = M E (stop) K  -> "MEXK", one continuous frame,
        // stop recorded at peptide index 2.
        let dna = b"ATGGAATGAAAA";
        let frames = translate_frame(dna, 0, false);
        assert_eq!(frames.len(), 1, "stop must not split the frame");
        assert_eq!(frames[0].sequence, "MEXK");
        assert_eq!(frames[0].stop_positions, vec![2]);
    }

    #[test]
    fn ambiguous_codon_is_x_but_not_flagged_as_stop() {
        // ATG GAN AAA -> M X K ; the N-containing codon is unknown, not a stop.
        let dna = b"ATGGANAAA";
        let frames = translate_frame(dna, 0, false);
        assert_eq!(frames[0].sequence, "MXK");
        assert!(frames[0].stop_positions.is_empty());
    }

    #[test]
    fn parse_fasta_strips_alignment_gap_markers() {
        let fasta = ">seq1\nac-t.g\n>seq2\nA.C-G\n";
        let parsed = parse_fasta(fasta);

        assert_eq!(
            parsed,
            vec![
                ("seq1".to_string(), "ACTG".to_string()),
                ("seq2".to_string(), "ACG".to_string()),
            ]
        );
    }
}
