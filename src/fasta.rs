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
}

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
    const AMINO_ACIDS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
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

    let mut frames = Vec::new();
    let mut aa = String::new();
    let mut current_nt_offset = offset;

    for codon_start in (offset..=sequence.len().saturating_sub(3)).step_by(3) {
        let codon = &sequence[codon_start..codon_start + 3];
        match translate_codon(codon) {
            Some(residue) => aa.push(residue),
            None => {
                if !aa.is_empty() {
                    frames.push(TranslatedFrame {
                        label: label.clone(),
                        sequence: aa.clone(),
                        nt_offset: current_nt_offset,
                        is_reverse,
                    });
                    aa.clear();
                }
                current_nt_offset = codon_start + 3;
            }
        }
    }

    if !aa.is_empty() {
        frames.push(TranslatedFrame {
            label,
            sequence: aa,
            nt_offset: current_nt_offset,
            is_reverse,
        });
    }

    frames
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

fn translate_codon(codon: &[u8]) -> Option<char> {
    let mut normalized = [b'N'; 3];
    for (i, base) in codon.iter().enumerate().take(3) {
        normalized[i] = match base.to_ascii_uppercase() {
            b'U' => b'T',
            other => other,
        };
    }

    match &normalized {
        b"TTT" | b"TTC" => Some('F'),
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => Some('L'),
        b"ATT" | b"ATC" | b"ATA" => Some('I'),
        b"ATG" => Some('M'),
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => Some('V'),
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => Some('S'),
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => Some('P'),
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => Some('T'),
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => Some('A'),
        b"TAT" | b"TAC" => Some('Y'),
        b"TAA" | b"TAG" | b"TGA" => None,
        b"CAT" | b"CAC" => Some('H'),
        b"CAA" | b"CAG" => Some('Q'),
        b"AAT" | b"AAC" => Some('N'),
        b"AAA" | b"AAG" => Some('K'),
        b"GAT" | b"GAC" => Some('D'),
        b"GAA" | b"GAG" => Some('E'),
        b"TGT" | b"TGC" => Some('C'),
        b"TGG" => Some('W'),
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => Some('R'),
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => Some('G'),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::parse_fasta;

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
