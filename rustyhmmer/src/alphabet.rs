











pub const K: usize = 20;


pub const KP: usize = 29;




pub const SYM: &[u8; KP] = b"ACDEFGHIKLMNPQRSTVWY-BJZOUX*~";


pub const GAP: u8 = 20; 
pub const AMINO_B: u8 = 21;
pub const AMINO_J: u8 = 22;
pub const AMINO_Z: u8 = 23;
pub const AMINO_O: u8 = 24;
pub const AMINO_U: u8 = 25;
pub const AMINO_X: u8 = 26; 
pub const NONRESIDUE: u8 = 27; 
pub const MISSING: u8 = 28; 


pub const ILLEGAL: u8 = 255;







pub struct AminoMap {
    inmap: [u8; 256],
}

impl AminoMap {
    pub fn new() -> Self {
        let mut inmap = [ILLEGAL; 256];
        for (code, &c) in SYM.iter().enumerate() {
            inmap[c as usize] = code as u8;
            inmap[c.to_ascii_lowercase() as usize] = code as u8;
        }
        AminoMap { inmap }
    }

    
    #[inline]
    pub fn code(&self, c: u8) -> u8 {
        self.inmap[c as usize]
    }

    
    
    
    
    pub fn digitize(&self, seq: &[u8]) -> Vec<u8> {
        let mut dsq = Vec::with_capacity(seq.len() + 2);
        dsq.push(255u8); 
        for &c in seq {
            if c.is_ascii_whitespace() {
                continue;
            }
            let code = self.inmap[c as usize];
            dsq.push(if code == ILLEGAL { AMINO_X } else { code });
        }
        dsq.push(255u8); 
        dsq
    }
}

impl Default for AminoMap {
    fn default() -> Self {
        Self::new()
    }
}





pub const AMINO_FREQ: [f64; K] = [
    0.0787945, 
    0.0151600, 
    0.0535222, 
    0.0668298, 
    0.0397062, 
    0.0695071, 
    0.0229198, 
    0.0590092, 
    0.0594422, 
    0.0963728, 
    0.0237718, 
    0.0414386, 
    0.0482904, 
    0.0395639, 
    0.0540978, 
    0.0683364, 
    0.0540687, 
    0.0673417, 
    0.0114135, 
    0.0304133, 
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symbol_table_len() {
        assert_eq!(SYM.len(), KP);
        assert_eq!(&SYM[..K], b"ACDEFGHIKLMNPQRSTVWY");
    }

    #[test]
    fn canonical_codes() {
        let m = AminoMap::new();
        assert_eq!(m.code(b'A'), 0);
        assert_eq!(m.code(b'Y'), 19);
        assert_eq!(m.code(b'a'), 0); 
        assert_eq!(m.code(b'X'), AMINO_X);
        assert_eq!(m.code(b'@'), ILLEGAL);
    }

    #[test]
    fn digitize_sentinels() {
        let m = AminoMap::new();
        let dsq = m.digitize(b"ACDY");
        assert_eq!(dsq, vec![255, 0, 1, 2, 19, 255]);
    }

    #[test]
    fn background_sums_to_one() {
        let s: f64 = AMINO_FREQ.iter().sum();
        assert!((s - 1.0).abs() < 1e-6, "sum={s}");
    }
}
