use crate::alphabets::Alphabet;
use std::sync::LazyLock;

pub fn alphabet() -> Alphabet {
    Alphabet::new(b"ACGTacgt")
}

pub fn n_alphabet() -> Alphabet {
    Alphabet::new(b"ACGTNacgtn")
}

pub fn iupac_alphabet() -> Alphabet {
    Alphabet::new(b"ACGTRYSWKMBDHVNZacgtryswkmbdhvnz")
}

static COMPLEMENT: LazyLock<[u8; 256]> = LazyLock::new(|| {
    let mut comp = [0; 256];
    comp.iter_mut().enumerate().for_each(|(v, a)| {
        *a = v as u8;
    });
    b"AGCTYRWSKMDVHBN"
        .iter()
        .zip(b"TCGARYWSMKHBDVN".iter())
        .for_each(|(&a, &b)| {
            comp[a as usize] = b;
            comp[a as usize + 32] = b + 32;
        });
    comp
});

#[inline]
pub fn complement(a: u8) -> u8 {
    COMPLEMENT[a as usize]
}

pub fn reverse_complement(text: &[u8]) -> Vec<u8> {
    text.iter().rev().map(|&a| complement(a)).collect()
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_word() {
        assert!(alphabet().is_word(b"GATTACA"));
    }

    #[test]
    fn is_no_word() {
        assert!(!alphabet().is_word(b"gaUUaca"));
    }

    #[test]
    fn symbol_is_no_word() {
        assert!(!alphabet().is_word(b"#"));
    }

    #[test]
    fn number_is_no_word() {
        assert!(!alphabet().is_word(b"42"));
    }
}