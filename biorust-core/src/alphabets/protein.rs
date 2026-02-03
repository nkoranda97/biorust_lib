use crate::alphabets::Alphabet;

pub fn alphabet() -> Alphabet {
    Alphabet::new(&b"ARNDCEQGHILKMFPSTWYVarndceqghilkmfpstwyv"[..])
}

pub fn iupac_alphabet() -> Alphabet {
    Alphabet::new(b"ABCDEFGHIKLMNPQRSTVWXYZ*abcdefghiklmnpqrstvwxyz")
}
