pub mod dna;
pub mod protein;

use bit_set::BitSet;
use std::borrow::Borrow;
use vector_map::VecMap;

pub type SymbolRanks = VecMap<usize, u8>;

#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Alphabet {
    pub symbols: BitSet,
}

impl Alphabet {
    pub fn new<C, T>(symbols: T) -> Self
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        let mut s = BitSet::new();
        s.extend(symbols.into_iter().map(|c| *c.borrow() as usize));

        Alphabet { symbols: s }
    }

    pub fn insert(&mut self, a: u8) {
        self.symbols.insert(a as usize);
    }

    pub fn is_word<C, T>(&self, text: T) -> bool
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        text.into_iter()
            .all(|c| self.symbols.contains(*c.borrow() as usize))
    }

    pub fn max_symbol(&self) -> Option<u8> {
        self.symbols.iter().max().map(|a| a as u8)
    }

    pub fn len(&self) -> usize {
        self.symbols.len()
    }

    pub fn is_empty(&self) -> bool {
        self.symbols.is_empty()
    }

    pub fn intersection(&self, other: &Alphabet) -> Self {
        return Alphabet {
            symbols: self.symbols.intersection(&other.symbols).collect(),
        };
    }

    pub fn difference(&self, others: &Alphabet) -> Self {
        return Alphabet {
            symbols: self.symbols.difference(&others.symbols).collect(),
        };
    }

    pub fn union(&self, others: &Alphabet) -> Self {
        return Alphabet {
            symbols: self.symbols.union(&others.symbols).collect(),
        };
    }
}

#[derive(Default, Clone, Debug)]
pub struct RankTransform {
    pub ranks: SymbolRanks,
}

impl RankTransform {
    pub fn new(alphabet: &Alphabet) -> Self {
        let mut ranks = VecMap::new();
        for (r, c) in alphabet.symbols.iter().enumerate() {
            ranks.insert(c as usize, r as u8);
        }
        RankTransform { ranks }
    }

    #[inline]
    pub fn get(&self, a: u8) -> u8 {
        *self.ranks.get(&(a as usize)).expect("Unexpected character.")
    }

    pub fn transform<C, T>(&self, text: T) -> Vec<u8>
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        text.into_iter()
            .map(|c| {
                *self
                    .ranks
                    .get(&(*c.borrow() as usize))
                    .expect("Unexpected character in text.")
            })
            .collect()
    }

    pub fn qgrams<C, T>(&self, q: u32, text: T) -> QGrams<'_, C, T::IntoIter>
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        assert!(q > 0, "Expecting q-gram length a to be larger than 0.");
        let bits = (self.ranks.len() as f32).log2().ceil() as u32;
        assert!(
            bits * q <= usize::BITS,
            "Expecting q to be smaller than usize /log2(|A|)"
        );

        let mut qgrams = QGrams {
            text: text.into_iter(),
            ranks: self,
            bits,
            mask: 1usize.checked_shl(q * bits).unwrap_or(0).wrapping_sub(1),
            qgram: 0,
        };

        for _ in 0..q - 1 {
            qgrams.next();
        }

        qgrams
    }

    pub fn rev_qgrams<C, IT, T>(&self, q: u32, text: IT) -> RevQGrams<'_, C, T>
    where
        C: Borrow<u8>,
        T: DoubleEndedIterator<Item = C>,
        IT: IntoIterator<IntoIter = T>,
    {
        assert!(q > 0, "Expecting q-gram length q to be larger than 0.");
        let bits = (self.ranks.len() as f32).log2().ceil() as u32;
        assert!(
            bits * q <= usize::BITS,
            "Expecting q to be smaller than usize / log2(|A|)"
        );

        let mut rev_qgrams = RevQGrams {
            text: text.into_iter(),
            ranks: self,
            bits,
            left_shift: (q - 1) * bits,
            qgram: 0,
        };

        for _ in 0..q - 1 {
            rev_qgrams.next();
        }

        rev_qgrams
    }

    pub fn alphabet(&self) -> Alphabet {
        let mut symbols = BitSet::with_capacity(self.ranks.len());
        symbols.extend(self.ranks.keys().copied());
        Alphabet { symbols }
    }

    pub fn get_width(&self) -> usize {
        (self.ranks.len() as f32).log2().ceil() as usize
    }
}

#[derive(Copy, Clone, Debug)]
pub struct QGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    text: T,
    ranks: & 'a RankTransform,
    bits: u32,
    mask: usize,
    qgram: usize,
}

impl<'a, C, T> QGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    #[inline]
    fn qgram_push(&mut self, a: u8) {
        self.qgram <<= self.bits;
        self.qgram |= a as usize;
        self.qgram &= self.mask;
    }
}

impl<'a, C, T> Iterator for QGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<usize> {
        match self.text.next() {
            Some(a) => {
                let b = self.ranks.get(*a.borrow());
                self.qgram_push(b);
                Some(self.qgram)
            }
            None => None,
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.text.size_hint()
    }
}

impl<'a, C, T> ExactSizeIterator for QGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: ExactSizeIterator<Item = C>,
{
}

#[derive(Copy, Clone, Debug)]
pub struct RevQGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: DoubleEndedIterator<Item = C>,
{
    text: T,
    ranks: &'a RankTransform,
    bits: u32,
    left_shift: u32,
    qgram: usize,
}

impl<'a, C, T> RevQGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: DoubleEndedIterator<Item = C>,
{
    #[inline]
    fn qgram_push_rev(&mut self, a: u8) {
        self.qgram >>= self.bits;
        self.qgram |= (a as usize) << self.left_shift;
    }
}

impl<'a, C, T> Iterator for RevQGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: DoubleEndedIterator<Item = C>,
{
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<usize> {
        match self.text.next_back() {
            Some(a) => {
                let b = self.ranks.get(*a.borrow());
                self.qgram_push_rev(b);
                Some(self.qgram)
            }
            None => None,
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.text.size_hint()
    }
}

impl<'a, C, T> ExactSizeIterator for RevQGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: DoubleEndedIterator<Item = C> + ExactSizeIterator<Item = C>,
{
}

pub fn english_ascii_lower_alphabet() -> Alphabet {
    Alphabet::new(&b"abcdefghijklmnopqrstuvwxyz"[..])
}

pub fn english_ascii_upper_alphabet() -> Alphabet {
    Alphabet::new(&b"ABCDEFGHIJKLMNOPQRSTUVWXYZ"[..])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alphabet_eq() {
        assert_eq!(Alphabet::new(b"ATCG"), Alphabet::new(b"ATCG"));
        assert_eq!(Alphabet::new(b"ATCG"), Alphabet::new(b"TAGC"));
        assert_ne!(Alphabet::new(b"ATCG"), Alphabet::new(b"ATC"));
    }

    #[test]
    fn test_qgram_shiftleft_overflow() {
        let alphabet = Alphabet::new(b"ACTG");
        let transform = RankTransform::new(&alphabet);
        let text = b"ACTG".repeat(100);
        transform.qgrams(usize::BITS / 2, text);
    }

    #[test]
    fn test_rev_qgrams() {
        let alphabet = Alphabet::new(b"ACTG");
        let transform = RankTransform::new(&alphabet);
        let text = b"ACTGA";
        let mut rev_qgrams = transform.rev_qgrams(4, text);
        assert_eq!(rev_qgrams.next(), transform.qgrams(4, b"CTGA").next());
        assert_eq!(rev_qgrams.next(), transform.qgrams(4, b"ACTG").next());
        assert_eq!(rev_qgrams.next(), None);
    }

    #[test]
    fn test_exactsize_iterator() {
        let alphabet = Alphabet::new(b"ACTG");
        let transform = RankTransform::new(&alphabet);
        let text = b"ACTGACTG";

        let mut qgrams = transform.qgrams(4, text);
        assert_eq!(qgrams.len(), 5);
        qgrams.next();
        assert_eq!(qgrams.len(), 4);

        let mut rev_qgrams = transform.rev_qgrams(4, text);
        assert_eq!(rev_qgrams.len(), 5);
        rev_qgrams.next();
        assert_eq!(rev_qgrams.len(), 4);

        let qgrams = transform.qgrams(4, b"AC");
        assert_eq!(qgrams.len(), 0);
        let rev_qgrams = transform.rev_qgrams(4, b"AC");
        assert_eq!(rev_qgrams.len(), 0);
    }
}
