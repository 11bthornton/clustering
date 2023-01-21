use serde::{Serialize, Deserialize};


pub mod sequence {
    // Fixed array sizes help things go a bit quicker.
    // Also we only need to deal with bytes and not strings/characters.

    use std::{
        fmt::{
            Debug, Display
        },
        ops::Deref, array::TryFromSliceError, collections::HashMap
    };

    use lazy_static::lazy_static;

    lazy_static! {
        static ref QUAL_LOOKUP: HashMap<char, f64> = {
            {
                let mut lookup : HashMap<char, f64> = HashMap::new();

                for i in 0u8..43 {

                    let c: char = (i + 33) as char;

                    lookup.insert(c, 10f64.powf(i as f64 / 10f64) as f64);
                }

                lookup
            }
        };
    }

    use bio::alignment::{
        pairwise::Aligner, 
        AlignmentOperation::*, Alignment
    };
    use itertools::Itertools;
    use serde::{Serialize, Deserialize};

    #[derive(Clone, Eq, Hash, PartialEq, Serialize, Deserialize)]
    pub struct Sequence<const N: usize>(
        #[serde(with = "serde_arrays")]
        [u8; N],

        #[serde(with = "serde_arrays")]
        [u8; N]
    );

    const CDR_LENGTH: usize = 63;
    const BARCODE_LENGTH: usize = 28;

    pub type CDR = Sequence<CDR_LENGTH>;
    pub type UMI = Sequence<BARCODE_LENGTH>;

    // impl<const N: usize> From<&[u8]> for Sequence<N> {
    //     fn from(v: &[u8]) -> Self {
    //         Sequence(v.try_into().unwrap())
    //     }
    // }

    impl<const N: usize> TryFrom<(&[u8], &[u8])> for Sequence<N> {
        type Error = TryFromSliceError;

        fn try_from(value: (&[u8], &[u8])) -> Result<Self, Self::Error> {
            Ok(
                Sequence(value.0.try_into()?, value.0.try_into()?)
            )
        }
    }

    impl<const N: usize> Display for Sequence<N> {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "{}", std::str::from_utf8(&self.0).unwrap())
        }
    }

    impl<const N: usize> Debug for Sequence<N> {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "{}", std::str::from_utf8(&self.0).unwrap())
        }
    }

    impl<const N: usize> Deref for Sequence<N> {
        type Target = [u8];

        fn deref(&self) -> &Self::Target {
            &self.0
        }
    }

    impl<const N: usize> Sequence<N> {

        pub fn ratio_ns(&self) -> f32 {

            self
                .iter()
                .filter(|char| **char == b'N')
                .collect_vec()
                .len() as f32 / (N as f32)
        }

        pub fn calculate_quality(&self) -> f64 {
            let mut total : f64 = 0.0;

            for thing in self.0 {

                let as_char = (thing + 33) as char;

                total += QUAL_LOOKUP[&as_char];
            }

            total
        }

        pub fn to_protein(&self) -> Sequence<{N / 3}> {

            let seq = self
                .array_chunks::<3>()
                .map(|codon| -> u8 {
                    match codon {
                        b"ATA" => b'I', b"ATC" => b'I', b"ATT" => b'I', b"ATG" => b'M',
                        b"ACA" => b'T', b"ACC" => b'T', b"ACG" => b'T', b"ACT" => b'T',
                        b"AAC" => b'N', b"AAT" => b'N', b"AAA" => b'K', b"AAG" => b'K',
                        b"AGC" => b'S', b"AGT" => b'S', b"AGA" => b'R', b"AGG" => b'R',
                        b"CTA" => b'L', b"CTC" => b'L', b"CTG" => b'L', b"CTT" => b'L',
                        b"CCA" => b'P', b"CCC" => b'P', b"CCG" => b'P', b"CCT" => b'P',
                        b"CAC" => b'H', b"CAT" => b'H', b"CAA" => b'Q', b"CAG" => b'Q',
                        b"CGA" => b'R', b"CGC" => b'R', b"CGG" => b'R', b"CGT" => b'R',
                        b"GTA" => b'V', b"GTC" => b'V', b"GTG" => b'V', b"GTT" => b'V',
                        b"GCA" => b'A', b"GCC" => b'A', b"GCG" => b'A', b"GCT" => b'A',
                        b"GAC" => b'D', b"GAT" => b'D', b"GAA" => b'E', b"GAG" => b'E',
                        b"GGA" => b'G', b"GGC" => b'G', b"GGG" => b'G', b"GGT" => b'G',
                        b"TCA" => b'S', b"TCC" => b'S', b"TCG" => b'S', b"TCT" => b'S',
                        b"TTC" => b'F', b"TTT" => b'F', b"TTA" => b'L', b"TTG" => b'L',
                        b"TAC" => b'Y', b"TAT" => b'Y', b"TAA" => b'*', b"TAG" => b'?',
                        b"TGC" => b'C', b"TGT" => b'C', b"TGA" => b'*', b"TGG" => b'W',
                        _ => b'X'
                    }
                }).collect_vec();
            
            Sequence::try_from((seq.as_slice(), self.1.as_slice())).unwrap()
        }


        

        fn align(&self, other: &Self) -> Alignment {

            let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

            let mut aligner =
                    Aligner::with_capacity(
                        other.len(), 
                        self.len(), 
                        -5, 
                        -1, 
                        &score
                    );

            aligner.global(other, self)
        }

        pub fn pretty_align(&self, other: &Self) -> (isize, i32) {
            let mut index: isize = 0;

            let alignment = self.align(other);

            for op in alignment.operations.iter() {
                match *op {
                    Match => print!("-"),
                    Subst => print!(
                        "{}",
                        std::str::from_utf8(&[other[index as usize]]).unwrap()
                    ),
                    Del => {
                        index += -1;
                        print!("X");
                    }
                    Ins => {
                        index += -1;
                        print!("^");
                    }
                    _ => {}
                }
                index += 1;
            }

            (index, alignment.score)
        }

        pub fn normal_align(&self, other: &Self) {

            for (one, two) in self.iter().zip(other.0) {

                if *one != two {
                    print!("{}", std::str::from_utf8(&[two]).unwrap());
                } else {
                    print!("-");
                }
            }

        }

        pub fn hamming_distance(&self, other: &Self) -> u8 {

            self
                .iter()
                .zip(other.0)
                .fold(0 as u8, |mut accumulator, (one, two)| {

                    if *one != two {
                        accumulator += 1;
                    }

                    accumulator
                });

            0
        }
    }
}

// Type aliases for clusterings.
pub mod clustering_types {
    #![allow(dead_code)]
    
    // Some more of the same to help represent clusterings.
    // Sometimes I want to cluster based on CDR which has a
    // different type to UMI hence the need for the generic
    // type variables.

    use serde::{Serialize, Deserialize};

    use super::sequence::*;
    use serde_json_any_key::*;

    use std::{
        ops::{Deref, DerefMut},
        collections::HashMap, hash::Hash
    };

    pub type UMItoCDRClustering = Clustering<UMI, CDR>;
    pub type CDRtoUMIClustering = Clustering<CDR, UMI>;


    pub struct Clustering<T : Serialize, U: Serialize >(
       
        pub HashMap<T, HashMap<U, usize>>
    );

    impl<T : Serialize, U: Serialize> Deref for Clustering<T, U> {
        type Target = HashMap<T, HashMap<U, usize>>;
        fn deref(&self) -> &Self::Target {
            &self.0
        }   
    }

    impl<T : Serialize, U: Serialize> DerefMut for Clustering<T, U> {
        fn deref_mut(&mut self) -> &mut Self::Target {
            &mut self.0
        }
    }

    // Helper constructor method.
    impl<T : Serialize, U: Serialize> From<HashMap<T, HashMap<U , usize>>> for Clustering<T, U> {
        fn from(map: HashMap<T, HashMap<U, usize>>) -> Self {
            Self(map)
        }
    }

    impl<T: Serialize, U: Serialize> Clustering<T, U> 
    {
        // If only see the umi < threshold. Remove it.
        pub fn filter_threshold(&mut self, threshold: usize) {
            self.retain(|_, value| value.values().sum::<usize>() > threshold);
        }
    }
    
}

