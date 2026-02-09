pub mod distance;
pub mod newick;
pub mod tree;

pub use distance::{
    dna_distance_matrix, protein_distance_matrix, DistanceMatrix, DnaDistanceModel,
    ProteinDistanceModel,
};
pub use newick::to_newick;
pub use tree::{neighbor_joining, upgma, PhyloNode, PhyloTree};

#[cfg(test)]
mod tests;
