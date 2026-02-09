use super::*;

fn labels(names: &[&str]) -> Vec<Box<str>> {
    names
        .iter()
        .map(|s| s.to_string().into_boxed_str())
        .collect()
}

// ─── p-distance ─────────────────────────────────────────────

#[test]
fn pdist_identical() {
    let seqs: Vec<&[u8]> = vec![b"ACGT", b"ACGT"];
    let dm = dna_distance_matrix(&seqs, labels(&["a", "b"]), DnaDistanceModel::PDistance).unwrap();
    assert_eq!(dm.get(0, 1), 0.0);
    assert_eq!(dm.get(1, 0), 0.0);
    assert_eq!(dm.get(0, 0), 0.0);
}

#[test]
fn pdist_known() {
    // 2 mismatches out of 4
    let seqs: Vec<&[u8]> = vec![b"ACGT", b"ATAT"];
    let dm = dna_distance_matrix(&seqs, labels(&["a", "b"]), DnaDistanceModel::PDistance).unwrap();
    assert!((dm.get(0, 1) - 0.5).abs() < 1e-10);
}

#[test]
fn pdist_three_seqs() {
    let seqs: Vec<&[u8]> = vec![b"AAAA", b"AAAT", b"AATT"];
    let dm =
        dna_distance_matrix(&seqs, labels(&["a", "b", "c"]), DnaDistanceModel::PDistance).unwrap();
    assert!((dm.get(0, 1) - 0.25).abs() < 1e-10);
    assert!((dm.get(0, 2) - 0.50).abs() < 1e-10);
    assert!((dm.get(1, 2) - 0.25).abs() < 1e-10);
}

// ─── gap pairwise deletion ──────────────────────────────────

#[test]
fn gap_pairwise_deletion() {
    // Position 1 has gap in seq b, should be skipped for (a,b) pair
    let seqs: Vec<&[u8]> = vec![b"ACGT", b"A-GT", b"ACGT"];
    let dm =
        dna_distance_matrix(&seqs, labels(&["a", "b", "c"]), DnaDistanceModel::PDistance).unwrap();
    // a vs b: compare positions 0,2,3 -> all match -> 0.0
    assert!((dm.get(0, 1) - 0.0).abs() < 1e-10);
    // a vs c: all 4 positions valid, all match
    assert!((dm.get(0, 2) - 0.0).abs() < 1e-10);
}

#[test]
fn gap_dot_is_gap() {
    let seqs: Vec<&[u8]> = vec![b"ACGT", b"A.GT"];
    let dm = dna_distance_matrix(&seqs, labels(&["a", "b"]), DnaDistanceModel::PDistance).unwrap();
    assert!((dm.get(0, 1) - 0.0).abs() < 1e-10);
}

// ─── Jukes-Cantor ───────────────────────────────────────────

#[test]
fn jc69_small() {
    // p = 0.1 -> JC = -3/4 * ln(1 - 4*0.1/3)
    let seqs: Vec<&[u8]> = vec![b"AAAAAAAAAA", b"TAAAAAAAAA"];
    let dm =
        dna_distance_matrix(&seqs, labels(&["a", "b"]), DnaDistanceModel::JukesCantor).unwrap();
    let p = 0.1;
    let expected = -0.75 * (1.0 - 4.0 * p / 3.0_f64).ln();
    assert!((dm.get(0, 1) - expected).abs() < 1e-10);
}

#[test]
fn jc69_saturated() {
    // p = 0.8 >= 3/4 -> should fail
    let seqs: Vec<&[u8]> = vec![b"AAAAA", b"TTTTT"];
    let result = dna_distance_matrix(&seqs, labels(&["a", "b"]), DnaDistanceModel::JukesCantor);
    assert!(result.is_err());
}

// ─── Kimura 2-parameter ─────────────────────────────────────

#[test]
fn k2p_known() {
    // 1 transition (A->G), 0 transversions out of 4 valid
    let seqs: Vec<&[u8]> = vec![b"ACGT", b"GCGT"];
    let dm = dna_distance_matrix(&seqs, labels(&["a", "b"]), DnaDistanceModel::Kimura2P).unwrap();
    let p: f64 = 0.25; // transition proportion
    let q: f64 = 0.0; // transversion proportion
    let expected = -0.5 * (1.0 - 2.0 * p - q).ln() - 0.25 * (1.0 - 2.0 * q).ln();
    assert!((dm.get(0, 1) - expected).abs() < 1e-10);
}

#[test]
fn k2p_skips_ambiguity() {
    // N is not ACGT, so that column is skipped
    let seqs: Vec<&[u8]> = vec![b"ACGTN", b"GCGTN"];
    let dm = dna_distance_matrix(&seqs, labels(&["a", "b"]), DnaDistanceModel::Kimura2P).unwrap();
    // Same as without the N column
    let seqs2: Vec<&[u8]> = vec![b"ACGT", b"GCGT"];
    let dm2 = dna_distance_matrix(&seqs2, labels(&["a", "b"]), DnaDistanceModel::Kimura2P).unwrap();
    assert!((dm.get(0, 1) - dm2.get(0, 1)).abs() < 1e-10);
}

// ─── Protein distances ──────────────────────────────────────

#[test]
fn protein_pdist() {
    let seqs: Vec<&[u8]> = vec![b"ACDE", b"ACDF"];
    let dm = protein_distance_matrix(&seqs, labels(&["a", "b"]), ProteinDistanceModel::PDistance)
        .unwrap();
    assert!((dm.get(0, 1) - 0.25).abs() < 1e-10);
}

#[test]
fn poisson_known() {
    // p = 0.25 -> -ln(0.75)
    let seqs: Vec<&[u8]> = vec![b"ACDE", b"ACDF"];
    let dm =
        protein_distance_matrix(&seqs, labels(&["a", "b"]), ProteinDistanceModel::Poisson).unwrap();
    let expected = -(0.75_f64).ln();
    assert!((dm.get(0, 1) - expected).abs() < 1e-10);
}

#[test]
fn poisson_saturated() {
    // All different
    let seqs: Vec<&[u8]> = vec![b"ACDE", b"FGHI"];
    let result = protein_distance_matrix(&seqs, labels(&["a", "b"]), ProteinDistanceModel::Poisson);
    assert!(result.is_err());
}

// ─── no valid sites ─────────────────────────────────────────

#[test]
fn no_valid_sites() {
    let seqs: Vec<&[u8]> = vec![b"----", b"ACGT"];
    let result = dna_distance_matrix(&seqs, labels(&["a", "b"]), DnaDistanceModel::PDistance);
    assert!(result.is_err());
}

// ─── too few sequences ──────────────────────────────────────

#[test]
fn too_few_seqs() {
    let seqs: Vec<&[u8]> = vec![b"ACGT"];
    let result = dna_distance_matrix(&seqs, labels(&["a"]), DnaDistanceModel::PDistance);
    assert!(result.is_err());
}

#[test]
fn label_count_mismatch() {
    let seqs: Vec<&[u8]> = vec![b"ACGT", b"ACGT"];
    let result = dna_distance_matrix(&seqs, labels(&["a"]), DnaDistanceModel::PDistance);
    assert!(result.is_err());
}

#[test]
fn unequal_sequence_lengths() {
    let seqs: Vec<&[u8]> = vec![b"ACGT", b"ACG"];
    let result = dna_distance_matrix(&seqs, labels(&["a", "b"]), DnaDistanceModel::PDistance);
    assert!(result.is_err());
}

// ─── NJ tree ────────────────────────────────────────────────

fn simple_4taxa_dm() -> DistanceMatrix {
    // Additive distance matrix for tree: ((A:1,B:1):1,(C:1,D:1):1)
    // A-B=2, A-C=4, A-D=4, B-C=4, B-D=4, C-D=2
    let labels = labels(&["A", "B", "C", "D"]);
    let data = vec![
        0.0, 2.0, 4.0, 4.0, //
        2.0, 0.0, 4.0, 4.0, //
        4.0, 4.0, 0.0, 2.0, //
        4.0, 4.0, 2.0, 0.0, //
    ];
    DistanceMatrix::new(labels, data)
}

#[test]
fn nj_basic_topology() {
    let dm = simple_4taxa_dm();
    let tree = neighbor_joining(&dm).unwrap();
    assert_eq!(tree.num_leaves(), 4);
    assert!(tree.root().is_none()); // NJ is unrooted
    let ll = tree.leaf_labels();
    assert!(ll.contains(&"A".to_string()));
    assert!(ll.contains(&"B".to_string()));
    assert!(ll.contains(&"C".to_string()));
    assert!(ll.contains(&"D".to_string()));
}

#[test]
fn nj_branch_lengths() {
    let dm = simple_4taxa_dm();
    let tree = neighbor_joining(&dm).unwrap();
    // For additive tree ((A:1,B:1):1,(C:1,D:1):1), NJ should recover branch lengths
    // Leaves should have branch length 1.0
    for leaf in tree.leaves() {
        let bl = tree.node(leaf).branch_length.unwrap();
        assert!(
            (bl - 1.0).abs() < 1e-10,
            "leaf {} has branch_length {}",
            tree.node(leaf).label.as_deref().unwrap_or("?"),
            bl
        );
    }
}

#[test]
fn nj_two_taxa() {
    let labels = labels(&["X", "Y"]);
    let data = vec![0.0, 3.0, 3.0, 0.0];
    let dm = DistanceMatrix::new(labels, data);
    let tree = neighbor_joining(&dm).unwrap();
    assert_eq!(tree.num_leaves(), 2);
    // Each leaf should have branch length 1.5
    for leaf in tree.leaves() {
        let bl = tree.node(leaf).branch_length.unwrap();
        assert!((bl - 1.5).abs() < 1e-10);
    }
}

// ─── UPGMA tree ─────────────────────────────────────────────

#[test]
fn upgma_ultrametric() {
    // Ultrametric distances: ((A:1,B:1):1,(C:2,D:2))
    // A-B=2, A-C=4, A-D=4, B-C=4, B-D=4, C-D=4
    let labels = labels(&["A", "B", "C", "D"]);
    let data = vec![
        0.0, 2.0, 4.0, 4.0, //
        2.0, 0.0, 4.0, 4.0, //
        4.0, 4.0, 0.0, 4.0, //
        4.0, 4.0, 4.0, 0.0, //
    ];
    let dm = DistanceMatrix::new(labels, data);
    let tree = upgma(&dm).unwrap();

    assert_eq!(tree.num_leaves(), 4);
    assert!(tree.root().is_some()); // UPGMA is rooted
}

#[test]
fn upgma_two_taxa() {
    let labels = labels(&["X", "Y"]);
    let data = vec![0.0, 6.0, 6.0, 0.0];
    let dm = DistanceMatrix::new(labels, data);
    let tree = upgma(&dm).unwrap();
    assert_eq!(tree.num_leaves(), 2);
    assert!(tree.root().is_some());
    // Each leaf should have branch length 3.0
    for leaf in tree.leaves() {
        let bl = tree.node(leaf).branch_length.unwrap();
        assert!((bl - 3.0).abs() < 1e-10);
    }
}

// ─── Newick ─────────────────────────────────────────────────

#[test]
fn newick_format() {
    let dm = simple_4taxa_dm();
    let tree = neighbor_joining(&dm).unwrap();
    let nwk = to_newick(&tree);
    assert!(nwk.starts_with('('));
    assert!(nwk.ends_with(';'));
    // Should contain all leaf labels
    assert!(nwk.contains('A'));
    assert!(nwk.contains('B'));
    assert!(nwk.contains('C'));
    assert!(nwk.contains('D'));
}

#[test]
fn newick_quotes_labels() {
    let labels = labels(&["A B", "C:D", "E'F", "G"]);
    let data = vec![
        0.0, 1.0, 2.0, 3.0, //
        1.0, 0.0, 2.0, 3.0, //
        2.0, 2.0, 0.0, 3.0, //
        3.0, 3.0, 3.0, 0.0, //
    ];
    let dm = DistanceMatrix::new(labels, data);
    let tree = neighbor_joining(&dm).unwrap();
    let nwk = to_newick(&tree);
    assert!(nwk.contains("'A B'"));
    assert!(nwk.contains("'C:D'"));
    assert!(nwk.contains("'E''F'"));
}

#[test]
fn newick_upgma() {
    let labels = labels(&["X", "Y"]);
    let data = vec![0.0, 4.0, 4.0, 0.0];
    let dm = DistanceMatrix::new(labels, data);
    let tree = upgma(&dm).unwrap();
    let nwk = to_newick(&tree);
    assert!(nwk.starts_with('('));
    assert!(nwk.ends_with(';'));
    assert!(nwk.contains('X'));
    assert!(nwk.contains('Y'));
}

// ─── DistanceMatrix accessors ───────────────────────────────

#[test]
fn dm_accessors() {
    let labels = labels(&["a", "b"]);
    let data = vec![0.0, 1.5, 1.5, 0.0];
    let dm = DistanceMatrix::new(labels, data);
    assert_eq!(dm.n(), 2);
    assert_eq!(dm.labels().len(), 2);
    assert_eq!(dm.data().len(), 4);
    assert!((dm.get(0, 1) - 1.5).abs() < 1e-10);
}

#[test]
fn dm_set_symmetric() {
    let labels = labels(&["a", "b", "c"]);
    let data = vec![0.0; 9];
    let mut dm = DistanceMatrix::new(labels, data);
    dm.set(0, 2, 5.0);
    assert!((dm.get(0, 2) - 5.0).abs() < 1e-10);
    assert!((dm.get(2, 0) - 5.0).abs() < 1e-10);
}

// ─── node counts ────────────────────────────────────────────

#[test]
fn nj_node_count() {
    let dm = simple_4taxa_dm();
    let tree = neighbor_joining(&dm).unwrap();
    // 4 leaves + (4-2) internal + 1 pseudo-root = 4+2+1 = 7? No:
    // NJ: n leaves + (n-2) internal from main loop + 1 final = 2n-1 for n>2
    // For n=4: 4 leaves + 2 internal + 1 final = 7
    // But arena is 2n-2 = 6, plus one final -> 7
    assert_eq!(tree.num_leaves(), 4);
    assert_eq!(tree.num_nodes(), 7); // 4 + 3 internal (2 from loop + 1 final)
}

#[test]
fn upgma_node_count() {
    let labels = labels(&["A", "B", "C"]);
    let data = vec![
        0.0, 2.0, 4.0, //
        2.0, 0.0, 4.0, //
        4.0, 4.0, 0.0, //
    ];
    let dm = DistanceMatrix::new(labels, data);
    let tree = upgma(&dm).unwrap();
    // 3 leaves + 2 internal = 5
    assert_eq!(tree.num_nodes(), 5);
    assert_eq!(tree.num_leaves(), 3);
}
