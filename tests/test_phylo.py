import math
import pytest
from biorust import (
    AlignmentDNA,
    AlignmentProtein,
    GappedDNA,
    GappedProtein,
    distance_matrix,
    build_tree,
)


def make_dna_alignment(pairs):
    """Helper: list of (id, seq_str) -> AlignmentDNA."""
    return AlignmentDNA([(id, GappedDNA(seq)) for id, seq in pairs])


def make_protein_alignment(pairs):
    """Helper: list of (id, seq_str) -> AlignmentProtein."""
    return AlignmentProtein([(id, GappedProtein(seq)) for id, seq in pairs])


# ─── distance_matrix ────────────────────────────────────────


class TestDistanceMatrixDNA:
    def test_pdistance_identical(self):
        aln = make_dna_alignment([("a", "ACGT"), ("b", "ACGT")])
        dm = distance_matrix(aln, model="p-distance")
        assert dm.get(0, 1) == pytest.approx(0.0)

    def test_pdistance_known(self):
        aln = make_dna_alignment([("a", "ACGT"), ("b", "ATAT")])
        dm = distance_matrix(aln, model="p-distance")
        assert dm.get(0, 1) == pytest.approx(0.5)

    def test_jc69(self):
        aln = make_dna_alignment([("a", "AAAAAAAAAA"), ("b", "TAAAAAAAAA")])
        dm = distance_matrix(aln, model="jc69")
        p = 0.1
        expected = -0.75 * math.log(1.0 - 4.0 * p / 3.0)
        assert dm.get(0, 1) == pytest.approx(expected)

    def test_k2p(self):
        # 1 transition (A->G), 0 transversions out of 4
        aln = make_dna_alignment([("a", "ACGT"), ("b", "GCGT")])
        dm = distance_matrix(aln, model="k2p")
        p = 0.25
        q = 0.0
        expected = -0.5 * math.log(1.0 - 2.0 * p - q) - 0.25 * math.log(1.0 - 2.0 * q)
        assert dm.get(0, 1) == pytest.approx(expected)

    def test_gap_handling(self):
        aln = make_dna_alignment([("a", "ACGT"), ("b", "A-GT"), ("c", "ACGT")])
        dm = distance_matrix(aln, model="p-distance")
        # a vs b: gap at pos 1 skipped, 3 valid sites, all match -> 0
        assert dm.get(0, 1) == pytest.approx(0.0)
        # a vs c: 4 sites, all match -> 0
        assert dm.get(0, 2) == pytest.approx(0.0)

    def test_invalid_model(self):
        aln = make_dna_alignment([("a", "ACGT"), ("b", "ACGT")])
        with pytest.raises(ValueError, match="unknown DNA distance model"):
            distance_matrix(aln, model="invalid")


class TestDistanceMatrixProtein:
    def test_pdistance(self):
        aln = make_protein_alignment([("a", "ACDE"), ("b", "ACDF")])
        dm = distance_matrix(aln, model="p-distance")
        assert dm.get(0, 1) == pytest.approx(0.25)

    def test_poisson(self):
        aln = make_protein_alignment([("a", "ACDE"), ("b", "ACDF")])
        dm = distance_matrix(aln, model="poisson")
        expected = -math.log(0.75)
        assert dm.get(0, 1) == pytest.approx(expected)

    def test_gap_handling(self):
        aln = make_protein_alignment([("a", "ACD-E"), ("b", "ACDFE")])
        dm = distance_matrix(aln, model="p-distance")
        # gap at pos 3 skipped for pair (a,b), 4 valid, all match -> 0
        assert dm.get(0, 1) == pytest.approx(0.0)

    def test_invalid_model(self):
        aln = make_protein_alignment([("a", "ACDE"), ("b", "ACDF")])
        with pytest.raises(ValueError, match="unknown protein distance model"):
            distance_matrix(aln, model="jc69")


# ─── DistanceMatrix accessors ───────────────────────────────


class TestDistanceMatrixAccessors:
    def test_n_and_len(self):
        aln = make_dna_alignment([("a", "AA"), ("b", "AA"), ("c", "AA")])
        dm = distance_matrix(aln)
        assert dm.n == 3
        assert len(dm) == 3

    def test_labels(self):
        aln = make_dna_alignment([("seq1", "ACGT"), ("seq2", "ACGT")])
        dm = distance_matrix(aln)
        assert dm.labels() == ["seq1", "seq2"]

    def test_to_list(self):
        aln = make_dna_alignment([("a", "AAAA"), ("b", "AAAT")])
        dm = distance_matrix(aln)
        lst = dm.to_list()
        assert len(lst) == 4  # 2x2
        assert lst[0] == pytest.approx(0.0)  # a-a
        assert lst[1] == pytest.approx(0.25)  # a-b
        assert lst[2] == pytest.approx(0.25)  # b-a
        assert lst[3] == pytest.approx(0.0)  # b-b

    def test_to_list_of_lists(self):
        aln = make_dna_alignment([("a", "AAAA"), ("b", "AAAT")])
        dm = distance_matrix(aln)
        lol = dm.to_list_of_lists()
        assert len(lol) == 2
        assert len(lol[0]) == 2
        assert lol[0][1] == pytest.approx(0.25)

    def test_repr(self):
        aln = make_dna_alignment([("a", "ACGT"), ("b", "ACGT")])
        dm = distance_matrix(aln)
        assert "DistanceMatrix" in repr(dm)

    def test_str(self):
        aln = make_dna_alignment([("a", "ACGT"), ("b", "ACGT")])
        dm = distance_matrix(aln)
        s = str(dm)
        assert "a" in s
        assert "b" in s

    def test_get_out_of_range(self):
        aln = make_dna_alignment([("a", "ACGT"), ("b", "ACGT")])
        dm = distance_matrix(aln)
        with pytest.raises(ValueError):
            dm.get(5, 0)


# ─── build_tree ──────────────────────────────────────────────


class TestBuildTree:
    def _simple_dm(self):
        aln = make_dna_alignment(
            [("A", "AAAA"), ("B", "AAAT"), ("C", "AATT"), ("D", "ATTT")]
        )
        return distance_matrix(aln)

    def test_nj_unrooted(self):
        dm = self._simple_dm()
        tree = build_tree(dm, method="nj")
        assert not tree.is_rooted()
        assert tree.num_leaves() == 4
        assert set(tree.leaf_labels()) == {"A", "B", "C", "D"}

    def test_upgma_rooted(self):
        dm = self._simple_dm()
        tree = build_tree(dm, method="upgma")
        assert tree.is_rooted()
        assert tree.num_leaves() == 4
        assert set(tree.leaf_labels()) == {"A", "B", "C", "D"}

    def test_invalid_method(self):
        dm = self._simple_dm()
        with pytest.raises(ValueError, match="unknown tree method"):
            build_tree(dm, method="parsimony")

    def test_default_method_is_nj(self):
        dm = self._simple_dm()
        tree = build_tree(dm)
        assert not tree.is_rooted()

    def test_newick_ends_with_semicolon(self):
        dm = self._simple_dm()
        tree = build_tree(dm)
        nwk = tree.to_newick()
        assert nwk.startswith("(")
        assert nwk.endswith(";")

    def test_newick_contains_labels(self):
        dm = self._simple_dm()
        tree = build_tree(dm)
        nwk = tree.to_newick()
        for label in ["A", "B", "C", "D"]:
            assert label in nwk

    def test_tree_repr(self):
        dm = self._simple_dm()
        tree = build_tree(dm, method="nj")
        assert "PhyloTree" in repr(tree)
        assert "unrooted" in repr(tree)

    def test_tree_str_is_newick(self):
        dm = self._simple_dm()
        tree = build_tree(dm)
        assert str(tree) == tree.to_newick()

    def test_two_taxa(self):
        aln = make_dna_alignment([("X", "AAAA"), ("Y", "TTTT")])
        dm = distance_matrix(aln)
        tree = build_tree(dm, method="nj")
        assert tree.num_leaves() == 2
        assert set(tree.leaf_labels()) == {"X", "Y"}

    def test_ascii_diagram(self):
        dm = self._simple_dm()
        tree = build_tree(dm, method="nj")
        diagram = tree.ascii_diagram()
        assert "A" in diagram
        assert "B" in diagram
        assert "C" in diagram
        assert "D" in diagram
        assert "|--" in diagram or "`--" in diagram


class TestAlignmentDiagram:
    def test_dot_gap_not_conserved(self):
        aln = AlignmentDNA([("a", GappedDNA("A.C")), ("b", GappedDNA("A-C"))])
        diagram = aln.alignment_diagram()
        cons_line = diagram.splitlines()[-1]
        cons = cons_line[-aln.width :]
        assert cons == "* *"

    def test_empty_slice(self):
        aln = AlignmentDNA([("a", GappedDNA("ACGT")), ("b", GappedDNA("ACGT"))])
        empty = aln[:0]
        assert len(empty) == 0
        assert empty.width == 0
        assert empty.aligned_strings() == []
        assert str(empty) == ""
