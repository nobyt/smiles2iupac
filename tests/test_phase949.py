"""Phase 949: isotopic labeling (IUPAC 2013 P-82) on simple heteroaryl substituents.

Phase 939 extended isotope labeling to phenyl/naphthyl substituents, but the
separate _name_aryl_substituent code path used for simple retained-name
heteroaromatic substituents (pyridyl, furyl, thiophenyl, ...) had no
format_isotope_descriptor call at all -- isotope labels on such substituents
were silently and completely dropped.

All expected names below are OPSIN round-trip verified.
"""
from smiles2iupac import smiles_to_iupac


class TestPhase949IsotopePyridylSubstituent:
    def test_unlabeled_pyridyl_unaffected(self):
        assert smiles_to_iupac("CC(=O)Nc1ccncc1") == "N-(pyridin-4-yl)acetamide"

    def test_meta_deuterated_pyridyl(self):
        assert smiles_to_iupac("CC(=O)Nc1cc([2H])ncc1") == "N-[(6-2H)pyridin-4-yl]acetamide"

    def test_ortho_to_n_deuterated_pyridyl(self):
        assert smiles_to_iupac("CC(=O)Nc1ccnc([2H])c1") == "N-[(2-2H)pyridin-4-yl]acetamide"


class TestPhase949IsotopeFurylThiophenylSubstituent:
    def test_unlabeled_furyl_unaffected(self):
        assert smiles_to_iupac("CC(=O)Nc1ccoc1") == "N-(furan-3-yl)acetamide"

    def test_deuterated_furyl(self):
        assert smiles_to_iupac("CC(=O)Nc1cc([2H])oc1") == "N-[(5-2H)furan-3-yl]acetamide"

    def test_unlabeled_thiophenyl_unaffected(self):
        assert smiles_to_iupac("CC(=O)Nc1cccs1") == "N-(thiophen-2-yl)acetamide"

    def test_deuterated_thiophenyl(self):
        assert smiles_to_iupac("CC(=O)Nc1ccc([2H])s1") == "N-[(5-2H)thiophen-2-yl]acetamide"
