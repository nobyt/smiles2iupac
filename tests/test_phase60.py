"""
Phase 60 テスト: クロロホルメート (chloroformate)

対象 (IUPAC P-65.1.2.6):
  Cl-C(=O)-O-R パターン → {alkyl} carbonochloridate
  例: ClC(=O)OCC → ethyl carbonochloridate
"""
from smiles2iupac import smiles_to_iupac


class TestChloroformate:

    def test_methyl_carbonochloridate(self):
        assert smiles_to_iupac("ClC(=O)OC") == "methyl carbonochloridate"

    def test_ethyl_carbonochloridate(self):
        assert smiles_to_iupac("ClC(=O)OCC") == "ethyl carbonochloridate"

    def test_propyl_carbonochloridate(self):
        assert smiles_to_iupac("ClC(=O)OCCC") == "propyl carbonochloridate"


class TestChloroformateVsEster:

    def test_ester_unchanged(self):
        assert smiles_to_iupac("CC(=O)OCC") == "ethyl acetate"

    def test_acid_chloride_unchanged(self):
        assert smiles_to_iupac("CC(=O)Cl") == "acetyl chloride"
