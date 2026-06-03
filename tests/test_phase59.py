"""
Phase 59 テスト: スルホニルクロライド (sulfonyl chloride)

対象 (IUPAC P-65.3.1):
  R-S(=O)₂-Cl パターン → {alkane}sulfonyl chloride
  例: CS(=O)(=O)Cl → methanesulfonyl chloride
"""
from smiles2iupac import smiles_to_iupac


class TestSulfonylChloride:

    def test_methanesulfonyl_chloride(self):
        assert smiles_to_iupac("CS(=O)(=O)Cl") == "methanesulfonyl chloride"

    def test_ethanesulfonyl_chloride(self):
        assert smiles_to_iupac("CCS(=O)(=O)Cl") == "ethanesulfonyl chloride"

    def test_propanesulfonyl_chloride(self):
        assert smiles_to_iupac("CCCS(=O)(=O)Cl") == "propane-1-sulfonyl chloride"


class TestSulfonylChlorideVsSulfonate:

    def test_sulfonic_acid_unchanged(self):
        assert smiles_to_iupac("CS(=O)(=O)O") == "methanesulfonic acid"

    def test_sulfonamide_unchanged(self):
        assert smiles_to_iupac("CS(=O)(=O)N") == "methanesulfonamide"
