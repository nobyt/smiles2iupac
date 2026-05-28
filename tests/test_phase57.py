"""
Phase 57 テスト: 有機ペルオキシド (organic peroxide)

対象 (IUPAC P-65.1.4.4):
  R-O-O-R' パターン (H なし) → {dialkyl} peroxide
  例: CCOOCC → diethyl peroxide
"""
from smiles2iupac import smiles_to_iupac


class TestPeroxide:

    def test_dimethyl_peroxide(self):
        assert smiles_to_iupac("COOC") == "dimethyl peroxide"

    def test_diethyl_peroxide(self):
        assert smiles_to_iupac("CCOOCC") == "diethyl peroxide"

    def test_ethyl_methyl_peroxide(self):
        assert smiles_to_iupac("CCOOC") == "ethyl methyl peroxide"


class TestPeroxideVsHydroperoxide:

    def test_hydroperoxide_unchanged(self):
        assert smiles_to_iupac("CCOO") == "ethyl hydroperoxide"

    def test_methyl_hydroperoxide_unchanged(self):
        assert smiles_to_iupac("COO") == "methyl hydroperoxide"
