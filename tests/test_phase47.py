"""
Phase 47 テスト: 芳香族酸エステル (aryl acid ester)

対象 (IUPAC P-65.1.1.4):
  R-O-C(=O)-Ph パターン → {alkyl} benzoate
  例: CCOC(=O)c1ccccc1 → ethyl benzoate
"""
from smiles2iupac import smiles_to_iupac


class TestBenzoate:

    def test_methyl_benzoate(self):
        assert smiles_to_iupac("COC(=O)c1ccccc1") == "methyl benzoate"

    def test_ethyl_benzoate(self):
        assert smiles_to_iupac("CCOC(=O)c1ccccc1") == "ethyl benzoate"

    def test_propyl_benzoate(self):
        assert smiles_to_iupac("CCCOC(=O)c1ccccc1") == "propyl benzoate"


class TestBenzoateVsAlkanoate:

    def test_ethyl_ethanoate_unchanged(self):
        assert smiles_to_iupac("CCOC(=O)C") == "ethyl acetate"

    def test_ethyl_propanoate_unchanged(self):
        assert smiles_to_iupac("CCOC(=O)CC") == "ethyl propanoate"
