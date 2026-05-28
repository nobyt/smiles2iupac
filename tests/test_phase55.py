"""
Phase 55 テスト: チオエステル (thioester)

対象 (IUPAC P-65.1.2.3):
  R-C(=O)-S-R' パターン → S-{alkyl} {acid}thioate
  例: CCSC(=O)C → S-ethyl ethanethioate
"""
from smiles2iupac import smiles_to_iupac


class TestThioester:

    def test_s_methyl_ethanethioate(self):
        assert smiles_to_iupac("CSC(=O)C") == "S-methyl ethanethioate"

    def test_s_ethyl_ethanethioate(self):
        assert smiles_to_iupac("CCSC(=O)C") == "S-ethyl ethanethioate"

    def test_s_methyl_propanethioate(self):
        assert smiles_to_iupac("CSC(=O)CC") == "S-methyl propanethioate"


class TestThioesterVsSulfide:

    def test_dimethyl_sulfide_unchanged(self):
        assert smiles_to_iupac("CSC") == "dimethyl sulfide"

    def test_thioester_not_sulfide(self):
        result = smiles_to_iupac("CCSC(=O)C")
        assert "sulfide" not in result
        assert "thioate" in result
