"""
Phase 51 テスト: ジエステル (diester)

対象 (IUPAC P-65.1.1.4.2):
  R-OC(=O)-chain-C(=O)OR パターン → {dialkyl} {stem}anedioate
  例: CCOC(=O)CC(=O)OCC → diethyl propanedioate
"""
from smiles2iupac import smiles_to_iupac


class TestSymmetricDiester:

    def test_dimethyl_propanedioate(self):
        assert smiles_to_iupac("COC(=O)CC(=O)OC") == "dimethyl propanedioate"

    def test_diethyl_propanedioate(self):
        assert smiles_to_iupac("CCOC(=O)CC(=O)OCC") == "diethyl propanedioate"

    def test_diethyl_butanedioate(self):
        assert smiles_to_iupac("CCOC(=O)CCC(=O)OCC") == "diethyl butanedioate"

    def test_dimethyl_pentanedioate(self):
        assert smiles_to_iupac("COC(=O)CCCC(=O)OC") == "dimethyl pentanedioate"


class TestDiesterVsMonoester:

    def test_monoester_unchanged(self):
        assert smiles_to_iupac("CCOC(=O)C") == "ethyl acetate"

    def test_monoester_propanoate_unchanged(self):
        assert smiles_to_iupac("CCOC(=O)CC") == "ethyl propanoate"
