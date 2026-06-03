"""
Phase 51 テスト: ジエステル (diester)

対象 (IUPAC P-65.1.1.4.2):
  R-OC(=O)-chain-C(=O)OR パターン → {dialkyl} {retained} (retained name preferred)
  例: CCOC(=O)CC(=O)OCC → diethyl malonate
"""
from smiles2iupac import smiles_to_iupac


class TestSymmetricDiester:

    def test_dimethyl_malonate(self):
        assert smiles_to_iupac("COC(=O)CC(=O)OC") == "dimethyl malonate"

    def test_diethyl_malonate(self):
        assert smiles_to_iupac("CCOC(=O)CC(=O)OCC") == "diethyl malonate"

    def test_diethyl_succinate(self):
        assert smiles_to_iupac("CCOC(=O)CCC(=O)OCC") == "diethyl succinate"

    def test_dimethyl_glutarate(self):
        assert smiles_to_iupac("COC(=O)CCCC(=O)OC") == "dimethyl glutarate"


class TestDiesterVsMonoester:

    def test_monoester_unchanged(self):
        assert smiles_to_iupac("CCOC(=O)C") == "ethyl acetate"

    def test_monoester_propanoate_unchanged(self):
        assert smiles_to_iupac("CCOC(=O)CC") == "ethyl propanoate"
