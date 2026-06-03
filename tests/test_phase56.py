"""
Phase 56 テスト: ジスルフィド (disulfide)

対象 (IUPAC P-63.7.1):
  R-S-S-R' パターン → {dialkyl} disulfide
  例: CSSC → dimethyl disulfide
"""
from smiles2iupac import smiles_to_iupac


class TestDisulfide:

    def test_dimethyl_disulfide(self):
        assert smiles_to_iupac("CSSC") == "dimethyl disulfide"

    def test_diethyl_disulfide(self):
        assert smiles_to_iupac("CCSSCC") == "diethyl disulfide"


class TestDisulfideVsSulfide:

    def test_sulfide_unchanged(self):
        assert smiles_to_iupac("CSC") == "dimethyl sulfide"

    def test_thiol_unchanged(self):
        assert smiles_to_iupac("CCS") == "ethanethiol"
