"""
Phase 56 テスト: ジスルフィド (disulfide)

対象 (IUPAC P-63.7.1):
  R-S-S-R' パターン → {dialkyl} disulfide
  例: CSSC → (methyldisulfanyl)methane
"""
from smiles2iupac import smiles_to_iupac


class TestDisulfide:

    def test_dimethyl_disulfide(self):
        assert smiles_to_iupac("CSSC") == "(methyldisulfanyl)methane"

    def test_diethyl_disulfide(self):
        assert smiles_to_iupac("CCSSCC") == "(ethyldisulfanyl)ethane"


class TestDisulfideVsSulfide:

    def test_sulfide_unchanged(self):
        assert smiles_to_iupac("CSC") == "(methylsulfanyl)methane"

    def test_thiol_unchanged(self):
        assert smiles_to_iupac("CCS") == "ethanethiol"
