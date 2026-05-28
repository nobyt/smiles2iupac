"""
Phase 38 テスト: スルフィン酸 (sulfinic acid)

対象 (IUPAC P-65.3.2):
  R-S(=O)-OH パターン。
  命名規則: {stem}anesulfinic acid
  例: CH3-S(=O)-OH → methanesulfinic acid
"""
from smiles2iupac import smiles_to_iupac


class TestSulfinicAcid:

    def test_methanesulfinic_acid(self):
        assert smiles_to_iupac("CS(=O)O") == "methanesulfinic acid"

    def test_ethanesulfinic_acid(self):
        assert smiles_to_iupac("CCS(=O)O") == "ethanesulfinic acid"

    def test_propanesulfinic_acid(self):
        assert smiles_to_iupac("CCCS(=O)O") == "propanesulfinic acid"

    def test_butanesulfinic_acid(self):
        assert smiles_to_iupac("CCCCS(=O)O") == "butanesulfinic acid"


class TestSulfinicAcidVsSulfoxide:

    def test_dimethyl_sulfoxide_unchanged(self):
        # S(=O) + 2C (両側に炭素) → sulfoxide ではなく sulfoxide 命名
        assert smiles_to_iupac("CS(=O)C") == "dimethyl sulfoxide"
