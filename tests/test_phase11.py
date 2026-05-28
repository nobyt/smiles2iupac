"""
Phase 11 テスト: アミド (amide)

対象:
  - 第一級アミド (-CONH₂) → -amide suffix
  - 鎖状・環状・芳香族
"""
from smiles2iupac import smiles_to_iupac


class TestChainAmides:

    def test_formamide(self):
        # HCONH₂
        assert smiles_to_iupac("NC=O") == "formamide"

    def test_acetamide(self):
        # CH₃CONH₂ (acetamide)
        assert smiles_to_iupac("CC(=O)N") == "acetamide"

    def test_propanamide(self):
        assert smiles_to_iupac("CCC(=O)N") == "propanamide"

    def test_butanamide(self):
        assert smiles_to_iupac("CCCC(=O)N") == "butanamide"

    def test_2_methylpropanamide(self):
        # (CH₃)₂CHCONH₂
        assert smiles_to_iupac("CC(C)C(=O)N") == "2-methylpropanamide"


class TestAromaticAmides:

    def test_benzamide(self):
        # Ph-CONH₂
        assert smiles_to_iupac("NC(=O)c1ccccc1") == "benzamide"
