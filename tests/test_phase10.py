"""
Phase 10 テスト: エステル・酸ハライド

対象:
  - エステル (-COOR) → functional class name: (alkyl) (acid)oate
  - 酸ハライド (-COX) → functional class name: (acid)oyl (halide)
"""
from smiles2iupac import smiles_to_iupac


# ─── エステル ────────────────────────────────────────────────────────────

class TestEsters:

    def test_methyl_ethanoate(self):
        # CH₃COOCH₃ (酢酸メチル)
        assert smiles_to_iupac("COC(=O)C") == "methyl acetate"

    def test_ethyl_ethanoate(self):
        # CH₃COOC₂H₅ (酢酸エチル)
        assert smiles_to_iupac("CCOC(=O)C") == "ethyl acetate"

    def test_methyl_propanoate(self):
        # CH₃CH₂COOCH₃
        assert smiles_to_iupac("COC(=O)CC") == "methyl propanoate"

    def test_ethyl_propanoate(self):
        assert smiles_to_iupac("CCOC(=O)CC") == "ethyl propanoate"

    def test_methyl_methanoate(self):
        # HCOOCH₃ (ギ酸メチル)
        assert smiles_to_iupac("COC=O") == "methyl formate"


# ─── 酸ハライド ──────────────────────────────────────────────────────────

class TestAcidHalides:

    def test_ethanoyl_chloride(self):
        # CH₃COCl (塩化アセチル)
        assert smiles_to_iupac("CC(=O)Cl") == "acetyl chloride"

    def test_acetyl_bromide(self):
        assert smiles_to_iupac("CC(=O)Br") == "acetyl bromide"

    def test_propanoyl_chloride(self):
        assert smiles_to_iupac("CCC(=O)Cl") == "propanoyl chloride"

    def test_methanoyl_chloride(self):
        # HCOCl (塩化ホルミル)
        assert smiles_to_iupac("ClC=O") == "formyl chloride"

    def test_butanoyl_chloride(self):
        assert smiles_to_iupac("CCCC(=O)Cl") == "butanoyl chloride"
