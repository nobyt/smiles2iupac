"""
Phase 44 テスト: ヒドロペルオキシド (hydroperoxide)

対象 (IUPAC P-65.1.4.3):
  R-O-O-H パターン（過酸化水素の一置換体）。
  命名規則 (functional class): {alkyl} hydroperoxide
  例: CH₃-O-O-H → methyl hydroperoxide
      CH₃CH₂-O-O-H → ethyl hydroperoxide
"""
from smiles2iupac import smiles_to_iupac


class TestHydroperoxide:

    def test_methyl_hydroperoxide(self):
        assert smiles_to_iupac("COO") == "methyl hydroperoxide"

    def test_methyl_hydroperoxide_reversed(self):
        assert smiles_to_iupac("OOC") == "methyl hydroperoxide"

    def test_ethyl_hydroperoxide(self):
        assert smiles_to_iupac("CCOO") == "ethyl hydroperoxide"

    def test_propyl_hydroperoxide(self):
        assert smiles_to_iupac("CCCOO") == "propyl hydroperoxide"


class TestHydroperoxideVsAlcohol:

    def test_alcohol_unchanged(self):
        # アルコール (O-H, C-O-H) は変わらない → ethanol (IUPAC 保留名)
        assert smiles_to_iupac("CCO") == "ethanol"

    def test_hydroperoxide_not_oxy(self):
        result = smiles_to_iupac("CCOO")
        assert result == "ethyl hydroperoxide"
        assert "oxy" not in result
