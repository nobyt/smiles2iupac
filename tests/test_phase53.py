"""
Phase 53 テスト: アジド (azide)

対象 (IUPAC P-68.3.2.3):
  R-N=[N+]=[N-] パターン。命名: {alkyl} azide (functional class)
  例: CCN=[N+]=[N-] → ethyl azide
"""
from smiles2iupac import smiles_to_iupac


class TestAzide:

    def test_methyl_azide(self):
        assert smiles_to_iupac("CN=[N+]=[N-]") == "methyl azide"

    def test_ethyl_azide(self):
        assert smiles_to_iupac("CCN=[N+]=[N-]") == "ethyl azide"

    def test_propyl_azide(self):
        assert smiles_to_iupac("CCCN=[N+]=[N-]") == "propyl azide"


class TestAzideVsAmine:

    def test_amine_unchanged(self):
        assert smiles_to_iupac("CCN") == "ethanamine"

    def test_azide_not_amine(self):
        result = smiles_to_iupac("CCN=[N+]=[N-]")
        assert "amine" not in result
        assert "azide" in result
