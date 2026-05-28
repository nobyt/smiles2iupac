"""
Phase 54 テスト: イソシアネート (isocyanate)

対象 (IUPAC P-65.2.2):
  R-N=C=O パターン。命名: {alkyl} isocyanate (functional class)
  例: CCN=C=O → ethyl isocyanate
"""
from smiles2iupac import smiles_to_iupac


class TestIsocyanate:

    def test_methyl_isocyanate(self):
        assert smiles_to_iupac("CN=C=O") == "methyl isocyanate"

    def test_ethyl_isocyanate(self):
        assert smiles_to_iupac("CCN=C=O") == "ethyl isocyanate"

    def test_propyl_isocyanate(self):
        assert smiles_to_iupac("CCCN=C=O") == "propyl isocyanate"


class TestIsocyanateVsAmide:

    def test_amide_unchanged(self):
        assert smiles_to_iupac("CC(=O)N") == "acetamide"

    def test_isocyanate_not_amide(self):
        result = smiles_to_iupac("CCN=C=O")
        assert "amide" not in result
        assert "isocyanate" in result
