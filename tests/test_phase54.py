"""
Phase 54 テスト: isocyanate substitutive PIN (IUPAC 2013 P-65.3.1)

R-N=C=O → isocyanato{alkane} (PIN; "{alkyl} isocyanate" is retained acceptable).
  例: CCN=C=O → isocyanatoethane
"""
from smiles2iupac import smiles_to_iupac


class TestIsocyanate:

    def test_methyl_isocyanate(self):
        assert smiles_to_iupac("CN=C=O") == "isocyanatomethane"

    def test_ethyl_isocyanate(self):
        assert smiles_to_iupac("CCN=C=O") == "isocyanatoethane"

    def test_propyl_isocyanate(self):
        assert smiles_to_iupac("CCCN=C=O") == "isocyanatopropane"


class TestIsocyanateVsAmide:

    def test_amide_unchanged(self):
        assert smiles_to_iupac("CC(=O)N") == "acetamide"

    def test_isocyanate_not_amide(self):
        result = smiles_to_iupac("CCN=C=O")
        assert "amide" not in result
        assert "isocyanato" in result
