"""
Phase 53 テスト: azide substitutive PIN (IUPAC 2013 P-65.3.1)

R-N3 → azido{alkane} (PIN; "{alkyl} azide" is retained acceptable).
  例: CCN=[N+]=[N-] → azidoethane
"""
from smiles2iupac import smiles_to_iupac


class TestAzide:

    def test_methyl_azide(self):
        assert smiles_to_iupac("CN=[N+]=[N-]") == "azidomethane"

    def test_ethyl_azide(self):
        assert smiles_to_iupac("CCN=[N+]=[N-]") == "azidoethane"

    def test_propyl_azide(self):
        assert smiles_to_iupac("CCCN=[N+]=[N-]") == "azidopropane"


class TestAzideVsAmine:

    def test_amine_unchanged(self):
        assert smiles_to_iupac("CCN") == "ethanamine"

    def test_azide_not_amine(self):
        result = smiles_to_iupac("CCN=[N+]=[N-]")
        assert "amine" not in result
        assert "azido" in result
