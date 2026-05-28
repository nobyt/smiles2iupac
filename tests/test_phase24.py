"""
Phase 24 テスト: 酸無水物 (acid anhydrides)

対象 (IUPAC P-65.1.2.2):
  - 対称無水物: ethanoic anhydride
  - 非対称無水物: butanoic ethanoic anhydride (アルファベット順)
"""
from smiles2iupac import smiles_to_iupac


class TestSymmetricAnhydrides:

    def test_ethanoic_anhydride(self):
        # (CH3CO)2O → ethanoic anhydride
        assert smiles_to_iupac("CC(=O)OC(=O)C") == "ethanoic anhydride"

    def test_methanoic_anhydride(self):
        # (HCOO)2 → methanoic anhydride
        assert smiles_to_iupac("O=COC=O") == "methanoic anhydride"

    def test_propanoic_anhydride(self):
        assert smiles_to_iupac("CCC(=O)OC(=O)CC") == "propanoic anhydride"


class TestAsymmetricAnhydrides:

    def test_butanoic_ethanoic_anhydride(self):
        # CH3CO-O-CO-C3H7 → butanoic ethanoic anhydride
        assert smiles_to_iupac("CC(=O)OC(=O)CCC") == "butanoic ethanoic anhydride"

    def test_ethanoic_methanoic_anhydride(self):
        # HCOO-CO-CH3
        assert smiles_to_iupac("C(=O)OC(=O)C") == "ethanoic methanoic anhydride"

    def test_ethanoic_propanoic_anhydride(self):
        assert smiles_to_iupac("CC(=O)OC(=O)CC") == "ethanoic propanoic anhydride"
