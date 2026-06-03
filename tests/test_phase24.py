"""
Phase 24 テスト: 酸無水物 (acid anhydrides)

対象 (IUPAC P-65.1.2.2, P-65.1.1.3.2):
  - 対称無水物: acetic anhydride (保留名使用), formic anhydride
  - 非対称無水物: acetic butanoic anhydride (アルファベット順)
"""
from smiles2iupac import smiles_to_iupac


class TestSymmetricAnhydrides:

    def test_acetic_anhydride(self):
        # (CH3CO)2O → acetic anhydride (IUPAC 2013 P-65.1.1.3.2 保留優先名)
        assert smiles_to_iupac("CC(=O)OC(=O)C") == "acetic anhydride"

    def test_formic_anhydride(self):
        # (HCOO)2 → formic anhydride (IUPAC 2013 formic acid 保留名)
        assert smiles_to_iupac("O=COC=O") == "formic anhydride"

    def test_propanoic_anhydride(self):
        assert smiles_to_iupac("CCC(=O)OC(=O)CC") == "propanoic anhydride"


class TestAsymmetricAnhydrides:

    def test_acetic_butanoic_anhydride(self):
        # CH3CO-O-CO-C3H7 → acetic butanoic anhydride (a<b alphabetical)
        assert smiles_to_iupac("CC(=O)OC(=O)CCC") == "acetic butanoic anhydride"

    def test_acetic_formic_anhydride(self):
        # HCOO-CO-CH3 → acetic formic anhydride (a<f alphabetical)
        assert smiles_to_iupac("C(=O)OC(=O)C") == "acetic formic anhydride"

    def test_acetic_propanoic_anhydride(self):
        assert smiles_to_iupac("CC(=O)OC(=O)CC") == "acetic propanoic anhydride"
