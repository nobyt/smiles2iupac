"""
Phase 26 テスト: スルホン酸・スルフィド

対象 (IUPAC P-63.6, P-63.5):
  スルホン酸: C-S(=O)₂-OH → {stem}anesulfonic acid
  スルフィド: C-S-C → {alkyl1} {alkyl2} sulfide (functional class naming)
"""
from smiles2iupac import smiles_to_iupac


class TestSulfonicAcids:

    def test_methanesulfonic_acid(self):
        assert smiles_to_iupac("CS(=O)(=O)O") == "methanesulfonic acid"

    def test_ethanesulfonic_acid(self):
        assert smiles_to_iupac("CCS(=O)(=O)O") == "ethanesulfonic acid"

    def test_propanesulfonic_acid(self):
        assert smiles_to_iupac("CCCS(=O)(=O)O") == "propane-1-sulfonic acid"

    def test_butanesulfonic_acid(self):
        assert smiles_to_iupac("CCCCS(=O)(=O)O") == "butane-1-sulfonic acid"


class TestSulfides:

    def test_dimethyl_sulfide(self):
        assert smiles_to_iupac("CSC") == "dimethyl sulfide"

    def test_ethyl_methyl_sulfide(self):
        assert smiles_to_iupac("CSCC") == "ethyl methyl sulfide"

    def test_diethyl_sulfide(self):
        assert smiles_to_iupac("CCSCC") == "diethyl sulfide"

    def test_methyl_phenyl_sulfide(self):
        assert smiles_to_iupac("CSc1ccccc1") == "methyl phenyl sulfide"
