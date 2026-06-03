"""
Phase 21 テスト: スルホキシド・スルホン・スルホンアミド

対象 (IUPAC P-63.6):
  - スルホキシド (S=O): dimethyl sulfoxide 形式
  - スルホン (S(=O)₂): dimethyl sulfone 形式
  - スルホンアミド (S(=O)₂-N): alkanesulfonamide 形式
"""
from smiles2iupac import smiles_to_iupac


class TestSulfoxides:

    def test_dimethyl_sulfoxide(self):
        # DMSO: (CH3)2S=O → dimethyl sulfoxide
        assert smiles_to_iupac("CS(=O)C") == "dimethyl sulfoxide"

    def test_diethyl_sulfoxide(self):
        assert smiles_to_iupac("CCS(=O)CC") == "diethyl sulfoxide"

    def test_ethyl_methyl_sulfoxide(self):
        assert smiles_to_iupac("CCS(=O)C") == "ethyl methyl sulfoxide"


class TestSulfones:

    def test_dimethyl_sulfone(self):
        # (CH3)2SO2 → dimethyl sulfone
        assert smiles_to_iupac("CS(=O)(=O)C") == "dimethyl sulfone"

    def test_diethyl_sulfone(self):
        assert smiles_to_iupac("CCS(=O)(=O)CC") == "diethyl sulfone"

    def test_ethyl_methyl_sulfone(self):
        assert smiles_to_iupac("CCS(=O)(=O)C") == "ethyl methyl sulfone"


class TestSulfonamides:

    def test_methanesulfonamide(self):
        # CH3SO2NH2 → methanesulfonamide
        assert smiles_to_iupac("CS(=O)(=O)N") == "methanesulfonamide"

    def test_ethanesulfonamide(self):
        assert smiles_to_iupac("CCS(=O)(=O)N") == "ethanesulfonamide"

    def test_n_methylmethanesulfonamide(self):
        # CH3SO2NHCH3 → N-methylmethanesulfonamide
        assert smiles_to_iupac("CS(=O)(=O)NC") == "N-methylmethanesulfonamide"
