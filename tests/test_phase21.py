"""
Phase 21 テスト: スルホキシド・スルホン・スルホンアミド

対象 (IUPAC P-63.6):
  - スルホキシド (S=O): (methylsulfinyl)methane 形式
  - スルホン (S(=O)₂): (methylsulfonyl)methane 形式
  - スルホンアミド (S(=O)₂-N): alkanesulfonamide 形式
"""
from smiles2iupac import smiles_to_iupac


class TestSulfoxides:

    def test_dimethyl_sulfoxide(self):
        # DMSO: (CH3)2S=O → (methylsulfinyl)methane
        assert smiles_to_iupac("CS(=O)C") == "(methylsulfinyl)methane"

    def test_diethyl_sulfoxide(self):
        assert smiles_to_iupac("CCS(=O)CC") == "(ethylsulfinyl)ethane"

    def test_ethyl_methyl_sulfoxide(self):
        assert smiles_to_iupac("CCS(=O)C") == "(methylsulfinyl)ethane"


class TestSulfones:

    def test_dimethyl_sulfone(self):
        # (CH3)2SO2 → (methylsulfonyl)methane
        assert smiles_to_iupac("CS(=O)(=O)C") == "(methylsulfonyl)methane"

    def test_diethyl_sulfone(self):
        assert smiles_to_iupac("CCS(=O)(=O)CC") == "(ethylsulfonyl)ethane"

    def test_ethyl_methyl_sulfone(self):
        assert smiles_to_iupac("CCS(=O)(=O)C") == "(methylsulfonyl)ethane"


class TestSulfonamides:

    def test_methanesulfonamide(self):
        # CH3SO2NH2 → methanesulfonamide
        assert smiles_to_iupac("CS(=O)(=O)N") == "methanesulfonamide"

    def test_ethanesulfonamide(self):
        assert smiles_to_iupac("CCS(=O)(=O)N") == "ethanesulfonamide"

    def test_n_methylmethanesulfonamide(self):
        # CH3SO2NHCH3 → N-methylmethanesulfonamide
        assert smiles_to_iupac("CS(=O)(=O)NC") == "N-methylmethanesulfonamide"
