"""
Phase 18 テスト: 二級・三級アミン

対象 (IUPAC P-62.2):
  - 二級アミン (N-alkyl-alkanamine)
  - 三級アミン (N,N-dialkyl-alkanamine)
"""
from smiles2iupac import smiles_to_iupac


class TestSecondaryAmines:

    def test_n_methylmethanamine(self):
        # (CH3)2NH → N-methylmethanamine
        assert smiles_to_iupac("CNC") == "N-methylmethanamine"

    def test_n_methylethanamine(self):
        # CH3NHCH2CH3 → N-methylethanamine (親: ethane > methane)
        assert smiles_to_iupac("CCNC") == "N-methylethanamine"

    def test_n_ethylethanamine(self):
        # (CH3CH2)2NH → N-ethylethanamine
        assert smiles_to_iupac("CCNCC") == "N-ethylethanamine"

    def test_n_propylpropanamine(self):
        # (C3H7)2NH → N-propylpropan-1-amine (IUPAC 2013 requires locant for 3+ C)
        assert smiles_to_iupac("CCCNCCC") == "N-propylpropan-1-amine"


class TestTertiaryAmines:

    def test_n_n_dimethylmethanamine(self):
        # (CH3)3N → N,N-dimethylmethanamine
        assert smiles_to_iupac("CN(C)C") == "N,N-dimethylmethanamine"

    def test_n_n_diethylethanamine(self):
        # (CH3CH2)3N → N,N-diethylethanamine
        assert smiles_to_iupac("CCN(CC)CC") == "N,N-diethylethanamine"

    def test_n_n_dimethylethanamine(self):
        # (CH3)2NCH2CH3 → N,N-dimethylethanamine (親: ethane > methane)
        assert smiles_to_iupac("CCN(C)C") == "N,N-dimethylethanamine"
