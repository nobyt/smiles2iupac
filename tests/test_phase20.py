"""
Phase 20 テスト: N-置換アミド

対象 (IUPAC P-66.4):
  - 二級アミド (N-alkyl alkanamide)
  - 三級アミド (N,N-dialkyl alkanamide)
"""
from smiles2iupac import smiles_to_iupac


class TestSecondaryAmides:

    def test_n_methylformamide(self):
        # HCONHCH3 → N-methylformamide
        assert smiles_to_iupac("CNC=O") == "N-methylformamide"

    def test_n_methylacetamide(self):
        # CH3CONHCH3 → N-methylacetamide
        assert smiles_to_iupac("CC(=O)NC") == "N-methylacetamide"

    def test_n_ethylacetamide(self):
        # CH3CONHCH2CH3 → N-ethylacetamide
        assert smiles_to_iupac("CC(=O)NCC") == "N-ethylacetamide"

    def test_n_methylpropanamide(self):
        # CH3CH2CONHCH3 → N-methylpropanamide
        assert smiles_to_iupac("CCC(=O)NC") == "N-methylpropanamide"

    def test_n_propylpropanamide(self):
        assert smiles_to_iupac("CCC(=O)NCCC") == "N-propylpropanamide"


class TestTertiaryAmides:

    def test_n_n_dimethylformamide(self):
        # HCON(CH3)2 → N,N-dimethylformamide (DMF)
        assert smiles_to_iupac("CN(C)C=O") == "N,N-dimethylformamide"

    def test_n_n_dimethylacetamide(self):
        # CH3CON(CH3)2 → N,N-dimethylacetamide
        assert smiles_to_iupac("CC(=O)N(C)C") == "N,N-dimethylacetamide"

    def test_n_n_diethylacetamide(self):
        # CH3CON(CH2CH3)2 → N,N-diethylacetamide
        assert smiles_to_iupac("CC(=O)N(CC)CC") == "N,N-diethylacetamide"


class TestPrimaryAmidesRegression:
    """Phase 11 の一級アミドが退行しないことを確認。"""

    def test_formamide(self):
        assert smiles_to_iupac("NC=O") == "formamide"

    def test_acetamide(self):
        assert smiles_to_iupac("CC(=O)N") == "acetamide"
