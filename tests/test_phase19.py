"""
Phase 19 テスト: イミン (C=N-H)

対象 (IUPAC P-66.1.5):
  - 第一級イミン (=NH): alkan-N-imine
"""
from smiles2iupac import smiles_to_iupac


class TestImines:

    def test_methan_1_imine(self):
        # CH2=NH → methan-1-imine
        assert smiles_to_iupac("C=N") == "methan-1-imine"

    def test_ethan_1_imine(self):
        # CH3-CH=NH → ethan-1-imine (N は C1 に付く)
        assert smiles_to_iupac("CC=N") == "ethan-1-imine"

    def test_propan_2_imine(self):
        # CH3-C(=NH)-CH3 → propan-2-imine
        assert smiles_to_iupac("CC(=N)C") == "propan-2-imine"

    def test_pentan_3_imine(self):
        # CCC(=NH)CC → pentan-3-imine
        assert smiles_to_iupac("CCC(=N)CC") == "pentan-3-imine"

    def test_propan_1_imine(self):
        # CH3-CH2-CH=NH → propan-1-imine (N is at C1)
        assert smiles_to_iupac("CCC=N") == "propan-1-imine"

    def test_butan_2_imine(self):
        # CH3-C(=NH)-CH2-CH3 → butan-2-imine
        assert smiles_to_iupac("CC(=N)CC") == "butan-2-imine"
