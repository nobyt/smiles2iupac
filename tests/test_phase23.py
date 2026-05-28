"""
Phase 23 テスト: オキシム (oximes)

対象 (IUPAC P-68.3):
  - ケトオキシム (C=N-OH, C が 2 個の C 隣接): {ketone name} oxime
  - アルドキシム (C=N-OH, C が 0/1 個の C 隣接): {aldehyde name} oxime
"""
from smiles2iupac import smiles_to_iupac


class TestKetoximes:

    def test_propan_2_one_oxime(self):
        # (CH3)2C=NOH → propan-2-one oxime
        assert smiles_to_iupac("CC(=NO)C") == "propan-2-one oxime"

    def test_butan_2_one_oxime(self):
        # CH3C(=NOH)CH2CH3
        assert smiles_to_iupac("CC(=NO)CC") == "butan-2-one oxime"

    def test_pentan_3_one_oxime(self):
        assert smiles_to_iupac("CCC(=NO)CC") == "pentan-3-one oxime"


class TestAldoximes:

    def test_methanal_oxime(self):
        # CH2=NOH → methanal oxime
        assert smiles_to_iupac("C=NO") == "methanal oxime"

    def test_ethanal_oxime(self):
        # CH3CH=NOH → ethanal oxime
        assert smiles_to_iupac("CC=NO") == "ethanal oxime"

    def test_propanal_oxime(self):
        # CH3CH2CH=NOH → propanal oxime
        assert smiles_to_iupac("CCC=NO") == "propanal oxime"
