"""
Phase 32 テスト: フェニルアルカン酸・フェニルアルカノール

対象 (IUPAC P-65.1.1):
  ベンゼン環を置換基として持つ直鎖酸・アルコール・アルデヒド。
  命名: {n}-phenyl{alkanoic acid/alkanol/alkanal}
"""
from smiles2iupac import smiles_to_iupac


class TestPhenylalkanoicAcids:

    def test_2_phenylethanoic_acid(self):
        # PhCH2COOH → 2-phenylacetic acid (acetic acid 保留名ベース)
        assert smiles_to_iupac("c1ccc(CC(=O)O)cc1") == "2-phenylacetic acid"

    def test_3_phenylpropanoic_acid(self):
        # PhCH2CH2COOH
        assert smiles_to_iupac("c1ccc(CCC(=O)O)cc1") == "3-phenylpropanoic acid"

    def test_4_phenylbutanoic_acid(self):
        # PhCH2CH2CH2COOH
        assert smiles_to_iupac("c1ccc(CCCC(=O)O)cc1") == "4-phenylbutanoic acid"


class TestPhenylalkanols:

    def test_phenylmethanol(self):
        # PhCH2OH (benzyl alcohol)
        assert smiles_to_iupac("c1ccc(CO)cc1") == "phenylmethanol"

    def test_2_phenylethan_1_ol(self):
        # PhCH2CH2OH → 2-phenylethanol (ethanol 保留名ベース)
        assert smiles_to_iupac("c1ccc(CCO)cc1") == "2-phenylethanol"


class TestPhenylaldehydes:

    def test_2_phenylethanal(self):
        # PhCH2CHO (phenylacetaldehyde)
        assert smiles_to_iupac("c1ccc(CC=O)cc1") == "2-phenylethanal"
