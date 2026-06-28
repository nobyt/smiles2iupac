"""
Phase 37 テスト: 1,1'-biphenyl (biphenyl は保留名; PIN は 1,1'-biphenyl)

対象 (IUPAC 2013 P-31.1.3):
  2 つのベンゼン環が単結合で繋がる構造の PIN は '1,1'-biphenyl'。
  biphenyl は保留名 (PIN ではない)。
"""
from smiles2iupac import smiles_to_iupac


class TestBiphenylRetainedName:

    def test_biphenyl_kekulized(self):
        assert smiles_to_iupac("c1ccccc1c1ccccc1") == "1,1'-biphenyl"

    def test_biphenyl_dash_notation(self):
        assert smiles_to_iupac("c1ccc(-c2ccccc2)cc1") == "1,1'-biphenyl"

    def test_biphenyl_kekulized_explicit(self):
        assert smiles_to_iupac("C1=CC=CC=C1-C1=CC=CC=C1") == "1,1'-biphenyl"


class TestBiphenylDoesNotBreakBenzene:

    def test_plain_benzene_unchanged(self):
        assert smiles_to_iupac("c1ccccc1") == "benzene"

    def test_chlorobenzene_unchanged(self):
        assert smiles_to_iupac("Clc1ccccc1") == "chlorobenzene"
