"""
Phase 37 テスト: ビフェニル (biphenyl)

対象 (IUPAC P-31.1.6.2):
  2 つのベンゼン環が単結合で繋がる構造に保留名 "biphenyl" を使用。
  一方のリングを親名、他方を substituent として収集し、
  assemble_ring_name で "biphenyl" にまとめる。
"""
from smiles2iupac import smiles_to_iupac


class TestBiphenylRetainedName:

    def test_biphenyl_kekulized(self):
        assert smiles_to_iupac("c1ccccc1c1ccccc1") == "biphenyl"

    def test_biphenyl_dash_notation(self):
        assert smiles_to_iupac("c1ccc(-c2ccccc2)cc1") == "biphenyl"

    def test_biphenyl_kekulized_explicit(self):
        assert smiles_to_iupac("C1=CC=CC=C1-C1=CC=CC=C1") == "biphenyl"


class TestBiphenylDoesNotBreakBenzene:

    def test_plain_benzene_unchanged(self):
        assert smiles_to_iupac("c1ccccc1") == "benzene"

    def test_chlorobenzene_unchanged(self):
        assert smiles_to_iupac("Clc1ccccc1") == "chlorobenzene"
