"""
Phase 48 テスト: ジアルデヒド (dial / dialdehyde)

対象 (IUPAC P-66.6.4.1):
  C1 と Cn の両端がアルデヒドの鎖状化合物。
  命名規則: {stem}anedial
  例: O=CCC=O → propanedial
      O=CCCC=O → butanedial
"""
from smiles2iupac import smiles_to_iupac


class TestDial:

    def test_ethanedial(self):
        assert smiles_to_iupac("O=CC=O") == "ethanedial"

    def test_propanedial(self):
        assert smiles_to_iupac("O=CCC=O") == "propanedial"

    def test_butanedial(self):
        assert smiles_to_iupac("O=CCCC=O") == "butanedial"

    def test_pentanedial(self):
        assert smiles_to_iupac("O=CCCCC=O") == "pentanedial"


class TestDialVsAlkanal:

    def test_propanal_unchanged(self):
        assert smiles_to_iupac("CCC=O") == "propanal"

    def test_dial_not_oxyalkanal(self):
        result = smiles_to_iupac("O=CCC=O")
        assert "oxy" not in result
        assert result == "propanedial"
