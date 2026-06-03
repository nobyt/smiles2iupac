"""
Phase 52 テスト: ニトロソ化合物 (nitroso)

対象 (IUPAC P-61.7.1):
  R-N=O パターン。命名: nitroso{alkane}
  例: CN=O → nitrosomethane
"""
from smiles2iupac import smiles_to_iupac


class TestNitroso:

    def test_nitrosomethane(self):
        assert smiles_to_iupac("CN=O") == "nitrosomethane"

    def test_nitrosoethane(self):
        assert smiles_to_iupac("CCN=O") == "nitrosoethane"

    def test_nitrosopropane(self):
        assert smiles_to_iupac("CCCN=O") == "1-nitrosopropane"


class TestNitrosoVsNitro:

    def test_nitro_not_nitroso(self):
        # ニトロ化合物はニトロソではない
        result = smiles_to_iupac("CC[N+](=O)[O-]")
        assert "nitroso" not in result

    def test_nitrobenzene_unchanged(self):
        result = smiles_to_iupac("c1ccc([N+](=O)[O-])cc1")
        assert result == "nitrobenzene"
