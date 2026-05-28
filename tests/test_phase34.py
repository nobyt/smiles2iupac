"""
Phase 34 テスト: エテノール・エチノール (2炭素不飽和アルコール)

対象 (IUPAC P-63.6):
  2 炭素鎖の ene+ol または yne+ol は、ロカントを省略した形で命名。
  "eth-1-en-1-ol" → "ethenol"
  "eth-1-yn-1-ol" → "ethynol"
"""
from smiles2iupac import smiles_to_iupac


class TestEthenol:

    def test_ethenol(self):
        assert smiles_to_iupac("OC=C") == "ethenol"

    def test_ethenol_alt_smiles(self):
        # 同じ分子の別 SMILES
        assert smiles_to_iupac("C=CO") == "ethenol"


class TestEthynol:

    def test_ethynol(self):
        assert smiles_to_iupac("OC#C") == "ethynol"

    def test_ethynol_alt_smiles(self):
        assert smiles_to_iupac("C#CO") == "ethynol"


class TestPropenol:

    def test_prop_2_en_1_ol(self):
        # allyl alcohol: 3炭素鎖なのでロカントあり
        assert smiles_to_iupac("OCC=C") == "prop-2-en-1-ol"

    def test_prop_1_en_1_ol(self):
        # propenol with OH at C1
        assert smiles_to_iupac("OC=CC") == "prop-1-en-1-ol"
