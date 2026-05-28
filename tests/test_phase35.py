"""
Phase 35 テスト: 一置換エテン (monosubstituted ethene)

対象 (IUPAC P-31.1.3.1):
  2 炭素のアルケン (エテン) に 1 つの置換基が付く場合、
  ロカントは自明なので省略する。
  例: "1-chloroethene" → "chloroethene"
"""
from smiles2iupac import smiles_to_iupac


class TestMonosubstitutedEthene:

    def test_chloroethene(self):
        assert smiles_to_iupac("ClC=C") == "chloroethene"

    def test_chloroethene_alt_smiles(self):
        # 同じ分子の別 SMILES
        assert smiles_to_iupac("C=CCl") == "chloroethene"

    def test_bromoethene(self):
        assert smiles_to_iupac("BrC=C") == "bromoethene"

    def test_fluoroethene(self):
        assert smiles_to_iupac("FC=C") == "fluoroethene"

    def test_iodoethene(self):
        assert smiles_to_iupac("IC=C") == "iodoethene"


class TestDisubstitutedEtheneKeepsLocants:

    def test_1_2_dichloroethene(self):
        # 2 置換は引き続きロカントを保持
        assert smiles_to_iupac("ClC=CCl") == "1,2-dichloroethene"
