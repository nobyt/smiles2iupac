"""
Phase 40 テスト: 環上のエテニル・エチニル置換基 (vinyl/ethynyl on ring)

対象 (IUPAC P-31.1.6):
  環 (benzene / cycloalkane) に C=C または C≡C の 2 炭素鎖が付く場合、
  環が主構造となり、側鎖は "ethenyl" / "ethynyl" 置換基として命名する。
  規則: 環の炭素数 > 非環炭素数 のとき ring path を選択。
"""
from smiles2iupac import smiles_to_iupac


class TestEthenylBenzene:

    def test_ethenylbenzene_smiles1(self):
        # IUPAC 2013 P-31.1.3.4: styrene が保留優先名
        assert smiles_to_iupac("C=Cc1ccccc1") == "styrene"

    def test_ethenylbenzene_smiles2(self):
        # 別 SMILES 表記 (same compound, same retained name)
        assert smiles_to_iupac("c1ccc(C=C)cc1") == "styrene"


class TestEthynylBenzene:

    def test_ethynylbenzene(self):
        assert smiles_to_iupac("C#Cc1ccccc1") == "ethynylbenzene"

    def test_ethynylbenzene_reversed(self):
        assert smiles_to_iupac("c1ccc(C#C)cc1") == "ethynylbenzene"


class TestEthenylCycloalkane:

    def test_ethenylcyclohexane(self):
        assert smiles_to_iupac("C=CC1CCCCC1") == "ethenylcyclohexane"

    def test_ethenylcyclopentane(self):
        assert smiles_to_iupac("C=CC1CCCC1") == "ethenylcyclopentane"


class TestAcyclicAlkeneUnchanged:

    def test_propene_unchanged(self):
        # 非環式のアルケンは引き続き acyclic path
        assert smiles_to_iupac("C=CC") == "prop-1-ene"

    def test_but1ene_unchanged(self):
        assert smiles_to_iupac("C=CCC") == "but-1-ene"
