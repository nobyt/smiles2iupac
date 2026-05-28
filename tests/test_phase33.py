"""
Phase 33 テスト: シクロアルカンポリオール

対象 (IUPAC P-63.3):
  シクロアルカン環上の複数 OH 基 (diol, triol)。
  命名: cycloalkan-{loc1,loc2}-diol / cyclohexane-{loc1,loc2,loc3}-triol
  ロカント: 全 12 回転で最小ロカントセットを選択。
"""
from smiles2iupac import smiles_to_iupac


class TestCyclohexaneDiols:

    def test_cyclohexane_1_2_diol(self):
        assert smiles_to_iupac("OC1CCCCC1O") == "cyclohexane-1,2-diol"

    def test_cyclohexane_1_3_diol(self):
        assert smiles_to_iupac("OC1CCCC(O)C1") == "cyclohexane-1,3-diol"

    def test_cyclohexane_1_4_diol(self):
        assert smiles_to_iupac("OC1CCC(O)CC1") == "cyclohexane-1,4-diol"

    def test_cyclohexane_1_3_diol_alt_smiles(self):
        # 同じ 1,3-diol の別 SMILES
        assert smiles_to_iupac("OC1CC(O)CCC1") == "cyclohexane-1,3-diol"


class TestCyclohexaneTriols:

    def test_cyclohexane_1_3_5_triol(self):
        assert smiles_to_iupac("OC1CC(O)CC(O)C1") == "cyclohexane-1,3,5-triol"


class TestCyclopentaneDiol:

    def test_cyclopentane_1_2_diol(self):
        assert smiles_to_iupac("OC1CCCC1O") == "cyclopentane-1,2-diol"

    def test_cyclopentane_1_3_diol(self):
        assert smiles_to_iupac("OC1CCC(O)C1") == "cyclopentane-1,3-diol"
