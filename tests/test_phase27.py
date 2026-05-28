"""
Phase 27 テスト: ベンゼンポリオール・ポリオン

対象 (IUPAC P-63.1.3):
  ベンゼン環上の複数 OH / C=O → benzene-X,Y-diol 等
  命名: benzene-{locs}-{diol|triol|dione}
  ロカント最小化: 最低ロカントセットを選択
"""
from smiles2iupac import smiles_to_iupac


class TestBenzeneDiols:

    def test_benzene_1_2_diol(self):
        # 隣接ジオール
        assert smiles_to_iupac("Oc1ccccc1O") == "benzene-1,2-diol"

    def test_benzene_1_3_diol(self):
        # メタジオール
        assert smiles_to_iupac("Oc1cccc(O)c1") == "benzene-1,3-diol"

    def test_benzene_1_4_diol(self):
        # パラジオール (hydroquinone)
        assert smiles_to_iupac("Oc1ccc(O)cc1") == "benzene-1,4-diol"


class TestBenzeneTriols:

    def test_benzene_1_2_3_triol(self):
        # pyrogallol
        assert smiles_to_iupac("Oc1ccccc1O.Oc1ccccc1O") != "benzene-1,2,3-triol"  # skip multi-mol
        # 単分子テスト
        assert smiles_to_iupac("Oc1cccc(O)c1O") == "benzene-1,2,3-triol"

    def test_benzene_1_2_4_triol(self):
        assert smiles_to_iupac("Oc1ccc(O)c(O)c1") == "benzene-1,2,4-triol"

    def test_benzene_1_3_5_triol(self):
        # phloroglucinol
        assert smiles_to_iupac("Oc1cc(O)cc(O)c1") == "benzene-1,3,5-triol"
