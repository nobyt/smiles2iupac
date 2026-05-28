"""
Phase 16 テスト: 4 環以上の縮合多環芳香族

対象:
  - pyrene (4 環, C16H10)
  - chrysene (4 環, C18H12)
  - triphenylene (4 環, C18H12)
  - perylene (5 環, C20H12)
"""
from smiles2iupac import smiles_to_iupac


class TestFourRingPolycyclics:

    def test_pyrene(self):
        # C16H10: 4 環縮合, 非線形配置
        assert smiles_to_iupac("c1cc2ccc3cccc4ccc(c1)c2c34") == "pyrene"

    def test_triphenylene(self):
        # C18H12: スター型 4 環 (星形環グラフ)
        assert smiles_to_iupac("c1ccc2c(c1)c1ccccc1c1ccccc21") == "triphenylene"

    def test_chrysene(self):
        # C18H12: 線形 4 環 (パス型環グラフ)
        assert smiles_to_iupac("c1ccc2cc3ccc4ccccc4c3cc2c1") == "chrysene"


class TestFiveRingPolycyclics:

    def test_perylene(self):
        # C20H12: 2 ナフタレンのペリ縮合
        assert smiles_to_iupac("c1cc2cccc3c4cccc5cccc(c(c1)c23)c54") == "perylene"
