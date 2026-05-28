"""
Phase 46 テスト: N-置換ヘテロ芳香環 クラッシュ修正

対象: 1H-pyrrole / 1H-imidazole の N 上の置換基により NH が消え、
      シグネチャが一致せず命名クラッシュしていた問題を修正。

命名規則: N-置換 pyrrole / imidazole はロカント 1 に置換基を付与。
  Cn1cccc1   → 1-methylpyrrole
  Cn1ccnc1   → 1-methylimidazole
"""
from smiles2iupac import smiles_to_iupac


class TestNSubstitutedPyrrole:

    def test_1_methylpyrrole(self):
        assert smiles_to_iupac("Cn1cccc1") == "1-methylpyrrole"

    def test_1_ethylpyrrole(self):
        assert smiles_to_iupac("CCn1cccc1") == "1-ethylpyrrole"

    def test_pyrrole_unchanged(self):
        assert smiles_to_iupac("c1cc[nH]c1") == "1H-pyrrole"


class TestNSubstitutedImidazole:

    def test_1_methylimidazole(self):
        assert smiles_to_iupac("Cn1ccnc1") == "1-methylimidazole"

    def test_imidazole_unchanged(self):
        assert smiles_to_iupac("c1cnc[nH]1") == "1H-imidazole"
