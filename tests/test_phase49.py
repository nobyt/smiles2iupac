"""
Phase 49 テスト: 置換基命名の改善

修正内容:
  49-a: アシル置換基 (C=O at root) → acetyl, propanoyl 等
  49-b: シクロアルキル置換基 → cyclohexyl, cyclopentyl 等
  49-c: 尿素 (NC(=O)N) → "urea" 保留名
"""
from smiles2iupac import smiles_to_iupac


class TestAcylSubstituent:

    def test_n_acetylacetamide(self):
        # CC(=O)NC(=O)C: N 上に acetyl 基
        assert smiles_to_iupac("CC(=O)NC(=O)C") == "N-acetylacetamide"

    def test_n_acetyl_label(self):
        result = smiles_to_iupac("CC(=O)NC(=O)C")
        assert "acetyl" in result
        assert "ethyl" not in result


class TestCycloalkylSubstituent:

    def test_1_cyclohexylethan_1_one(self):
        # CC(=O)C1CCCCC1: cyclohexyl ketone
        assert smiles_to_iupac("CC(=O)C1CCCCC1") == "1-cyclohexylethan-1-one"

    def test_cyclopentyl_label(self):
        # CC(=O)C1CCCC1: cyclopentyl ketone
        result = smiles_to_iupac("CC(=O)C1CCCC1")
        assert "cyclopentyl" in result


class TestUrea:

    def test_urea(self):
        assert smiles_to_iupac("NC(=O)N") == "urea"

    def test_urea_not_amide(self):
        result = smiles_to_iupac("NC(=O)N")
        assert result == "urea"
        assert "amide" not in result

    def test_acetamide_unchanged(self):
        # 普通のアミドは変わらない
        assert smiles_to_iupac("CC(=O)N") == "acetamide"
