"""
Phase 22 テスト: ヘテロ芳香環上の外環官能基

対象 (IUPAC P-31.1.6):
  - カルボン酸: pyridine-3-carboxylic acid (nicotinic acid) 形式
  - アルデヒド: pyridine-3-carbaldehyde 形式
  - アミド: furan-2-carboxamide 形式
  - ニトリル: furan-2-carbonitrile 形式
"""
from smiles2iupac import smiles_to_iupac


class TestPyridineCarboxylicAcid:

    def test_nicotinic_acid_smiles_correct_direction(self):
        # COOH が N から 2 ステップ（直接表記）
        assert smiles_to_iupac("O=C(O)c1cccnc1") == "pyridine-3-carboxylic acid"

    def test_nicotinic_acid_smiles_reversed_direction(self):
        # COOH が逆方向 SMILES（ロカント修正が必要）
        assert smiles_to_iupac("c1ccc(cn1)C(=O)O") == "pyridine-3-carboxylic acid"

    def test_pyridine_4_carboxylic_acid(self):
        # COOH が N と para 位置（4 位）
        assert smiles_to_iupac("OC(=O)c1ccncc1") == "pyridine-4-carboxylic acid"


class TestFuranDerivatives:

    def test_furan_2_carboxylic_acid(self):
        assert smiles_to_iupac("c1ccoc1C(=O)O") == "furan-2-carboxylic acid"

    def test_furan_2_carbaldehyde(self):
        assert smiles_to_iupac("O=Cc1ccco1") == "furan-2-carbaldehyde"

    def test_furan_2_carbonitrile(self):
        assert smiles_to_iupac("c1ccoc1C#N") == "furan-2-carbonitrile"

    def test_furan_2_carboxamide(self):
        assert smiles_to_iupac("c1ccoc1C(N)=O") == "furan-2-carboxamide"


class TestPyridineOtherFG:

    def test_pyridine_3_carbaldehyde(self):
        assert smiles_to_iupac("c1ccc(cn1)C=O") == "pyridine-3-carbaldehyde"

    def test_pyridine_3_carbonitrile(self):
        assert smiles_to_iupac("c1ccc(cn1)C#N") == "pyridine-3-carbonitrile"


class TestPyrroleDerivative:

    def test_pyrrole_2_carboxylic_acid(self):
        assert smiles_to_iupac("c1cc[nH]c1C(=O)O") == "1H-pyrrole-2-carboxylic acid"
