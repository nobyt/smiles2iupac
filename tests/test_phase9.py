"""
Phase 9 テスト: 複数官能基・ポリオール・二酸・ジケトン

対象:
  - ジオール (diol): ethane-1,2-diol, propane-1,2-diol, propane-1,3-diol
  - トリオール (triol): propane-1,2,3-triol
  - ジオン (dione): pentane-2,4-dione, butane-2,3-dione
  - 二酸 (dioic acid): propanedioic acid, butanedioic acid
  - 混合官能基 (highest priority wins): 2-hydroxypropanoic acid, 3-hydroxypropanoic acid
"""
from smiles2iupac import smiles_to_iupac


# ─── ジオール ────────────────────────────────────────────────────────────

class TestDiols:

    def test_ethane_1_2_diol(self):
        assert smiles_to_iupac("OCCO") == "ethane-1,2-diol"

    def test_propane_1_2_diol(self):
        assert smiles_to_iupac("OCC(O)C") == "propane-1,2-diol"

    def test_propane_1_3_diol(self):
        assert smiles_to_iupac("OCCCO") == "propane-1,3-diol"

    def test_butane_1_4_diol(self):
        assert smiles_to_iupac("OCCCCO") == "butane-1,4-diol"

    def test_butane_2_3_diol(self):
        assert smiles_to_iupac("CC(O)C(O)C") == "butane-2,3-diol"


# ─── トリオール ──────────────────────────────────────────────────────────

class TestTriols:

    def test_propane_1_2_3_triol(self):
        # グリセロール
        assert smiles_to_iupac("OCC(O)CO") == "propane-1,2,3-triol"


# ─── ジオン ──────────────────────────────────────────────────────────────

class TestDiones:

    def test_pentane_2_4_dione(self):
        # アセチルアセトン
        assert smiles_to_iupac("CC(=O)CC(=O)C") == "pentane-2,4-dione"

    def test_butane_2_3_dione(self):
        # ジアセチル
        assert smiles_to_iupac("CC(=O)C(=O)C") == "butane-2,3-dione"


# ─── 二酸 ────────────────────────────────────────────────────────────────

class TestDioicAcids:

    def test_propanedioic_acid(self):
        # マロン酸: IUPAC 2013 P-65.1.1.4 保留名
        assert smiles_to_iupac("OC(=O)CC(=O)O") == "malonic acid"

    def test_butanedioic_acid(self):
        # コハク酸: IUPAC 2013 P-65.1.1.4 保留名
        assert smiles_to_iupac("OC(=O)CCC(=O)O") == "succinic acid"

    def test_pentanedioic_acid(self):
        # グルタル酸: IUPAC 2013 P-65.1.1.4 保留名
        assert smiles_to_iupac("OC(=O)CCCC(=O)O") == "glutaric acid"


# ─── 混合官能基 ──────────────────────────────────────────────────────────

class TestMixedFunctionalGroups:

    def test_2_hydroxypropanoic_acid(self):
        # 乳酸: IUPAC 2013 P-65.1.1.4 保留名 "lactic acid"
        assert smiles_to_iupac("CC(O)C(=O)O") == "lactic acid"

    def test_3_hydroxypropanoic_acid(self):
        assert smiles_to_iupac("OCCC(=O)O") == "3-hydroxypropanoic acid"

    def test_2_aminopropanoic_acid(self):
        # アラニン: IUPAC 2013 P-12.1 保留名 "alanine"
        assert smiles_to_iupac("CC(N)C(=O)O") == "alanine"
