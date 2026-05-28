"""
Phase 15 テスト: 架橋多環系・スピロ化合物

対象:
  (a) スピロ化合物: spiro[m.n]alkane
  (b) 架橋二環 (von Baeyer): bicyclo[l.m.n]alkane
"""
from smiles2iupac import smiles_to_iupac


# ─── スピロ化合物 ───────────────────────────────────────────────────────

class TestSpiroCompounds:

    def test_spiro_4_5_decane(self):
        # 5員環 + 6員環をスピロ原子で連結 → spiro[4.5]decane
        assert smiles_to_iupac("C1CCC2(CC1)CCCC2") == "spiro[4.5]decane"

    def test_spiro_2_2_pentane(self):
        # 3員環 + 3員環 → spiro[2.2]pentane
        assert smiles_to_iupac("C1CC12CC2") == "spiro[2.2]pentane"

    def test_spiro_5_5_undecane(self):
        # 6員環 + 6員環 → spiro[5.5]undecane
        assert smiles_to_iupac("C1CCCCC12CCCCC2") == "spiro[5.5]undecane"

    def test_spiro_4_4_nonane(self):
        # 5員環 + 5員環 → spiro[4.4]nonane
        assert smiles_to_iupac("C1CCCC12CCCC2") == "spiro[4.4]nonane"


# ─── 架橋二環 (von Baeyer) ─────────────────────────────────────────────

class TestBridgedBicyclics:

    def test_bicyclo_2_2_1_heptane(self):
        # ノルボルナン: bicyclo[2.2.1]heptane
        assert smiles_to_iupac("C1CC2CCC1C2") == "bicyclo[2.2.1]heptane"

    def test_bicyclo_2_2_2_octane(self):
        # bicyclo[2.2.2]octane
        assert smiles_to_iupac("C1CC2CCC1CC2") == "bicyclo[2.2.2]octane"

    def test_bicyclo_3_1_1_heptane(self):
        # bicyclo[3.1.1]heptane
        assert smiles_to_iupac("C12CCCC(C1)C2") == "bicyclo[3.1.1]heptane"

    def test_bicyclo_2_1_0_pentane(self):
        # bicyclo[2.1.0]pentane
        assert smiles_to_iupac("C12CC1CC2") == "bicyclo[2.1.0]pentane"
