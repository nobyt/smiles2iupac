"""
Phase 6 テスト: 三環縮合多環芳香族

対象:
  - アントラセン (anthracene)     … 直線状 3 環縮合 (14 C)
  - フェナントレン (phenanthrene) … 角型 3 環縮合 (14 C)
  - 置換アントラセン
    - 9-位 / 10-位 = meso 位 (中央環の非橋頭炭素)
    - 1–8 位 = 外環周辺位

SMILES の原子番号対応 (c1ccc2cc3ccccc3cc2c1, 0-indexed):
  左外環周辺: C13(α,1), C0(β,2), C1(β,3), C2(α,4)
  右外環周辺: C6(α,5), C7(β,6), C8(β,7), C9(α,8)
  meso 位    : C4(9), C11(10)
  橋頭       : C3, C5, C10, C12 (ロカントなし)
"""
import pytest
from smiles2iupac import smiles_to_iupac


# ─── 非置換三環系 ────────────────────────────────────────────────────

class TestTricyclicBase:

    def test_anthracene(self):
        assert smiles_to_iupac("c1ccc2cc3ccccc3cc2c1") == "anthracene"

    def test_phenanthrene(self):
        assert smiles_to_iupac("c1ccc2ccc3ccccc3c2c1") == "phenanthrene"


# ─── 置換アントラセン ────────────────────────────────────────────────

class TestSubstitutedAnthracene:

    def test_9_methylanthracene(self):
        # C4 (meso) にメチル → 9-methylanthracene
        assert smiles_to_iupac("c1ccc2c(C)c3ccccc3cc2c1") == "9-methylanthracene"

    def test_9_chloroanthracene(self):
        # C4 (meso) に塩素 → 9-chloroanthracene
        assert smiles_to_iupac("c1ccc2c(Cl)c3ccccc3cc2c1") == "9-chloroanthracene"

    def test_9_10_dimethylanthracene(self):
        # 両 meso 位 (C4, C11) にメチル → 9,10-dimethylanthracene
        assert smiles_to_iupac("c1ccc2c(C)c3ccccc3c(C)c2c1") == "9,10-dimethylanthracene"

    def test_2_methylanthracene(self):
        # C0 (β 位, locant=2) にメチル → 2-methylanthracene
        assert smiles_to_iupac("Cc1ccc2cc3ccccc3cc2c1") == "2-methylanthracene"

    def test_2_chloroanthracene(self):
        # C0 に塩素 → 2-chloroanthracene
        assert smiles_to_iupac("Clc1ccc2cc3ccccc3cc2c1") == "2-chloroanthracene"
