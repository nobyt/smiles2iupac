"""
Phase 5 テスト: ナフタレン誘導体

対象:
  - ナフタレン (naphthalene)
  - 一置換ナフタレン (1-/2- 位の区別)
  - 多置換ナフタレン
"""
import pytest
from smiles2iupac import smiles_to_iupac


# ─── 非置換 ──────────────────────────────────────────────────────────

class TestNaphthalene:

    def test_naphthalene(self):
        assert smiles_to_iupac("c1ccc2ccccc2c1") == "naphthalene"


# ─── 一置換ナフタレン ────────────────────────────────────────────────

class TestMonosubstitutedNaphthalene:

    def test_1_chloronaphthalene(self):
        # alpha 位 (1-位)
        assert smiles_to_iupac("Clc1cccc2ccccc12") == "1-chloronaphthalene"

    def test_2_chloronaphthalene(self):
        # beta 位 (2-位)
        assert smiles_to_iupac("Clc1ccc2ccccc2c1") == "2-chloronaphthalene"

    def test_1_methylnaphthalene(self):
        assert smiles_to_iupac("Cc1cccc2ccccc12") == "1-methylnaphthalene"

    def test_2_methylnaphthalene(self):
        assert smiles_to_iupac("Cc1ccc2ccccc2c1") == "2-methylnaphthalene"

    def test_2_bromonaphthalene(self):
        assert smiles_to_iupac("Brc1ccc2ccccc2c1") == "2-bromonaphthalene"


# ─── 多置換ナフタレン ────────────────────────────────────────────────

class TestPolysubstitutedNaphthalene:

    def test_1_4_dimethylnaphthalene(self):
        assert smiles_to_iupac("Cc1ccc(C)c2ccccc12") == "1,4-dimethylnaphthalene"

    def test_2_6_dimethylnaphthalene(self):
        assert smiles_to_iupac("Cc1ccc2cc(C)ccc2c1") == "2,6-dimethylnaphthalene"
