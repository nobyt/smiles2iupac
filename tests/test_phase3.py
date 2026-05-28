"""
Phase 3 テスト: シクロアルカン・芳香族 (ベンゼン誘導体)
"""
import pytest
from smiles2iupac import smiles_to_iupac


# ─── 非置換シクロアルカン ────────────────────────────────────────────

class TestSimpleCycloalkanes:

    def test_cyclopropane(self):
        assert smiles_to_iupac("C1CC1") == "cyclopropane"

    def test_cyclobutane(self):
        assert smiles_to_iupac("C1CCC1") == "cyclobutane"

    def test_cyclopentane(self):
        assert smiles_to_iupac("C1CCCC1") == "cyclopentane"

    def test_cyclohexane(self):
        assert smiles_to_iupac("C1CCCCC1") == "cyclohexane"

    def test_cycloheptane(self):
        assert smiles_to_iupac("C1CCCCCC1") == "cycloheptane"


# ─── 一置換シクロアルカン ────────────────────────────────────────────

class TestMonosubstitutedCycloalkanes:

    def test_methylcyclopropane(self):
        assert smiles_to_iupac("C1CC1C") == "methylcyclopropane"

    def test_methylcyclopentane(self):
        assert smiles_to_iupac("C1CCCC1C") == "methylcyclopentane"

    def test_methylcyclohexane(self):
        assert smiles_to_iupac("C1CCCCC1C") == "methylcyclohexane"

    def test_ethylcyclohexane(self):
        assert smiles_to_iupac("C1CCCCC1CC") == "ethylcyclohexane"

    def test_chlorocyclohexane(self):
        assert smiles_to_iupac("C1CCCCC1Cl") == "chlorocyclohexane"


# ─── 多置換シクロアルカン ────────────────────────────────────────────

class TestPolysubstitutedCycloalkanes:

    def test_1_2_dimethylcyclopentane(self):
        assert smiles_to_iupac("CC1CCCC1C") == "1,2-dimethylcyclopentane"

    def test_1_3_dimethylcyclohexane(self):
        assert smiles_to_iupac("CC1CCCC(C)C1") == "1,3-dimethylcyclohexane"

    def test_1_4_dimethylcyclohexane(self):
        assert smiles_to_iupac("CC1CCC(C)CC1") == "1,4-dimethylcyclohexane"

    def test_1_bromo_3_methylcyclopentane(self):
        assert smiles_to_iupac("BrC1CC(C)CC1") == "1-bromo-3-methylcyclopentane"


# ─── シクロアルカンの官能基 ─────────────────────────────────────────

class TestCycloalkaneWithFunctionalGroups:

    def test_cyclohexanol(self):
        assert smiles_to_iupac("OC1CCCCC1") == "cyclohexan-1-ol"

    def test_cyclohexanone(self):
        assert smiles_to_iupac("O=C1CCCCC1") == "cyclohexan-1-one"

    def test_cyclopentanol(self):
        assert smiles_to_iupac("OC1CCCC1") == "cyclopentan-1-ol"


# ─── ベンゼン ────────────────────────────────────────────────────────

class TestBenzene:

    def test_benzene(self):
        assert smiles_to_iupac("c1ccccc1") == "benzene"


# ─── 一置換ベンゼン ─────────────────────────────────────────────────

class TestMonosubstitutedBenzene:

    def test_methylbenzene(self):
        # toluene の IUPAC 系統名
        assert smiles_to_iupac("Cc1ccccc1") == "methylbenzene"

    def test_chlorobenzene(self):
        assert smiles_to_iupac("Clc1ccccc1") == "chlorobenzene"

    def test_bromobenzene(self):
        assert smiles_to_iupac("Brc1ccccc1") == "bromobenzene"

    def test_ethylbenzene(self):
        assert smiles_to_iupac("CCc1ccccc1") == "ethylbenzene"


# ─── 多置換ベンゼン ─────────────────────────────────────────────────

class TestPolysubstitutedBenzene:

    def test_1_2_dimethylbenzene(self):
        # o-xylene
        assert smiles_to_iupac("Cc1ccccc1C") == "1,2-dimethylbenzene"

    def test_1_3_dimethylbenzene(self):
        # m-xylene
        assert smiles_to_iupac("Cc1cccc(C)c1") == "1,3-dimethylbenzene"

    def test_1_4_dimethylbenzene(self):
        # p-xylene
        assert smiles_to_iupac("Cc1ccc(C)cc1") == "1,4-dimethylbenzene"

    def test_1_bromo_4_chlorobenzene(self):
        assert smiles_to_iupac("Brc1ccc(Cl)cc1") == "1-bromo-4-chlorobenzene"

    def test_1_chloro_2_methylbenzene(self):
        assert smiles_to_iupac("Clc1ccccc1C") == "1-chloro-2-methylbenzene"
