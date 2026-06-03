"""
Phase 14 テスト: シクロアルケン / シクロジエン

対象:
  - 単純シクロアルケン (cyclopent-1-ene, cyclohex-1-ene, ...)
  - シクロジエン (cyclohexa-1,3-diene)
  - 置換シクロアルケン (1-methylcyclohex-1-ene 等)
  - 官能基付きシクロアルケン (cyclohex-2-en-1-ol 等)
"""
from smiles2iupac import smiles_to_iupac


# ─── 単純シクロアルケン ─────────────────────────────────────────────────

class TestSimpleCycloalkenes:

    def test_cyclopent_1_ene(self):
        # 5 員環に 1 二重結合
        assert smiles_to_iupac("C1=CCCC1") == "cyclopentene"

    def test_cyclohex_1_ene(self):
        # 6 員環に 1 二重結合
        assert smiles_to_iupac("C1=CCCCC1") == "cyclohexene"

    def test_cyclohept_1_ene(self):
        # 7 員環
        assert smiles_to_iupac("C1=CCCCCC1") == "cycloheptene"

    def test_cyclobut_1_ene(self):
        # 4 員環
        assert smiles_to_iupac("C1=CCC1") == "cyclobutene"


# ─── シクロジエン ──────────────────────────────────────────────────────

class TestCyclodienes:

    def test_cyclohexa_1_3_diene(self):
        # 共役ジエン
        assert smiles_to_iupac("C1=CC=CCC1") == "cyclohexa-1,3-diene"

    def test_cyclohexa_1_4_diene(self):
        # 非共役ジエン
        assert smiles_to_iupac("C1=CCC=CC1") == "cyclohexa-1,4-diene"


# ─── 置換シクロアルケン ────────────────────────────────────────────────

class TestSubstitutedCycloalkenes:

    def test_1_methylcyclohex_1_ene(self):
        # メチルが二重結合の C1 に付く
        assert smiles_to_iupac("CC1=CCCCC1") == "1-methylcyclohex-1-ene"

    def test_3_methylcyclohex_1_ene(self):
        # メチルが C3 に付く
        assert smiles_to_iupac("CC1CCCC=C1") == "3-methylcyclohex-1-ene"


# ─── 官能基付きシクロアルケン ──────────────────────────────────────────

class TestCycloalkeneWithFunctionalGroups:

    def test_cyclohex_2_en_1_ol(self):
        # allylic alcohol: OH at C1, double bond at C2=C3
        assert smiles_to_iupac("OC1C=CCCC1") == "cyclohex-2-en-1-ol"
