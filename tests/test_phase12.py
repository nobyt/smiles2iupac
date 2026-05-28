"""
Phase 12 テスト: ヘテロ環保留名

対象:
  - 芳香族ヘテロ環: pyridine, furan, thiophene, pyrrole, imidazole, pyrimidine 等
  - 飽和ヘテロ環: piperidine, oxane, morpholine 等
  - 縮合ヘテロ環: quinoline, indole 等
  - 置換ヘテロ環: 2-methylpyridine, 3-chloropyridine 等
"""
from smiles2iupac import smiles_to_iupac


# ─── 6 員芳香族ヘテロ環 ─────────────────────────────────────────────────

class TestSixMemberAromaticHeterocycles:

    def test_pyridine(self):
        assert smiles_to_iupac("c1ccncc1") == "pyridine"

    def test_pyrimidine(self):
        assert smiles_to_iupac("c1ccnc(n1)") == "pyrimidine"

    def test_pyrazine(self):
        assert smiles_to_iupac("c1cnccn1") == "pyrazine"

    def test_pyridazine(self):
        assert smiles_to_iupac("c1ccnnc1") == "pyridazine"


# ─── 5 員芳香族ヘテロ環 ─────────────────────────────────────────────────

class TestFiveMemberAromaticHeterocycles:

    def test_furan(self):
        assert smiles_to_iupac("c1ccoc1") == "furan"

    def test_thiophene(self):
        assert smiles_to_iupac("c1ccsc1") == "thiophene"

    def test_pyrrole(self):
        assert smiles_to_iupac("c1cc[nH]c1") == "1H-pyrrole"

    def test_imidazole(self):
        assert smiles_to_iupac("c1c[nH]cn1") == "1H-imidazole"


# ─── 飽和ヘテロ環 ───────────────────────────────────────────────────────

class TestSaturatedHeterocycles:

    def test_piperidine(self):
        assert smiles_to_iupac("C1CCNCC1") == "piperidine"

    def test_morpholine(self):
        assert smiles_to_iupac("C1CNCCO1") == "morpholine"


# ─── 置換ヘテロ環 ───────────────────────────────────────────────────────

class TestSubstitutedHeterocycles:

    def test_2_methylpyridine(self):
        assert smiles_to_iupac("Cc1ccccn1") == "2-methylpyridine"

    def test_3_chloropyridine(self):
        assert smiles_to_iupac("Clc1cccnc1") == "3-chloropyridine"
