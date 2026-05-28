"""
Phase 29 テスト: アレン・累積ジエン・共役ジエン

対象 (IUPAC P-31.1.3):
  累積ジエン (allene): C=C=C → propa-1,2-diene
  共役ジエン: C=CC=C → buta-1,3-diene
  累積トリエン: C=C=C=C → buta-1,2,3-triene

命名パターン:
  - 1つの C=C: {stem}-{loc}-ene
  - 2つの C=C: {stem}a-{locs}-diene  (allene, conjugated)
  - 3つの C=C: {stem}a-{locs}-triene
"""
from smiles2iupac import smiles_to_iupac


class TestAllenes:

    def test_propadiene(self):
        # C=C=C (allene)
        assert smiles_to_iupac("C=C=C") == "propa-1,2-diene"

    def test_buta_1_2_diene(self):
        # CH₂=C=CH-CH₃
        assert smiles_to_iupac("C=C=CC") == "buta-1,2-diene"

    def test_penta_2_3_diene(self):
        # CH₃-CH=C=CH-CH₃
        assert smiles_to_iupac("CC=C=CC") == "penta-2,3-diene"


class TestConjugatedDienes:

    def test_buta_1_3_diene(self):
        # CH₂=CH-CH=CH₂
        assert smiles_to_iupac("C=CC=C") == "buta-1,3-diene"

    def test_penta_1_4_diene(self):
        # CH₂=CH-CH₂-CH=CH₂
        assert smiles_to_iupac("C=CCC=C") == "penta-1,4-diene"

    def test_hexa_1_3_diene(self):
        assert smiles_to_iupac("C=CC=CCC") == "hexa-1,3-diene"


class TestTrienes:

    def test_buta_1_2_3_triene(self):
        # C=C=C=C (butatriene)
        assert smiles_to_iupac("C=C=C=C") == "buta-1,2,3-triene"
