"""
Phase 7 テスト: 含窒素官能基 + 既存改善

対象:
  - 第一級アミン (primary amine) → -amine suffix
  - ニトリル (nitrile)           → -nitrile suffix
  - ニトロ基 (nitro)             → nitro- prefix
  - aniline / benzonitrile (芳香族保留名)
  - cyclohexanamine / cyclohexanecarbonitrile (環状)
  - 改善 #1: cyclohexanecarbaldehyde (旧 cyclohexanal を修正)
  - 改善 #2: bis/tris 接頭辞
"""
import pytest
from smiles2iupac import smiles_to_iupac


# ─── 第一級アミン (脂肪族) ──────────────────────────────────────────

class TestPrimaryAmines:

    def test_methanamine(self):
        # CH₃-NH₂
        assert smiles_to_iupac("CN") == "methanamine"

    def test_ethanamine(self):
        # CH₃-CH₂-NH₂
        assert smiles_to_iupac("CCN") == "ethanamine"

    def test_propan_1_amine(self):
        # CH₃-CH₂-CH₂-NH₂
        assert smiles_to_iupac("CCCN") == "propan-1-amine"

    def test_propan_2_amine(self):
        # CH₃-CH(NH₂)-CH₃
        assert smiles_to_iupac("CC(N)C") == "propan-2-amine"

    def test_butan_1_amine(self):
        assert smiles_to_iupac("CCCCN") == "butan-1-amine"

    def test_butan_2_amine(self):
        assert smiles_to_iupac("CCC(N)C") == "butan-2-amine"

    def test_2_methylpropan_1_amine(self):
        # (CH₃)₂CHCH₂NH₂ → 2-methylpropan-1-amine
        assert smiles_to_iupac("CC(C)CN") == "2-methylpropan-1-amine"


# ─── ニトリル (脂肪族) ─────────────────────────────────────────────

class TestNitriles:

    def test_methanenitrile(self):
        # H-C≡N (formonitrile / methanenitrile)
        assert smiles_to_iupac("C#N") == "methanenitrile"

    def test_acetonitrile(self):
        # CH₃-C≡N = acetonitrile (IUPAC 2013 retained name)
        assert smiles_to_iupac("CC#N") == "acetonitrile"

    def test_propanenitrile(self):
        # CH₃-CH₂-C≡N
        assert smiles_to_iupac("CCC#N") == "propanenitrile"

    def test_butanenitrile(self):
        assert smiles_to_iupac("CCCC#N") == "butanenitrile"

    def test_2_methylpropanenitrile(self):
        # (CH₃)₂CH-C≡N → 2-methylpropanenitrile
        assert smiles_to_iupac("CC(C)C#N") == "2-methylpropanenitrile"

    def test_3_methylbutanenitrile(self):
        # CH₃CH(CH₃)CH₂C≡N → 3-methylbutanenitrile
        assert smiles_to_iupac("CC(C)CC#N") == "3-methylbutanenitrile"


# ─── 芳香族アミン / ニトリル ────────────────────────────────────────

class TestAromaticNitrogenGroups:

    def test_aniline(self):
        # PhNH₂ → aniline (保留名)
        assert smiles_to_iupac("Nc1ccccc1") == "aniline"

    def test_benzonitrile(self):
        # Ph-C≡N → benzonitrile (保留名)
        assert smiles_to_iupac("N#Cc1ccccc1") == "benzonitrile"

    def test_4_methylaniline(self):
        # 4-CH₃-C₆H₄-NH₂ → 4-methylaniline
        assert smiles_to_iupac("Cc1ccc(N)cc1") == "4-methylaniline"

    def test_4_chloroaniline(self):
        assert smiles_to_iupac("Clc1ccc(N)cc1") == "4-chloroaniline"

    def test_4_methylbenzonitrile(self):
        assert smiles_to_iupac("Cc1ccc(C#N)cc1") == "4-methylbenzonitrile"


# ─── ニトロ基 (置換基接頭辞) ────────────────────────────────────────

class TestNitroGroup:

    def test_nitrobenzene(self):
        assert smiles_to_iupac("c1ccc([N+](=O)[O-])cc1") == "nitrobenzene"

    def test_1_nitropropane(self):
        assert smiles_to_iupac("CCC[N+](=O)[O-]") == "1-nitropropane"

    def test_2_nitropropane(self):
        assert smiles_to_iupac("CC([N+](=O)[O-])C") == "2-nitropropane"

    def test_4_nitrotoluene(self):
        # 4-nitrotoluene = 1-methyl-4-nitrobenzene
        assert smiles_to_iupac("Cc1ccc([N+](=O)[O-])cc1") == "1-methyl-4-nitrobenzene"


# ─── 環状アミン / ニトリル ──────────────────────────────────────────

class TestCyclicNitrogenGroups:

    def test_cyclohexanamine(self):
        # cyclohexylamine → cyclohexanamine (IUPAC PIN)
        assert smiles_to_iupac("NC1CCCCC1") == "cyclohexanamine"

    def test_cyclohexanecarbonitrile(self):
        # N#C-c6H11 → cyclohexanecarbonitrile
        assert smiles_to_iupac("N#CC1CCCCC1") == "cyclohexanecarbonitrile"


# ─── 改善 #1: carbaldehyde 形式 ─────────────────────────────────────

class TestCarbaldehyde:

    def test_cyclohexanecarbaldehyde(self):
        # O=C-c6H11 → cyclohexanecarbaldehyde (旧: cyclohexanal)
        assert smiles_to_iupac("O=CC1CCCCC1") == "cyclohexanecarbaldehyde"

    def test_cyclopentanecarbaldehyde(self):
        assert smiles_to_iupac("O=CC1CCCC1") == "cyclopentanecarbaldehyde"


# ─── 改善 #2: bis/tris 接頭辞 ────────────────────────────────────────

class TestBisTris:

    def test_bis_propan_2_yl(self):
        # C4 に propan-2-yl (isopropyl) が 2 つ付いた heptane
        # 4,4-bis(propan-2-yl)heptane: CCCC(C(C)C)(C(C)C)CCC
        result = smiles_to_iupac("CCCC(C(C)C)(C(C)C)CCC")
        # 置換基 "propan-2-yl" が 2 つ → bis(propan-2-yl)
        assert "bis(propan-2-yl)" in result
