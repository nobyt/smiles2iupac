"""
Phase 1 テスト: 直鎖・分岐アルカン、アルコール、ハロゲン化物、アルケン、アルキン

IUPAC 2013 Blue Book の系統名を基準にする。
  - アルコール: ethanol は保留名 (PIN), propan-2-ol 等はロカント付き
  - アルケン・アルキンは chain_length==2 の場合のみロカント省略 (ethene, ethyne)
"""

import pytest
from smiles2iupac import smiles_to_iupac


# ─── アルカン ────────────────────────────────────────────────────────

class TestAlkanes:

    def test_methane(self):
        assert smiles_to_iupac("C") == "methane"

    def test_ethane(self):
        assert smiles_to_iupac("CC") == "ethane"

    def test_propane(self):
        assert smiles_to_iupac("CCC") == "propane"

    def test_butane(self):
        assert smiles_to_iupac("CCCC") == "butane"

    def test_hexane(self):
        assert smiles_to_iupac("CCCCCC") == "hexane"

    def test_2_methylpropane(self):
        assert smiles_to_iupac("CC(C)C") == "2-methylpropane"

    def test_2_methylbutane(self):
        assert smiles_to_iupac("CC(C)CC") == "2-methylbutane"

    def test_3_methylpentane(self):
        assert smiles_to_iupac("CCC(C)CC") == "3-methylpentane"

    def test_2_2_dimethylpropane(self):
        assert smiles_to_iupac("CC(C)(C)C") == "2,2-dimethylpropane"

    def test_2_3_dimethylbutane(self):
        assert smiles_to_iupac("CC(C)C(C)C") == "2,3-dimethylbutane"


# ─── アルコール ─────────────────────────────────────────────────────

class TestAlcohols:

    def test_ethan_1_ol(self):
        # IUPAC 2013 保留名 (PIN): ethanol
        assert smiles_to_iupac("CCO") == "ethanol"

    def test_propan_1_ol(self):
        assert smiles_to_iupac("CCCO") == "propan-1-ol"

    def test_propan_2_ol(self):
        assert smiles_to_iupac("CC(O)C") == "propan-2-ol"

    def test_pentan_2_ol(self):
        assert smiles_to_iupac("CC(O)CCC") == "pentan-2-ol"

    def test_butan_1_ol(self):
        assert smiles_to_iupac("CCCCO") == "butan-1-ol"


# ─── ハロゲン化物 ────────────────────────────────────────────────────

class TestHalides:

    def test_1_chloroethane(self):
        # CCCl = C-C-Cl = 2炭素、単一置換 → ロカント省略 (IUPAC 2013)
        assert smiles_to_iupac("CCCl") == "chloroethane"

    def test_1_chloropropane(self):
        # CCCCl = C-C-C-Cl = 3炭素 → 1-chloropropane
        assert smiles_to_iupac("CCCCl") == "1-chloropropane"

    def test_2_chloropropane(self):
        assert smiles_to_iupac("CC(Cl)C") == "2-chloropropane"

    def test_1_bromopropane(self):
        # CCCBr = C-C-C-Br = 3炭素 → 1-bromopropane
        assert smiles_to_iupac("CCCBr") == "1-bromopropane"

    def test_1_fluoroethane(self):
        assert smiles_to_iupac("CCF") == "fluoroethane"

    def test_4_chloro_2_methylhexane(self):
        # 置換基アルファベット順: chloro < methyl
        assert smiles_to_iupac("CCC(Cl)CC(C)C") == "4-chloro-2-methylhexane"


# ─── アルケン ────────────────────────────────────────────────────────

class TestAlkenes:

    def test_ethene(self):
        # 2C: ロカント省略
        assert smiles_to_iupac("C=C") == "ethene"

    def test_prop_1_ene(self):
        assert smiles_to_iupac("C=CC") == "propene"

    def test_but_2_ene(self):
        assert smiles_to_iupac("CC=CC") == "but-2-ene"

    def test_hex_1_ene(self):
        assert smiles_to_iupac("C=CCCCC") == "hex-1-ene"


# ─── アルキン ────────────────────────────────────────────────────────

class TestAlkynes:

    def test_ethyne(self):
        # 2C: ロカント省略
        assert smiles_to_iupac("C#C") == "ethyne"

    def test_but_2_yne(self):
        assert smiles_to_iupac("CC#CC") == "but-2-yne"

    def test_prop_1_yne(self):
        assert smiles_to_iupac("C#CC") == "propyne"
