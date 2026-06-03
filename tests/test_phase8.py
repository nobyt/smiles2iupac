"""
Phase 8 テスト: チオール・エーテル・スルフィド

対象:
  - チオール (thiol) → -thiol suffix
  - エーテル (ether) → (R)oxy prefix
  - スルフィド (sulfide) → (R)sulfanyl prefix
  - 芳香族チオール / エーテル (benzenethiol, methoxybenzene)
"""
from smiles2iupac import smiles_to_iupac


# ─── チオール (脂肪族) ─────────────────────────────────────────────────

class TestAliphaticThiols:

    def test_methanethiol(self):
        # CH₃-SH
        assert smiles_to_iupac("CS") == "methanethiol"

    def test_ethanethiol(self):
        # CH₃-CH₂-SH
        assert smiles_to_iupac("CCS") == "ethanethiol"

    def test_propane_1_thiol(self):
        # CH₃-CH₂-CH₂-SH
        assert smiles_to_iupac("CCCS") == "propane-1-thiol"

    def test_propane_2_thiol(self):
        # CH₃-CH(SH)-CH₃
        assert smiles_to_iupac("CC(S)C") == "propane-2-thiol"

    def test_butane_1_thiol(self):
        assert smiles_to_iupac("CCCCS") == "butane-1-thiol"

    def test_butane_2_thiol(self):
        assert smiles_to_iupac("CC(S)CC") == "butane-2-thiol"

    def test_2_methylpropane_1_thiol(self):
        # (CH₃)₂CHCH₂SH
        assert smiles_to_iupac("CC(C)CS") == "2-methylpropane-1-thiol"


# ─── 芳香族チオール ────────────────────────────────────────────────────

class TestAromaticThiols:

    def test_benzenethiol(self):
        # Ph-SH → benzenethiol (IUPAC 保留名)
        assert smiles_to_iupac("Sc1ccccc1") == "benzenethiol"

    def test_4_methylbenzenethiol(self):
        assert smiles_to_iupac("Cc1ccc(S)cc1") == "4-methylbenzenethiol"


# ─── 環状チオール ──────────────────────────────────────────────────────

class TestCyclicThiols:

    def test_cyclohexane_1_thiol(self):
        assert smiles_to_iupac("SC1CCCCC1") == "cyclohexanethiol"

    def test_cyclopentane_1_thiol(self):
        assert smiles_to_iupac("SC1CCCC1") == "cyclopentanethiol"


# ─── エーテル (prefix: (R)oxy) ──────────────────────────────────────────

class TestEthers:

    def test_methoxymethane(self):
        # CH₃-O-CH₃ (dimethyl ether): 1C 鎖のためロカント付き
        assert smiles_to_iupac("COC") == "methoxymethane"

    def test_1_methoxyethane(self):
        # CH₃-O-CH₂-CH₃ → ロカント省略 (IUPAC 2013, 2炭素鎖単一置換)
        assert smiles_to_iupac("CCOC") == "methoxyethane"

    def test_1_methoxypropane(self):
        # CH₃-O-CH₂-CH₂-CH₃
        assert smiles_to_iupac("CCCOC") == "1-methoxypropane"

    def test_1_ethoxyethane(self):
        # CH₃-CH₂-O-CH₂-CH₃ → ロカント省略 (IUPAC 2013, 2炭素鎖単一置換)
        assert smiles_to_iupac("CCOCC") == "ethoxyethane"

    def test_methoxybenzene(self):
        # Ph-O-CH₃ (anisole)
        assert smiles_to_iupac("COc1ccccc1") == "methoxybenzene"

    def test_ethoxybenzene(self):
        # Ph-O-CH₂-CH₃ (phenetole)
        assert smiles_to_iupac("CCOc1ccccc1") == "ethoxybenzene"


# ─── スルフィド (prefix: (R)sulfanyl) ──────────────────────────────────

class TestSulfides:

    def test_dimethyl_sulfide(self):
        # CH₃-S-CH₃ → dimethyl sulfide (IUPAC 2013 P-63.6.1.1 substitutive PIN)
        assert smiles_to_iupac("CSC") == "dimethyl sulfide"

    def test_methyl_phenyl_sulfide(self):
        # Ph-S-CH₃ → methyl phenyl sulfide (substitutive PIN)
        assert smiles_to_iupac("CSc1ccccc1") == "methyl phenyl sulfide"
