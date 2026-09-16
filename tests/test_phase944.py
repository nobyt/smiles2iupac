"""Phase 944: general ylidene naming for phosphonium/sulfonium ylides (IUPAC 2013 P-74).

Phase 943 introduced ylide naming but its hand-rolled ``_name_ylidene`` only
recognised methylidene / ethylidene / propylidene / propan-2-ylidene and
silently collapsed every other ylidene carbon to "methylidene", dropping atoms
(e.g. a butylidene Wittig reagent was mis-named as the methylidene one).

Phase 944 replaces that with the shared, general ``_alkylidene_name`` so any
chain length, branching, aryl (benzylidene) or on-chain substituent is named
correctly.  All expected names below are OPSIN round-trip verified.
"""
from smiles2iupac import smiles_to_iupac


class TestPhase944LongChainYlidenes:
    def test_butylidene_triphenylphosphorane_double_bond(self):
        assert smiles_to_iupac("CCCC=P(c1ccccc1)(c2ccccc2)c3ccccc3") == "(butan-1-ylidene)triphenyl-lambda5-phosphane"

    def test_butylidene_triphenylphosphorane_zwitterionic(self):
        assert smiles_to_iupac("[P+](c1ccccc1)(c2ccccc2)(c3ccccc3)[CH-]CCC") == "(butan-1-ylidene)triphenyl-lambda5-phosphane"

    def test_propylidene_triphenylphosphorane(self):
        assert smiles_to_iupac("CCC=P(c1ccccc1)(c2ccccc2)c3ccccc3") == "triphenyl(propan-1-ylidene)-lambda5-phosphane"

    def test_butylidene_trimethylphosphorane(self):
        assert smiles_to_iupac("CCCC=P(C)(C)C") == "(butan-1-ylidene)trimethyl-lambda5-phosphane"


class TestPhase944BranchedYlidenes:
    def test_2_methylpropylidene_triphenylphosphorane(self):
        assert smiles_to_iupac("CC(C)C=P(c1ccccc1)(c2ccccc2)c3ccccc3") == "(2-methylpropan-1-ylidene)triphenyl-lambda5-phosphane"


class TestPhase944ArylYlidenes:
    def test_benzylidene_triphenylphosphorane_double_bond(self):
        assert smiles_to_iupac("c1ccccc1C=P(c2ccccc2)(c3ccccc3)c4ccccc4") == "(benzylidene)triphenyl-lambda5-phosphane"

    def test_benzylidene_triphenylphosphorane_zwitterionic(self):
        assert smiles_to_iupac("[P+](c1ccccc1)(c2ccccc2)(c3ccccc3)[CH-]c4ccccc4") == "(benzylidene)triphenyl-lambda5-phosphane"


class TestPhase944SulfoniumYlides:
    def test_butylidene_dimethylsulfane_double_bond(self):
        assert smiles_to_iupac("CCCC=S(C)C") == "(butan-1-ylidene)dimethyl-lambda4-sulfane"

    def test_butylidene_dimethylsulfane_zwitterionic(self):
        assert smiles_to_iupac("C[S+](C)[CH-]CCC") == "(butan-1-ylidene)dimethyl-lambda4-sulfane"

    def test_2_2_dimethylpropylidene_dimethylsulfane(self):
        assert smiles_to_iupac("CC(C)(C)C=S(C)C") == "(2,2-dimethylpropan-1-ylidene)dimethyl-lambda4-sulfane"
