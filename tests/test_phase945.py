"""Phase 945: disulfanediyl (-S-S-) multiplicative linker (IUPAC 2013 P-15.3/P-54).

Phase 940 handled single-atom multiplicative linkers (-O-, -S-, -SO2-, -NH-,
-CH2-, -CO-) and the two-carbon -CH2CH2- linker, but not the two-sulfur
disulfanediyl -S-S- linker.  A symmetric disulfide diacid such as
dithiodiglycolic acid (HOOC-CH2-S-S-CH2-COOH) previously fell through to the
generic chain namer and produced the confidently wrong, OPSIN-unparseable
"2-sulfanylethanedioic acid" (dropping half the molecule).

All expected names below are OPSIN round-trip verified.
"""
from smiles2iupac import smiles_to_iupac


class TestPhase945DisulfanediylDiacids:
    def test_dithiodiglycolic_acid(self):
        assert smiles_to_iupac("OC(=O)CSSCC(=O)O") == "2,2'-disulfanediyldiacetic acid"

    def test_dithiodipropanoic_acid(self):
        assert smiles_to_iupac("OC(=O)CCSSCCC(=O)O") == "3,3'-disulfanediyldipropanoic acid"


class TestPhase945DisulfanediylAromatics:
    def test_disulfanediyldianiline(self):
        assert smiles_to_iupac("Nc1ccc(SSc2ccc(N)cc2)cc1") == "4,4'-disulfanediyldianiline"


class TestPhase945NoRegression:
    def test_thiodiacetic_still_works(self):
        assert smiles_to_iupac("OC(=O)CSCC(=O)O") == "2,2'-thiodiacetic acid"

    def test_plain_dicarboxylic_acid_not_multiplicative(self):
        # -S-S- absent: a plain saturated diacid must NOT become multiplicative
        assert smiles_to_iupac("OC(=O)CCCCCCC(=O)O") == "octanedioic acid"
