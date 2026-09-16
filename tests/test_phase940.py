"""Phase 940: Multiplicative Nomenclature (IUPAC 2013 P-15.3, P-54).

Symmetrical difunctional compounds connected by a divalent linking group
(-O-, -S-, -SO2-, -NH-, -CH2-, -CH2CH2-) are named using multiplicative nomenclature.
"""
from smiles2iupac import smiles_to_iupac


class TestPhase940MultiplicativeDiacids:
    def test_oxydiacetic_acid(self):
        assert smiles_to_iupac("O(CC(=O)O)CC(=O)O") == "2,2'-oxydiacetic acid"

    def test_thiodiacetic_acid(self):
        assert smiles_to_iupac("S(CC(=O)O)CC(=O)O") == "2,2'-thiodiacetic acid"

    def test_iminodiacetic_acid(self):
        assert smiles_to_iupac("N(CC(=O)O)CC(=O)O") == "2,2'-iminodiacetic acid"
        assert smiles_to_iupac("O=C(O)CNCC(=O)O") == "2,2'-iminodiacetic acid"

    def test_oxydipropanoic_acid(self):
        assert smiles_to_iupac("O(CCC(=O)O)CCC(=O)O") == "3,3'-oxydipropanoic acid"

    def test_thiodipropanoic_acid(self):
        assert smiles_to_iupac("S(CCC(=O)O)CCC(=O)O") == "3,3'-thiodipropanoic acid"


class TestPhase940MultiplicativeAromatics:
    def test_methylenedianiline(self):
        assert smiles_to_iupac("Nc1ccc(Cc2ccc(N)cc2)cc1") == "4,4'-methylenedianiline"

    def test_oxydianiline(self):
        assert smiles_to_iupac("Nc1ccc(Oc2ccc(N)cc2)cc1") == "4,4'-oxydianiline"

    def test_thiodianiline(self):
        assert smiles_to_iupac("Nc1ccc(Sc2ccc(N)cc2)cc1") == "4,4'-thiodianiline"

    def test_sulfonyldianiline(self):
        assert smiles_to_iupac("Nc1ccc(S(=O)(=O)c2ccc(N)cc2)cc1") == "4,4'-sulfonyldianiline"

    def test_methylenediphenol(self):
        assert smiles_to_iupac("Oc1ccc(Cc2ccc(O)cc2)cc1") == "4,4'-methylenediphenol"

    def test_ethanediyl_dipyridine(self):
        assert smiles_to_iupac("n1ccccc1CCc2ccccn2") == "2,2'-(ethane-1,2-diyl)dipyridine"
