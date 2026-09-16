"""Phase 943: Phosphonium & Sulfonium Ylides (IUPAC 2013 P-74, P-67, P-68).

Wittig reagents and sulfonium ylides in both double-bond (hypervalent)
and zwitterionic (formal charge) representations are named using IUPAC lambda nomenclature.
"""
from smiles2iupac import smiles_to_iupac


class TestPhase943PhosphoniumYlides:
    def test_ethylidene_triphenylphosphorane_double_bond(self):
        assert smiles_to_iupac("CC=[P](c1ccccc1)(c2ccccc2)c3ccccc3") == "(ethylidene)triphenyl-lambda5-phosphane"

    def test_ethylidene_triphenylphosphorane_zwitterionic(self):
        assert smiles_to_iupac("[P+](c1ccccc1)(c2ccccc2)(c3ccccc3)[CH-]C") == "(ethylidene)triphenyl-lambda5-phosphane"

    def test_methylidene_triphenylphosphorane_double_bond(self):
        assert smiles_to_iupac("C=[P](c1ccccc1)(c2ccccc2)c3ccccc3") == "(methylidene)triphenyl-lambda5-phosphane"

    def test_methylidene_triphenylphosphorane_zwitterionic(self):
        assert smiles_to_iupac("[P+](c1ccccc1)(c2ccccc2)(c3ccccc3)[CH2-]") == "(methylidene)triphenyl-lambda5-phosphane"

    def test_isopropylidene_triphenylphosphorane_zwitterionic(self):
        assert smiles_to_iupac("[P+](c1ccccc1)(c2ccccc2)(c3ccccc3)[C-](C)C") == "triphenyl(propan-2-ylidene)-lambda5-phosphane"


class TestPhase943SulfoniumYlides:
    def test_dimethyl_isopropylidene_sulfane_zwitterionic(self):
        assert smiles_to_iupac("[S+](C)(C)[C-](C)C") == "dimethyl(propan-2-ylidene)-lambda4-sulfane"

    def test_dimethyl_methylidene_sulfane_zwitterionic(self):
        assert smiles_to_iupac("[S+](C)(C)[CH2-]") == "dimethyl(methylidene)-lambda4-sulfane"

    def test_dimethyl_isopropylidene_sulfane_double_bond(self):
        assert smiles_to_iupac("CC(C)=[S](C)C") == "dimethyl(propan-2-ylidene)-lambda4-sulfane"
