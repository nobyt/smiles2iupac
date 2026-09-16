"""Phase 939: isotopic labeling (IUPAC 2013 P-82) on branched & aryl substituents and multicomponent salt anions.

Extends Phase 938's substituent isotopic labeling to branched alkyl chains (e.g. propan-2-yl),
aryl rings (e.g. phenyl), and carboxylate anions in multicomponent salts.
"""
from smiles2iupac import smiles_to_iupac


class TestPhase939IsotopeBranchedSubstituent:
    def test_deuterated_isopropyl_at_root(self):
        assert smiles_to_iupac("CC(=O)NC([2H])(C)C") == "N-[(2-2H)propan-2-yl]acetamide"

    def test_perdeuterated_isopropyl(self):
        assert smiles_to_iupac("CC(=O)NC([2H])(C([2H])([2H])[2H])C([2H])([2H])[2H]") == "N-[(1,1,1,2,3,3,3-2H7)propan-2-yl]acetamide"

    def test_unlabeled_isopropyl_unaffected(self):
        assert smiles_to_iupac("CC(=O)NC(C)C") == "N-(propan-2-yl)acetamide"


class TestPhase939IsotopeArylSubstituent:
    def test_perdeuterated_phenyl_substituent(self):
        assert smiles_to_iupac("CC(=O)Nc1c([2H])c([2H])c([2H])c([2H])c1[2H]") == "N-[(2,3,4,5,6-2H5)phenyl]acetamide"

    def test_4_deuterated_phenyl_substituent(self):
        assert smiles_to_iupac("CC(=O)Nc1ccc([2H])cc1") == "N-[(4-2H)phenyl]acetamide"

    def test_unlabeled_phenyl_substituent_unaffected(self):
        assert smiles_to_iupac("CC(=O)Nc1ccccc1") == "N-phenylacetamide"


class TestPhase939IsotopeSaltAnion:
    def test_deuterated_sodium_acetate(self):
        assert smiles_to_iupac("[Na+].[2H]C([2H])([2H])C(=O)[O-]") == "sodium (2,2,2-2H3)acetate"

    def test_deuterated_potassium_formate(self):
        assert smiles_to_iupac("[K+].[2H]C(=O)[O-]") == "potassium (1-2H)formate"

    def test_deuterated_sodium_propanoate(self):
        assert smiles_to_iupac("[Na+].[2H]C([2H])([2H])CC(=O)[O-]") == "sodium (3,3,3-2H3)propanoate"

    def test_unlabeled_salt_unaffected(self):
        assert smiles_to_iupac("[Na+].CC(=O)[O-]") == "sodium acetate"
