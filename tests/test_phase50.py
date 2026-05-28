"""
Phase 50 テスト: フェニルエステル / アリールエステル (aryl ester)

対象 (IUPAC P-65.1.1.4):
  R-C(=O)-O-Ph パターン → {phenyl/aryl} {acid}anoate
  例: CC(=O)Oc1ccccc1 → phenyl acetate
"""
from smiles2iupac import smiles_to_iupac


class TestPhenylEster:

    def test_phenyl_ethanoate(self):
        assert smiles_to_iupac("CC(=O)Oc1ccccc1") == "phenyl acetate"

    def test_phenyl_propanoate(self):
        assert smiles_to_iupac("CCC(=O)Oc1ccccc1") == "phenyl propanoate"

    def test_phenyl_methanoate(self):
        assert smiles_to_iupac("O=COc1ccccc1") == "phenyl formate"


class TestSubstitutedPhenylEster:

    def test_4_methylphenyl_ethanoate(self):
        assert smiles_to_iupac("CC(=O)Oc1ccc(C)cc1") == "4-methylphenyl acetate"


class TestPhenylEsterVsAlkylEster:

    def test_ethyl_ethanoate_unchanged(self):
        assert smiles_to_iupac("CCOC(=O)C") == "ethyl acetate"

    def test_ethyl_benzoate_unchanged(self):
        # Phase 47 の回帰テスト: 芳香環が酸側
        assert smiles_to_iupac("CCOC(=O)c1ccccc1") == "ethyl benzoate"

    def test_phenylmethanol_unchanged(self):
        assert smiles_to_iupac("OCc1ccccc1") == "phenylmethanol"
