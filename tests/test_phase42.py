"""
Phase 42 テスト: カルバメート (carbamate / urethane)

対象 (IUPAC P-65.1.2.5):
  N-C(=O)-O-R パターン（アミドとエステルの複合官能基）。
  命名規則: {alkyl} carbamate / {alkyl} N-{sub}carbamate
  例: NH₂C(=O)OCH₃ → methyl carbamate
"""
from smiles2iupac import smiles_to_iupac


class TestPrimaryCarbamate:

    def test_methyl_carbamate(self):
        assert smiles_to_iupac("NC(=O)OC") == "methyl carbamate"

    def test_ethyl_carbamate(self):
        assert smiles_to_iupac("NC(=O)OCC") == "ethyl carbamate"

    def test_propyl_carbamate(self):
        assert smiles_to_iupac("NC(=O)OCCC") == "propyl carbamate"


class TestNSubstitutedCarbamate:

    def test_methyl_n_methylcarbamate(self):
        assert smiles_to_iupac("CNC(=O)OC") == "methyl N-methylcarbamate"

    def test_ethyl_n_methylcarbamate(self):
        assert smiles_to_iupac("CNC(=O)OCC") == "ethyl N-methylcarbamate"

    def test_methyl_nn_dimethylcarbamate(self):
        assert smiles_to_iupac("CN(C)C(=O)OC") == "methyl N,N-dimethylcarbamate"


class TestCarbamateVsEster:

    def test_plain_ester_unchanged(self):
        # ester (N なし) はカルバメートと区別される
        assert smiles_to_iupac("CC(=O)OC") == "methyl acetate"

    def test_carbamate_not_ester(self):
        result = smiles_to_iupac("NC(=O)OC")
        assert result == "methyl carbamate"
        assert "methanoate" not in result
