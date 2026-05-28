"""
Phase 41 テスト: チオアミド (thioamide)

対象 (IUPAC P-65.3.3):
  C(=S)-NR₂ パターン。
  命名規則: {stem}anethioamide (primary) / N-{sub}{stem}anethioamide (secondary/tertiary)
  例: CH₃C(=S)NH₂ → ethanethioamide
"""
from smiles2iupac import smiles_to_iupac


class TestPrimaryThioamide:

    def test_methanethioamide(self):
        assert smiles_to_iupac("NC(=S)") == "methanethioamide"

    def test_ethanethioamide(self):
        assert smiles_to_iupac("CC(=S)N") == "ethanethioamide"

    def test_propanethoamide(self):
        assert smiles_to_iupac("CCC(=S)N") == "propanethioamide"

    def test_butanethioamide(self):
        assert smiles_to_iupac("CCCC(=S)N") == "butanethioamide"


class TestNSubstitutedThioamide:

    def test_n_methylethanethioamide(self):
        assert smiles_to_iupac("CC(=S)NC") == "N-methylethanethioamide"

    def test_nn_dimethylethanethioamide(self):
        assert smiles_to_iupac("CC(=S)N(C)C") == "N,N-dimethylethanethioamide"

    def test_n_ethylethanethioamide(self):
        assert smiles_to_iupac("CC(=S)NCC") == "N-ethylethanethioamide"


class TestThioamideNotAmine:

    def test_does_not_produce_amine(self):
        result = smiles_to_iupac("CC(=S)N")
        assert "amine" not in result
        assert result == "ethanethioamide"
