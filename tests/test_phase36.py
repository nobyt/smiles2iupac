"""
Phase 36 テスト: N-アリールアミン (N-arylamine)

対象 (IUPAC P-62.2.3):
  N が芳香環の炭素に直結している二級・三級アミン。
  - 単一フェニル: aniline を親名とし N-置換基を列挙
  - N,N-二置換: N,N-dimethylaniline 等
"""
from smiles2iupac import smiles_to_iupac


class TestNMethylaniline:

    def test_n_methylaniline(self):
        assert smiles_to_iupac("CNc1ccccc1") == "N-methylaniline"

    def test_n_methylaniline_reversed(self):
        assert smiles_to_iupac("c1ccc(NC)cc1") == "N-methylaniline"


class TestNNDimethylaniline:

    def test_nn_dimethylaniline(self):
        assert smiles_to_iupac("CN(C)c1ccccc1") == "N,N-dimethylaniline"


class TestNEthylaniline:

    def test_n_ethylaniline(self):
        assert smiles_to_iupac("CCNc1ccccc1") == "N-ethylaniline"


class TestDiphenylamine:

    def test_diphenylamine(self):
        # IUPAC 2013 P-62.2.3.2: diphenylamine は保留優先名 (PIN)
        assert smiles_to_iupac("c1ccc(Nc2ccccc2)cc1") == "diphenylamine"
