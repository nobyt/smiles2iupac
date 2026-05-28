"""
Phase 58 テスト: 環状イミド (cyclic imide)

対象 (IUPAC P-66.4.2):
  N が両側の C(=O) に挟まれた環 → {base}-{loc1},{loc2}-dione
  例: O=C1NC(=O)CC1 → pyrrolidine-2,5-dione
"""
from smiles2iupac import smiles_to_iupac


class TestCyclicImide:

    def test_pyrrolidine_2_5_dione(self):
        assert smiles_to_iupac("O=C1NC(=O)CC1") == "pyrrolidine-2,5-dione"

    def test_piperidine_2_6_dione(self):
        assert smiles_to_iupac("O=C1NC(=O)CCC1") == "piperidine-2,6-dione"

    def test_azepane_2_7_dione(self):
        assert smiles_to_iupac("O=C1NC(=O)CCCC1") == "azepane-2,7-dione"


class TestCyclicImideVsLactam:

    def test_pyrrolidinone_unchanged(self):
        # lactam: C=O は N の片側のみ
        assert smiles_to_iupac("O=C1CCCN1") == "pyrrolidin-2-one"

    def test_piperidinone_unchanged(self):
        assert smiles_to_iupac("O=C1CCCCN1") == "piperidin-2-one"

    def test_azepanone_unchanged(self):
        assert smiles_to_iupac("O=C1CCCCCN1") == "azepan-2-one"
