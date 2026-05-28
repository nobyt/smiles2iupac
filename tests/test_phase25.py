"""
Phase 25 テスト: ラクタム (cyclic amides)

対象 (IUPAC P-66.4.1):
  - 4員環ラクタム: azetidin-2-one (β-lactam)
  - 5員環ラクタム: pyrrolidin-2-one (γ-butyrolactam)
  - 6員環ラクタム: piperidin-2-one (δ-valerolactam)
  - 7員環ラクタム: azepan-2-one (ε-caprolactam)
  命名パターン: {ring_name}-{loc}-one (末尾の "e" を除去してから "-one")
"""
from smiles2iupac import smiles_to_iupac


class TestLactams:

    def test_beta_lactam(self):
        # 4員環: azetidinone
        assert smiles_to_iupac("O=C1CCN1") == "azetidin-2-one"

    def test_gamma_butyrolactam(self):
        # 5員環: pyrrolidinone
        assert smiles_to_iupac("O=C1CCCN1") == "pyrrolidin-2-one"

    def test_delta_valerolactam(self):
        # 6員環: piperidinone
        assert smiles_to_iupac("O=C1CCCCN1") == "piperidin-2-one"

    def test_epsilon_caprolactam(self):
        # 7員環: azepanone
        assert smiles_to_iupac("O=C1CCCCCN1") == "azepan-2-one"

    def test_gamma_butyrolactam_alt_smiles(self):
        # 同じ分子を別 SMILES で
        assert smiles_to_iupac("N1CCCC1=O") == "pyrrolidin-2-one"

    def test_delta_valerolactam_alt_smiles(self):
        assert smiles_to_iupac("N1CCCCC1=O") == "piperidin-2-one"
