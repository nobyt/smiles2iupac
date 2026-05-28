"""
Phase 30 テスト: ラクトン (cyclic esters)

対象 (IUPAC P-65.1.1):
  ラクトン = cyclic ester で環内に O を含む。
  命名: Hantzsch-Widman 名 + "-{loc}-one"

  β-ラクトン (4員環): oxetan-2-one
  γ-ラクトン (5員環): oxolan-2-one
  δ-ラクトン (6員環): oxan-2-one
"""
from smiles2iupac import smiles_to_iupac


class TestLactones:

    def test_beta_lactone(self):
        # 4員環: oxetanone
        assert smiles_to_iupac("O=C1CCO1") == "oxetan-2-one"

    def test_gamma_butyrolactone(self):
        # 5員環: oxolanone (γ-butyrolactone)
        assert smiles_to_iupac("O=C1CCCO1") == "oxolan-2-one"

    def test_delta_valerolactone(self):
        # 6員環: oxanone (δ-valerolactone)
        assert smiles_to_iupac("O=C1CCCCO1") == "oxan-2-one"

    def test_alpha_lactone(self):
        # 3員環
        assert smiles_to_iupac("O=C1OC1") == "oxiran-2-one"

    def test_gamma_lactone_alt_smiles(self):
        # 同分子の別 SMILES
        assert smiles_to_iupac("C1COC(=O)C1") == "oxolan-2-one"

    def test_lactam_still_works(self):
        # ラクタム (N含有) は引き続き正常
        assert smiles_to_iupac("O=C1CCCN1") == "pyrrolidin-2-one"
        assert smiles_to_iupac("O=C1CCCCN1") == "piperidin-2-one"
