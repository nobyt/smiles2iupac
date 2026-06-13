"""
Phase 17 テスト: 縮合ヘテロ芳香族保留名

対象:
  - quinoline, isoquinoline, 1H-indole, 1H-benzimidazole
  - benzofuran, benzo[b]thiophene, purine, acridine
"""
from smiles2iupac import smiles_to_iupac


class TestFusedHeteroAromaticRetainedNames:

    def test_quinoline(self):
        assert smiles_to_iupac("c1ccc2ncccc2c1") == "quinoline"

    def test_isoquinoline(self):
        assert smiles_to_iupac("c1ccc2cnccc2c1") == "isoquinoline"

    def test_indole(self):
        assert smiles_to_iupac("c1ccc2[nH]ccc2c1") == "1H-indole"

    def test_benzimidazole(self):
        assert smiles_to_iupac("c1ccc2[nH]cnc2c1") == "1H-benzimidazole"

    def test_benzofuran(self):
        assert smiles_to_iupac("c1ccc2occc2c1") == "benzofuran"

    def test_benzothiophene(self):
        assert smiles_to_iupac("c1ccc2sccc2c1") == "benzo[b]thiophene"

    def test_purine(self):
        assert smiles_to_iupac("c1ncc2[nH]cnc2n1") == "7H-purine"

    def test_acridine(self):
        assert smiles_to_iupac("c1ccc2nc3ccccc3cc2c1") == "acridine"
