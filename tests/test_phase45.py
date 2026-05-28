"""
Phase 45 テスト: 検出バグ修正

修正内容:
  45-a: 環状炭酸エステル (O=C1OCCO1) → 誤って "carbonate" 検出 → 修正
  45-b: ペルオキシ酸 (CC(=O)OO) → 誤って "hydroperoxide" 検出 → 修正
  45-c: ジアゾメタン (C=[N+]=[N-]) → 誤って "hydrazone" 検出 → 修正
"""
from smiles2iupac import smiles_to_iupac


class TestCyclicCarbonateNotDetectedAsCarbonate:

    def test_cyclic_carbonate_not_labeled_carbonate(self):
        # 4員環状炭酸エステル: 炭酸エステルとして命名しない (lactone path → oxetan-2-one)
        result = smiles_to_iupac("O=C1OCC1")
        assert "carbonate" not in result

    def test_open_chain_carbonate_still_works(self):
        # 鎖状炭酸エステルは引き続き carbonate として命名される
        assert smiles_to_iupac("COC(=O)OC") == "dimethyl carbonate"


class TestPeroxoacidNotHydroperoxide:

    def test_peracetic_acid_not_hydroperoxide(self):
        # CH₃C(=O)OO → ペルオキシ酢酸; hydroperoxide ではない
        result = smiles_to_iupac("CC(=O)OO")
        assert "hydroperoxide" not in result

    def test_hydroperoxide_still_works(self):
        # 通常の hydroperoxide は変わらない
        assert smiles_to_iupac("CCOO") == "ethyl hydroperoxide"


class TestDiazomethaneNotHydrazone:

    def test_diazomethane_not_hydrazone(self):
        # C=[N+]=[N-] → ジアゾメタン; hydrazone ではない
        result = smiles_to_iupac("C=[N+]=[N-]")
        assert "hydrazone" not in result

    def test_hydrazone_still_works(self):
        # 通常の hydrazone は変わらない
        assert smiles_to_iupac("CC(=NN)C") == "propan-2-one hydrazone"

    def test_aldhydrazone_still_works(self):
        assert smiles_to_iupac("CC=NN") == "ethanal hydrazone"
