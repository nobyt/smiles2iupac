"""
Phase 43 テスト: ヒドラゾン (hydrazone)

対象 (IUPAC P-68.3.1):
  C=N-NH₂ または C=N-NHR パターン。
  命名規則: {parent ketone/aldehyde name} hydrazone (2 語)
  例: (CH₃)₂C=N-NH₂ → propan-2-one hydrazone
"""
from smiles2iupac import smiles_to_iupac


class TestKethydrazone:

    def test_propan_2_one_hydrazone(self):
        # (CH₃)₂C=NNH₂
        assert smiles_to_iupac("CC(=NN)C") == "propan-2-one hydrazone"

    def test_butan_2_one_hydrazone(self):
        assert smiles_to_iupac("CCC(=NN)C") == "butan-2-one hydrazone"


class TestAldhydrazone:

    def test_methanal_hydrazone(self):
        # CH₂=N-NH₂
        assert smiles_to_iupac("C=NN") == "methanal hydrazone"

    def test_ethanal_hydrazone(self):
        # CH₃CH=N-NH₂
        assert smiles_to_iupac("CC=NN") == "ethanal hydrazone"

    def test_propanal_hydrazone(self):
        assert smiles_to_iupac("CCC=NN") == "propanal hydrazone"


class TestHydrazoneVsOxime:

    def test_oxime_unchanged(self):
        # オキシムと区別できること
        assert smiles_to_iupac("CC(=NO)C") == "propan-2-one oxime"

    def test_hydrazone_not_alkane(self):
        result = smiles_to_iupac("CC(=NN)C")
        assert "hydrazone" in result
