"""
Phase 31 テスト: ナフタレン官能基誘導体

対象 (IUPAC P-31.1.2):
  ナフタレンに官能基 suffix (-ol, -amine, -carboxylic acid, -carbaldehyde 等) を付加。
  命名: naphthalen-{n}-ol, naphthalen-{n}-amine, naphthalene-{n}-carboxylic acid, etc.
"""
from smiles2iupac import smiles_to_iupac


class TestNaphthaleneOlAmine:

    def test_naphthalen_2_ol(self):
        assert smiles_to_iupac("Oc1ccc2ccccc2c1") == "naphthalen-2-ol"

    def test_naphthalen_1_ol(self):
        assert smiles_to_iupac("Oc1cccc2ccccc12") == "naphthalen-1-ol"

    def test_naphthalen_2_amine(self):
        assert smiles_to_iupac("Nc1ccc2ccccc2c1") == "naphthalen-2-amine"

    def test_naphthalen_1_amine(self):
        assert smiles_to_iupac("Nc1cccc2ccccc12") == "naphthalen-1-amine"


class TestNaphthaleneExocyclicFG:

    def test_naphthalene_2_carboxylic_acid(self):
        assert smiles_to_iupac("OC(=O)c1ccc2ccccc2c1") == "naphthalene-2-carboxylic acid"

    def test_naphthalene_1_carboxylic_acid(self):
        assert smiles_to_iupac("OC(=O)c1cccc2ccccc12") == "naphthalene-1-carboxylic acid"

    def test_naphthalene_2_carbaldehyde(self):
        assert smiles_to_iupac("O=Cc1ccc2ccccc2c1") == "naphthalene-2-carbaldehyde"

    def test_naphthalene_1_carbaldehyde(self):
        assert smiles_to_iupac("O=Cc1cccc2ccccc12") == "naphthalene-1-carbaldehyde"

    def test_naphthalene_2_carboxamide(self):
        assert smiles_to_iupac("NC(=O)c1ccc2ccccc2c1") == "naphthalene-2-carboxamide"
