"""
Phase 61 テスト: ヒドロキサム酸 / ジェミナルジオール

ヒドロキサム酸 (IUPAC P-66.4.1.4):
  R-C(=O)-NH-OH パターン → N-hydroxy{amide}
  例: CC(=O)NO → N-hydroxyacetamide

ジェミナルジオール locant 修正:
  OC(O)C → ethane-1,1-diol (1,2 ではなく 1,1)
"""
from smiles2iupac import smiles_to_iupac


class TestHydroxamicAcid:

    def test_n_hydroxyacetamide(self):
        assert smiles_to_iupac("CC(=O)NO") == "N-hydroxyacetamide"

    def test_n_hydroxymethnamide(self):
        assert smiles_to_iupac("C(=O)NO") == "N-hydroxyformamide"

    def test_n_hydroxypropanamide(self):
        assert smiles_to_iupac("CCC(=O)NO") == "N-hydroxypropanamide"


class TestGeminalDiol:

    def test_ethane_1_1_diol(self):
        assert smiles_to_iupac("OC(O)C") == "ethane-1,1-diol"

    def test_ethane_1_2_diol_unchanged(self):
        assert smiles_to_iupac("OCCO") == "ethane-1,2-diol"

    def test_propane_1_2_diol_unchanged(self):
        assert smiles_to_iupac("OC(O)CO") != "propane-1,1,2-triol"  # triol は別
        assert smiles_to_iupac("OCC(O)C") == "propane-1,2-diol"
