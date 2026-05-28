"""
Phase 13 テスト: Hantzsch-Widman ヘテロ環

対象:
  - 3員環: oxirane, thiirane, aziridine
  - 4員環: oxetane, thietane, azetidine
  - 5員飽和環: oxolane (THF), thiolane, pyrrolidine
  - 6員飽和環: oxane, thiane
"""
from smiles2iupac import smiles_to_iupac


# ─── 3員環 ─────────────────────────────────────────────────────────────

class TestThreeMemberedRings:

    def test_oxirane(self):
        assert smiles_to_iupac("C1CO1") == "oxirane"

    def test_thiirane(self):
        assert smiles_to_iupac("C1CS1") == "thiirane"

    def test_aziridine(self):
        assert smiles_to_iupac("C1CN1") == "aziridine"


# ─── 4員環 ─────────────────────────────────────────────────────────────

class TestFourMemberedRings:

    def test_oxetane(self):
        assert smiles_to_iupac("C1CCO1") == "oxetane"

    def test_thietane(self):
        assert smiles_to_iupac("C1CCS1") == "thietane"

    def test_azetidine(self):
        assert smiles_to_iupac("C1CCN1") == "azetidine"


# ─── 5員飽和環 ─────────────────────────────────────────────────────────

class TestFiveMemberedSaturatedRings:

    def test_oxolane(self):
        # テトラヒドロフラン (THF)
        assert smiles_to_iupac("C1CCCO1") == "oxolane"

    def test_thiolane(self):
        assert smiles_to_iupac("C1CCCS1") == "thiolane"

    def test_pyrrolidine(self):
        assert smiles_to_iupac("C1CCNC1") == "pyrrolidine"


# ─── 6員飽和環 ─────────────────────────────────────────────────────────

class TestSixMemberedSaturatedRings:

    def test_oxane(self):
        # テトラヒドロピラン (THP)
        assert smiles_to_iupac("C1CCCCO1") == "oxane"

    def test_thiane(self):
        assert smiles_to_iupac("C1CCCCS1") == "thiane"
