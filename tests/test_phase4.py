"""
Phase 4 テスト: ベンゼン誘導体 (官能基付き)

対象:
  - フェノール類 (phenol, 置換 phenol)
  - 安息香酸類 (benzoic acid, 置換 benzoic acid)
  - ベンズアルデヒド類 (benzaldehyde, 置換 benzaldehyde)
"""
import pytest
from smiles2iupac import smiles_to_iupac


# ─── フェノール ──────────────────────────────────────────────────────

class TestPhenol:

    def test_phenol(self):
        assert smiles_to_iupac("Oc1ccccc1") == "phenol"

    def test_4_chlorophenol(self):
        assert smiles_to_iupac("Oc1ccc(Cl)cc1") == "4-chlorophenol"

    def test_2_methylphenol(self):
        # o-cresol
        assert smiles_to_iupac("Oc1ccccc1C") == "2-methylphenol"

    def test_4_methylphenol(self):
        # p-cresol
        assert smiles_to_iupac("Oc1ccc(C)cc1") == "4-methylphenol"

    def test_4_bromophenol(self):
        assert smiles_to_iupac("Oc1ccc(Br)cc1") == "4-bromophenol"


# ─── 安息香酸 ────────────────────────────────────────────────────────

class TestBenzoicAcid:

    def test_benzoic_acid(self):
        assert smiles_to_iupac("OC(=O)c1ccccc1") == "benzoic acid"

    def test_4_methylbenzoic_acid(self):
        assert smiles_to_iupac("OC(=O)c1ccc(C)cc1") == "4-methylbenzoic acid"

    def test_4_chlorobenzoic_acid(self):
        assert smiles_to_iupac("OC(=O)c1ccc(Cl)cc1") == "4-chlorobenzoic acid"

    def test_2_chlorobenzoic_acid(self):
        assert smiles_to_iupac("OC(=O)c1ccccc1Cl") == "2-chlorobenzoic acid"


# ─── ベンズアルデヒド ────────────────────────────────────────────────

class TestBenzaldehyde:

    def test_benzaldehyde(self):
        assert smiles_to_iupac("O=Cc1ccccc1") == "benzaldehyde"

    def test_4_methylbenzaldehyde(self):
        assert smiles_to_iupac("O=Cc1ccc(C)cc1") == "4-methylbenzaldehyde"

    def test_4_chlorobenzaldehyde(self):
        assert smiles_to_iupac("O=Cc1ccc(Cl)cc1") == "4-chlorobenzaldehyde"
