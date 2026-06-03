"""
Phase 2 テスト: カルボン酸、アルデヒド、ケトン
"""

import pytest
from smiles2iupac import smiles_to_iupac


class TestCarboxylicAcids:

    def test_ethanoic_acid(self):
        # 酢酸: IUPAC 2013 保留名 (PIN): acetic acid
        assert smiles_to_iupac("CC(=O)O") == "acetic acid"

    def test_propanoic_acid(self):
        assert smiles_to_iupac("CCC(=O)O") == "propanoic acid"

    def test_butanoic_acid(self):
        assert smiles_to_iupac("CCCC(=O)O") == "butanoic acid"

    def test_pentanoic_acid(self):
        assert smiles_to_iupac("CCCCC(=O)O") == "pentanoic acid"


class TestAldehydes:

    def test_ethanal(self):
        # IUPAC 2013 P-31.1.3.4: acetaldehyde は保留 PIN
        assert smiles_to_iupac("CC=O") == "acetaldehyde"

    def test_propanal(self):
        assert smiles_to_iupac("CCC=O") == "propanal"

    def test_butanal(self):
        assert smiles_to_iupac("CCCC=O") == "butanal"


class TestKetones:

    def test_propan_2_one(self):
        # IUPAC 2013 P-31.1.3: acetone は保留 PIN
        assert smiles_to_iupac("CC(=O)C") == "acetone"

    def test_butan_2_one(self):
        assert smiles_to_iupac("CC(=O)CC") == "butan-2-one"

    def test_pentan_3_one(self):
        assert smiles_to_iupac("CCC(=O)CC") == "pentan-3-one"
