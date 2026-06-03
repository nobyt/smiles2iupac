"""Phase 270: benzene tri/tetra/pentacarboxylic acids (IUPAC 2013 P-65.1.2.3).

  OC(=O)c1cccc(C(=O)O)c1C(=O)O   → benzene-1,2,3-tricarboxylic acid
  OC(=O)c1ccc(C(=O)O)c(C(=O)O)c1 → benzene-1,2,4-tricarboxylic acid
  OC(=O)c1cc(C(=O)O)cc(C(=O)O)c1 → benzene-1,3,5-tricarboxylic acid
  ...etc.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # tricarboxylic isomers
    ("OC(=O)c1cccc(C(=O)O)c1C(=O)O",       "benzene-1,2,3-tricarboxylic acid"),
    ("OC(=O)c1ccc(C(=O)O)c(C(=O)O)c1",     "benzene-1,2,4-tricarboxylic acid"),
    ("OC(=O)c1cc(C(=O)O)cc(C(=O)O)c1",     "benzene-1,3,5-tricarboxylic acid"),
    # tetracarboxylic isomers
    ("OC(=O)c1ccc(C(=O)O)c(C(=O)O)c1C(=O)O", "benzene-1,2,3,4-tetracarboxylic acid"),
    ("OC(=O)c1cc(C(=O)O)c(C(=O)O)c(C(=O)O)c1", "benzene-1,2,3,5-tetracarboxylic acid"),
    ("OC(=O)c1cc(C(=O)O)c(C(=O)O)cc1C(=O)O",   "benzene-1,2,4,5-tetracarboxylic acid"),
    # pentacarboxylic
    ("OC(=O)c1c(C(=O)O)c(C(=O)O)c(C(=O)O)c(C(=O)O)c1", "benzene-1,2,3,4,5-pentacarboxylic acid"),
    # hexacarboxylic (mellitic acid, retained name)
    ("OC(=O)c1c(C(=O)O)c(C(=O)O)c(C(=O)O)c(C(=O)O)c1C(=O)O", "mellitic acid"),
    # regressions: dicarboxylic retained names unchanged
    ("OC(=O)c1ccccc1C(=O)O",     "phthalic acid"),
    ("OC(=O)c1cccc(C(=O)O)c1",   "isophthalic acid"),
    ("OC(=O)c1ccc(C(=O)O)cc1",   "terephthalic acid"),
    # regression: simple benzoic acid unchanged
    ("OC(=O)c1ccccc1",            "benzoic acid"),
])
def test_phase270_benzene_polycarboxylic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
