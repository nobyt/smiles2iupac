"""Phase 520: secondary / tertiary amino substituent naming (IUPAC 2013 P-62.2.3)

-NHR  → {R}amino
-NR₂  → di{R}amino  (same R) or ({R1},{R2})amino (different R)
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # secondary amino substituent on a chain
    ("CNCC(=O)O",       "2-(methylamino)acetic acid"),
    ("CCNCC(=O)O",      "2-(ethylamino)acetic acid"),
    # tertiary amino substituent (same groups)
    ("CN(C)CC(=O)O",    "2-(dimethylamino)acetic acid"),
    ("CCN(CC)CC(=O)O",  "2-(diethylamino)acetic acid"),
    # tertiary amino substituent (different groups) — square brackets b/c substituent contains comma
    ("CN(CC)CC(=O)O",   "2-[(ethyl,methyl)amino]acetic acid"),
    # regression: primary amino still works (NCC(=O)O is retained "glycine")
    ("NCCC(=O)O",       "3-aminopropanoic acid"),
])
def test_phase520_amino_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
