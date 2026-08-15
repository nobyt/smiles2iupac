"""Phase 520: secondary / tertiary amino substituent naming (IUPAC 2013 P-62.2.3)

-NHR  → {R}amino
-NR₂  → di{R}amino  (same R) or ({R1})({R2})amino (different R)

Phase 925 update: the different-R format below was originally
"({R1},{R2})amino" (comma-joined inside one paren pair). Confirmed via
OPSIN + RDKit canonical-SMILES that this comma-joined format does NOT
round-trip to the correct structure -- OPSIN silently merges the two
comma-separated substituent names into a single different alkyl group
instead of two independent N-substituents. Fixed to individually
bracket each substituent and concatenate with no separator:
"({R1})({R2})amino".
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
    # tertiary amino substituent (different groups) — square brackets since
    # the inner substituent name itself needs parens (Phase925: no comma)
    ("CN(CC)CC(=O)O",   "2-[(ethyl)(methyl)amino]acetic acid"),
    # regression: primary amino still works
    ("NCCC(=O)O",       "3-aminopropanoic acid"),
])
def test_phase520_amino_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
