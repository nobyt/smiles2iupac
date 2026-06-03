"""Phase 253: partial phosphonate ester naming (IUPAC 2013 P-62.5.2).

R-P(=O)(OR')(OH) → R' hydrogen Rphosphonate
Previously mis-classified as Rphosphinic acid.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # partial phosphonate esters
    ("CP(=O)(OC)O",    "methyl hydrogen methylphosphonate"),
    ("CCP(=O)(OC)O",   "methyl hydrogen ethylphosphonate"),
    ("CP(=O)(OCC)O",   "ethyl hydrogen methylphosphonate"),
    ("CCP(=O)(OCC)O",  "ethyl hydrogen ethylphosphonate"),
    # regression: full phosphonate ester unchanged
    ("CP(=O)(OC)OC",   "dimethyl methylphosphonate"),
    ("CCP(=O)(OC)OC",  "dimethyl ethylphosphonate"),
    # regression: phosphonic acid unchanged
    ("CP(=O)(O)O",     "methylphosphonic acid"),
    ("CCP(=O)(O)O",    "ethylphosphonic acid"),
    # regression: phosphinic acid (mono, no ester O) unchanged
    ("C[PH](=O)O",     "methylphosphinic acid"),
])
def test_phase253_phosphonate_halfester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
