"""Phase 250: carbonate half-ester naming (IUPAC 2013 P-65.1).

RO-C(=O)-OH → alkyl hydrogen carbonate
Distinguish from dialkyl carbonate (RO-C(=O)-OR) and methyl formate (H-C(=O)-OR).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # carbonate half-esters
    ("COC(=O)O",      "methyl hydrogen carbonate"),
    ("CCOC(=O)O",     "ethyl hydrogen carbonate"),
    ("CCCOC(=O)O",    "propyl hydrogen carbonate"),
    # regression: dialkyl carbonate unchanged
    ("COC(=O)OC",     "dimethyl carbonate"),
    ("COC(=O)OCC",    "ethyl methyl carbonate"),
    # regression: formate (H directly on carbonyl C) unchanged
    ("O=COC",         "methyl formate"),
    ("O=COCC",        "ethyl formate"),
    # regression: carbonic acid unchanged
    ("OC(=O)O",       "carbonic acid"),
])
def test_phase250_carbonate_halfester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
