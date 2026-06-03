"""Phase 359: Hydroxysulfoxide/hydroxysulfone naming (IUPAC 2013).

In IUPAC 2013, sulfoxide and sulfone are NOT principal characteristic groups.
When they co-occur with alcohol (or other higher-priority groups), alcohol wins
and the sulfinyl/sulfonyl group is expressed as a substituent prefix.

Previously CC(O)S(=O)C gave "(methylsulfinyl)1-hydroxyethane" because sulfoxide
(priority 62) beat alcohol (priority 50).  After lowering sulfoxide/sulfone
priorities, alcohol wins and the S=O group is named "(methylsulfinyl)".
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Hydroxysulfoxide: alcohol beats sulfoxide
    ("CC(O)S(=O)C",            "1-(methylsulfinyl)ethanol"),
    ("OCCS(=O)C",              "2-(methylsulfinyl)ethanol"),
    ("OCCCS(=O)C",             "3-(methylsulfinyl)propan-1-ol"),
    # Hydroxysulfone: alcohol beats sulfone
    ("CC(O)S(=O)(=O)C",        "1-(methylsulfonyl)ethanol"),
    ("OCCS(=O)(=O)C",          "2-(methylsulfonyl)ethanol"),
    # Hydroxysulfide: already worked (sulfide priority 30 < alcohol 50)
    ("CC(O)SC",                "1-(methylsulfanyl)ethanol"),
    ("OCCSC",                  "2-(methylsulfanyl)ethanol"),
    # Regressions: pure sulfoxide/sulfone still named correctly
    ("CS(=O)C",                "dimethyl sulfoxide"),
    ("CCS(=O)C",               "ethyl methyl sulfoxide"),
    ("CS(=O)(=O)C",            "dimethyl sulfone"),
    ("CCS(=O)(=O)C",           "ethyl methyl sulfone"),
    # Regression: vinyl sulfoxide/sulfone (Phase 357)
    ("C=CS(=O)C",              "ethenyl methyl sulfoxide"),
    ("C=CS(=O)(=O)C",          "ethenyl methyl sulfone"),
])
def test_phase359_hydroxysulfoxide_sulfone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
