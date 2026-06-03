"""Phase 315: benzenecarboximidamide and related exocyclic amidine naming (IUPAC 2013 P-66.6.1).

Amidine (imidamide) group attached to benzene ring → "benzenecarboximidamide".
Also adds nitrile/thioamide/selenoamide to polycyclic exo-suffix map.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzene
    ("NC(=N)c1ccccc1",         "benzenecarboximidamide"),
    ("NC(=N)c1ccc(Cl)cc1",    "4-chlorobenzenecarboximidamide"),
    ("NC(=N)c1ccccc1Cl",      "2-chlorobenzenecarboximidamide"),
    # cycloalkane
    ("NC(=N)C1CCCCC1",        "cyclohexanecarboximidamide"),
    ("NC(=N)C1CCCC1",         "cyclopentanecarboximidamide"),
    # regressions: chain amidines unchanged
    ("NC(=N)C",               "ethanimidamide"),
    ("NC(=N)CC",              "propanimidamide"),
    # regression: benzamide unchanged
    ("NC(=O)c1ccccc1",        "benzamide"),
])
def test_phase315_benzene_carboximidamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
