"""Phase 312: cycloalkane-1-carboxylic acid locant citation with substituents (IUPAC 2013 P-44.1).

When substituents are present on a cycloalkane, the C1 locant for exocyclic principal
characteristic groups (carboxylic acid, carbaldehyde, carboxamide, carbonitrile) must be cited.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("NC1CCCCC1C(=O)O",  "2-aminocyclohexane-1-carboxylic acid"),
    ("NC1CCCCC1C=O",     "2-aminocyclohexane-1-carbaldehyde"),
    ("NC1CCCCC1C(=O)N",  "2-aminocyclohexane-1-carboxamide"),
    ("NC1CCCCC1C#N",     "2-aminocyclohexane-1-carbonitrile"),
    # regressions: unsubstituted rings suppress C1 locant
    ("OC(=O)C1CCCCC1",  "cyclohexanecarboxylic acid"),
    ("O=CC1CCCCC1",     "cyclohexanecarbaldehyde"),
    ("NC(=O)C1CCCCC1",  "cyclohexanecarboxamide"),
    ("N#CC1CCCCC1",     "cyclohexanecarbonitrile"),
])
def test_phase312_exocyclic_locant_with_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
