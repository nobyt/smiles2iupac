"""Phase 848: 10H- indicated-H retained in phenoxazine/phenothiazine derivatives.

IUPAC 2013: 10H-phenoxazine and 10H-phenothiazine have '10H-' as an intrinsic part
of the retained ring name (not a tautomeric indicator). Unlike indole where '1H-'
is dropped when N-1 is substituted (P-14.7.1), phenoxazine/phenothiazine keep
'10H-' in ALL derivatives, including those with N-10 substituents.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N-10 methyl: 10H- must be retained even though N-10 has no H
    ("CN1c2ccccc2Oc2ccccc21",        "10-methyl-10H-phenoxazine"),
    ("CN1c2ccccc2Sc2ccccc21",        "10-methyl-10H-phenothiazine"),
    # C-substituted: 10H- retained (unaffected by Phase 844/845/846)
    ("Cc1cccc2c1Nc1ccccc1O2",        "1-methyl-10H-phenoxazine"),
    ("Cc1cccc2c1Nc1ccccc1S2",        "1-methyl-10H-phenothiazine"),
    # Unsubstituted parents: 10H- present
    ("c1ccc2c(c1)Nc1ccccc1O2",       "10H-phenoxazine"),
    ("c1ccc2c(c1)Nc1ccccc1S2",       "10H-phenothiazine"),
])
def test_phase848(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
