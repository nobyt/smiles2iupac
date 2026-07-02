"""Phase 740: drop indicated-H when N-H position is N-substituted in
2,3-dihydro-1H-indole and 2,3-dihydro-1H-isoindole (IUPAC 2013).

When the N at the indicated-H position bears a substituent the nH prefix
is dropped: 2,3-dihydro-1H-indole → 2,3-dihydroindole,
2,3-dihydro-1H-isoindole → 1,3-dihydroisoindole.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N1-substituted indoline: drop 1H
    ("CN1CCc2ccccc21",                "1-methyl-2,3-dihydroindole"),
    ("CCN1CCc2ccccc21",               "1-ethyl-2,3-dihydroindole"),
    ("Cc1ccccc1N1CCc2ccccc21",        "1-(2-methylphenyl)-2,3-dihydroindole"),
    # C-substituted indoline: keep 1H
    ("Cc1cccc2c1CCN2",                "4-methyl-2,3-dihydro-1H-indole"),
    # N2-substituted isoindoline: drop 1H
    ("CN1Cc2ccccc2C1",                "2-methyl-1,3-dihydroisoindole"),
    # C-substituted isoindoline: keep 1H
    ("Cc1ccc2c(c1)CNC2",              "5-methyl-2,3-dihydro-1H-isoindole"),
    # Regression: unsubstituted parents unchanged
    ("c1ccc2c(c1)CCN2",               "2,3-dihydro-1H-indole"),
    ("c1ccc2c(c1)CNC2",               "2,3-dihydro-1H-isoindole"),
])
def test_phase740_dihydroindole_indicated_h(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
