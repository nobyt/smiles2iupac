"""Phase 411: Chromone (4H-chromen-4-one) and indan-2-one retained names.

IUPAC 2013 P-31.1.3: chromone is a retained name for 4H-chromen-4-one;
indan-2-one is systematic for the 5-membered ring ketone fused at C2.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # chromone (4H-chromen-4-one) — aromatic pyranone
    ("O=c1ccoc2ccccc12",               "chromone"),
    # indan-2-one (5-membered ring with ketone at C2)
    ("O=C1Cc2ccccc2C1",                "indan-2-one"),
    # regression: coumarin unchanged
    ("O=c1ccc2ccccc2o1",               "coumarin"),
    # regression: isocoumarin unchanged
    ("O=c1occc2ccccc12",               "isocoumarin"),
    # regression: indan-1-one unchanged
    ("O=C1CCc2ccccc21",                "indan-1-one"),
    # regression: chromane unchanged
    ("C1CCc2ccccc2O1",                 "chromane"),
    # regression: chroman-4-one unchanged
    ("O=C1CCOc2ccccc21",               "chroman-4-one"),
    # regression: benzene unchanged
    ("c1ccccc1",                        "benzene"),
])
def test_phase411_chromone_indan2one(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
