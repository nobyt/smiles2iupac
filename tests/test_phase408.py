"""Phase 408: Fluorenone and thioxanthen-9-one retained names.

IUPAC 2013 P-31.1.3: 9H-fluoren-9-one (fluorenone) and thioxanthen-9-one
are retained names for tricyclic ketones derived from fluorene and thioxanthene.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 9H-fluoren-9-one (fluorenone)
    ("O=C1c2ccccc2-c2ccccc21",        "9H-fluoren-9-one"),
    # thioxanthen-9-one (thioxanthone)
    ("O=C1c2ccccc2Sc2ccccc21",        "thioxanthen-9-one"),
    # regression: fluorene unchanged
    ("c1ccc2c(c1)Cc1ccccc1-2",        "fluorene"),
    # regression: xanthen-9-one unchanged
    ("O=C1c2ccccc2Oc2ccccc21",        "xanthen-9-one"),
    # regression: thioxanthene unchanged
    ("c1ccc2c(c1)Cc1ccccc1S2",        "thioxanthene"),
    # regression: benzene unchanged
    ("c1ccccc1",                       "benzene"),
])
def test_phase408_tricyclic_ketones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
