"""Phase 413: Acridin-9(10H)-one and quinoxalin-2(1H)-one retained names.

IUPAC 2013 P-31.1.3: acridinone and quinoxalinone are retained/systematic names
for the tricyclic and bicyclic N-containing ketones derived from acridine
and quinoxaline.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # acridin-9(10H)-one (acridinone) — C9=O, N10-H
    ("O=c1c2ccccc2[nH]c2ccccc12",     "acridin-9(10H)-one"),
    # quinoxalin-2(1H)-one — C2=O, N1-H, N4
    ("O=c1cnc2ccccc2[nH]1",            "quinoxalin-2(1H)-one"),
    # regression: acridine unchanged
    ("c1ccc2nc3ccccc3cc2c1",            "acridine"),
    # regression: quinoxaline unchanged
    ("c1cnc2ccccc2n1",                  "quinoxaline"),
    # regression: xanthen-9-one unchanged
    ("O=C1c2ccccc2Oc2ccccc21",          "xanthen-9-one"),
    # regression: chromone unchanged
    ("O=c1ccoc2ccccc12",                "chromone"),
    # regression: benzene unchanged
    ("c1ccccc1",                         "benzene"),
])
def test_phase413_acridinone_quinoxalinone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
