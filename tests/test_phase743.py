"""Phase 743: phthalimide → isoindole-1,3(2H)-dione (IUPAC 2013 PIN).

Retained name 'phthalimide' is not a PIN; the systematic name
isoindole-1,3(2H)-dione is used instead.
N-substituted forms drop (2H) since N2 carries a substituent, not H.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Unsubstituted: phthalimide → isoindole-1,3(2H)-dione
    ("O=C1NC(=O)c2ccccc21",             "isoindole-1,3(2H)-dione"),
    # N-substituted: drop (2H) — N2 carries substituent, not H (Phase 843)
    ("CN1C(=O)c2ccccc2C1=O",            "2-methylisoindole-1,3-dione"),
    ("CCN1C(=O)c2ccccc2C1=O",           "2-ethylisoindole-1,3-dione"),
    # Regression: succinimide stays systematic (not a retained name)
    ("O=C1CCC(=O)N1",                   "pyrrolidine-2,5-dione"),
])
def test_phase743_phthalimide_pin(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
