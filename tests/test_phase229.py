"""Phase 229: toluene and acetic anhydride retained names (IUPAC 2013).

IUPAC 2013 P-31.1.3.4: "toluene" is a preferred IUPAC name for methylbenzene.
IUPAC 2013 P-65.1.1.3.2: "acetic anhydride" is a preferred IUPAC name.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # toluene
    ("Cc1ccccc1",         "toluene"),
    # acetic anhydride
    ("CC(=O)OC(=O)C",     "acetic anhydride"),
    # regression: other substituted benzenes still use systematic names
    ("CCc1ccccc1",        "ethylbenzene"),
    # regression: other anhydrides still work
    ("CCC(=O)OC(=O)CC",   "propanoic anhydride"),
])
def test_phase229_toluene_acetic_anhydride(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
