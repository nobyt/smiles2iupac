"""Phase 422: Pyrene and fluoranthene.

IUPAC 2013 P-31.1.3: retained names for the tetracyclic fused PAH ring
systems (pyrene 4×6; fluoranthene 5+6+6+6).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrene — C16H10, four fused 6-membered rings
    ("c1cc2ccc3cccc4ccc(c1)c2c34",     "pyrene"),
    # fluoranthene — C16H10, one 5-membered ring + three 6-membered rings
    ("c1ccc2c(c1)-c1cccc3cccc-2c13",   "fluoranthene"),
    # regression: acenaphthylene unchanged (Phase 421)
    ("C1=Cc2cccc3cccc1c23",             "acenaphthylene"),
    # regression: naphthalene unchanged
    ("c1ccc2ccccc2c1",                   "naphthalene"),
    # regression: benzene unchanged
    ("c1ccccc1",                          "benzene"),
])
def test_phase422_pyrene_fluoranthene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
