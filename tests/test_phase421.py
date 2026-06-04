"""Phase 421: Acenaphthylene and acenaphthene.

IUPAC 2013 P-31.1.3: retained names for the tricyclic 5+6+6 fused
hydrocarbon ring systems.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # acenaphthylene — C12H8, tricyclic (5+6+6), C=C bridge
    ("C1=Cc2cccc3cccc1c23",     "acenaphthylene"),
    # acenaphthene — C12H10, same skeleton with CH2-CH2 bridge
    ("c1cc2c3c(cccc3c1)CC2",    "acenaphthene"),
    # regression: naphthalene unchanged
    ("c1ccc2ccccc2c1",           "naphthalene"),
    # regression: fluorene unchanged (Phase 134)
    ("c1ccc2c(c1)Cc1ccccc1-2",  "fluorene"),
    # regression: benzene unchanged
    ("c1ccccc1",                  "benzene"),
    # regression: indane unchanged (Phase 133)
    ("c1ccc2c(c1)CCC2",          "indane"),
])
def test_phase421_acenaphthylene_acenaphthene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
