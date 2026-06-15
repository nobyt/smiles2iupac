"""Phase 460: phenanthroline isomers 1,6 / 2,6 / 2,7 / 3,5 / 3,6 / 4,5
(IUPAC 2013 P-31.1.3.4 retained names).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1cnc2ccc3ccncc3c2c1",     "2,7-phenanthroline"),
    ("c1ccc2c(c1)ncc1ccncc12",   "2,6-phenanthroline"),
    ("c1cnc2c(c1)cnc1ccccc12",   "1,6-phenanthroline"),
    ("c1ccc2c(c1)ncc1cnccc12",   "3,6-phenanthroline"),
    ("c1ccc2c(c1)cnc1cnccc12",   "3,5-phenanthroline"),
    ("c1ccc2c(c1)cnc1ncccc12",   "4,5-phenanthroline"),
])
def test_phase460_phenanthroline_isomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
