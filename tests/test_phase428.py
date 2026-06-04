"""Phase 428: Phenanthroline retained names (IUPAC 2013 P-31.1.3).

1,10-phenanthroline, 4,7-phenanthroline, and 1,7-phenanthroline are retained
names for the three tricyclic C12H8N2 diaza-phenanthrenes listed in
IUPAC 2013 P-31.1.3.4.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,10-phenanthroline — both N adjacent to C10a junction (chelating isomer)
    ("c1ccc2c(c1)cnc1ncccc12",       "1,10-phenanthroline"),
    # 4,7-phenanthroline — N at C4 (Ring A, adj to C4a) and C7 (Ring C)
    ("c1cnc2c(c1)ccc1cnccc12",       "4,7-phenanthroline"),
    # 1,7-phenanthroline — N at C1 (Ring A) and C7 (Ring C)
    ("c1cnc2ccc3cnccc3c2c1",         "1,7-phenanthroline"),
    # regression: phenanthridine unchanged (Phase 134)
    ("c1ccc2c(c1)cnc1ccccc12",       "phenanthridine"),
    # regression: acridine unchanged
    ("c1ccc2nc3ccccc3cc2c1",          "acridine"),
    # regression: 1,5-naphthyridine unchanged (Phase 132)
    ("c1cnc2cccnc2c1",               "1,5-naphthyridine"),
    # regression: benzene unchanged
    ("c1ccccc1",                       "benzene"),
])
def test_phase428_phenanthroline(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
