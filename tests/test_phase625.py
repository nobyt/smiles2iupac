"""Phase 625: methyl-substituted tricyclic PAH — phenanthrene, anthracene, fluorene (IUPAC 2013).

Covers all unique substitution positions.  Phenanthrene positions 9 and 10
(bay region) must receive locants 9/10, not 1/2 — this was previously broken.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phenanthrene (positions 1–4, 9 are unique; 5–8 ≡ 1–4 and 10 ≡ 9 by C2 symmetry)
    ("Cc1cccc2c1ccc1ccccc12",  "1-methylphenanthrene"),
    ("Cc1ccc2c(ccc3ccccc32)c1", "2-methylphenanthrene"),
    ("Cc1ccc2ccc3ccccc3c2c1",  "3-methylphenanthrene"),
    ("Cc1cccc2ccc3ccccc3c12",  "4-methylphenanthrene"),
    ("Cc1cc2ccccc2c2ccccc12",  "9-methylphenanthrene"),  # bay region
    # anthracene (positions 1, 2, 9 unique; 3,4 ≡ 1,2 and 10 ≡ 9 by C2h symmetry)
    ("Cc1cccc2cc3ccccc3cc12",  "1-methylanthracene"),
    ("Cc1ccc2cc3ccccc3cc2c1",  "2-methylanthracene"),
    ("Cc1c2ccccc2cc2ccccc12",  "9-methylanthracene"),   # meso position
    # fluorene (C9 is sp3; positions 1–4 and 9 unique; 5–8 ≡ 1–4 by C2v symmetry)
    ("Cc1cccc2c1Cc1ccccc1-2",  "1-methylfluorene"),
    ("Cc1ccc2c(c1)Cc1ccccc1-2", "2-methylfluorene"),
    ("Cc1ccc2c(c1)-c1ccccc1C2", "3-methylfluorene"),
    ("Cc1cccc2c1-c1ccccc1C2",  "4-methylfluorene"),
    ("CC1c2ccccc2-c2ccccc21",  "9-methylfluorene"),
])
def test_phase625_methyl_pah(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
