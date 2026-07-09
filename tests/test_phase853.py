"""Phase 853: N-methyl benzoxazol/benzothiazol-2-thione → -2-thione suffix.

Root cause: PGRP_DISPATCH at lines 512-518 of __init__.py was firing before
the early fused-hetero-retained check (line 621), misidentifying the ring C=S
adjacent to ring O (benzoxazolethione) or ring S (benzothiazolethione) as a
thioester/dithiocarboxylate functional group.

Fix (Phase 853): insert a _try_fused_hetero_retained guard immediately before
PGRP_DISPATCH, gated on the pgrp anchor being in a ring. This makes the fused
heterocycle path win when it can produce a name.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N-methyl benzoxazol-2-thione (C=S adjacent to ring O — was mis-named as thioester)
    ("Cn1c(=S)oc2ccccc21",          "3-methyl-1,3-benzoxazol-2-thione"),
    # N-methyl benzothiazol-2-thione (C=S adjacent to ring S — was mis-named as dithio)
    ("Cn1c(=S)sc2ccccc21",          "3-methyl-1,3-benzothiazol-2-thione"),
    # Parent heterocycles unaffected
    ("c1ccc2ocnc2c1",               "1,3-benzoxazole"),
    ("c1ccc2scnc2c1",               "1,3-benzothiazole"),
    # C-substituent on ring (not N-sub) is unaffected
    ("Cc1nc2ccccc2o1",              "2-methyl-1,3-benzoxazole"),
    ("Cc1nc2ccccc2s1",              "2-methyl-1,3-benzothiazole"),
])
def test_phase853(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
