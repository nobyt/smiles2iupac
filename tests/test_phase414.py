"""Phase 414: Phthalazin-1(2H)-one, cinnolin-4(1H)-one, quinazolin-4(3H)-one.

IUPAC 2013 P-31.1.3: retained/systematic names for benzo-fused pyridazinone
and pyrimidinone lactam ring systems derived from phthalazine, cinnoline,
and quinazoline.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phthalazin-1(2H)-one (C1=O, N2-H, N3=)
    ("O=c1[nH]ncc2ccccc12",            "phthalazin-1(2H)-one"),
    # cinnolin-4(1H)-one (C4=O, N1-H, N2=)
    ("O=c1cn[nH]c2ccccc12",            "cinnolin-4(1H)-one"),
    # quinazolin-4(3H)-one (C4=O, N3-H, N1=)
    ("O=c1[nH]cnc2ccccc12",            "quinazolin-4(3H)-one"),
    # regression: phthalazine unchanged
    ("c1ccc2cnncc2c1",                  "phthalazine"),
    # regression: cinnoline unchanged
    ("c1ccc2nnccc2c1",                  "cinnoline"),
    # regression: quinazoline unchanged
    ("c1ccc2ncncc2c1",                  "quinazoline"),
    # regression: quinoxalin-2(1H)-one (already Phase 413) unchanged
    ("O=c1cnc2ccccc2[nH]1",            "quinoxalin-2(1H)-one"),
    # regression: benzene unchanged
    ("c1ccccc1",                         "benzene"),
])
def test_phase414_benzodiazinones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
