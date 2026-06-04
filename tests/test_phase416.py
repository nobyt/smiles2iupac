"""Phase 416: Quinolin-4(1H)-one, isoquinolin-1(2H)-one, isoquinolin-3(2H)-one.

IUPAC 2013 P-31.1.3: retained/systematic names for benzo-fused pyridinone
lactam ring systems derived from quinoline and isoquinoline.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # quinolin-4(1H)-one (C4=O, N1-H)
    ("O=c1cc[nH]c2ccccc12",            "quinolin-4(1H)-one"),
    # isoquinolin-1(2H)-one (C1=O, N2-H)
    ("O=c1[nH]ccc2ccccc12",            "isoquinolin-1(2H)-one"),
    # isoquinolin-3(2H)-one (C3=O, N2-H)
    ("O=c1cc2ccccc2c[nH]1",            "isoquinolin-3(2H)-one"),
    # regression: quinolin-2(1H)-one unchanged (Phase 415)
    ("O=c1ccc2ccccc2[nH]1",            "quinolin-2(1H)-one"),
    # regression: quinoline unchanged
    ("c1ccc2ncccc2c1",                  "quinoline"),
    # regression: isoquinoline unchanged
    ("c1ccc2cnccc2c1",                  "isoquinoline"),
    # regression: quinoxalin-2(1H)-one unchanged (Phase 413)
    ("O=c1cnc2ccccc2[nH]1",            "quinoxalin-2(1H)-one"),
    # regression: benzene unchanged
    ("c1ccccc1",                         "benzene"),
])
def test_phase416_quinolinones_isoquinolinones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
