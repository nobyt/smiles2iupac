"""Phase 568: Substituted fluoranthene naming.
Fluoranthene: C2v-symmetric (16 atoms), 5 unique C environments:
{1,6}, {2,5}, {3,4}, {7,10}, {8,9}. IUPAC uses minimum locant.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("c1ccc2c(c1)-c1cccc3cccc-2c13",       "fluoranthene"),
    ("Cc1ccc2cccc3c2c1-c1ccccc1-3",        "1-methylfluoranthene"),
    ("Cc1cc2c3c(cccc3c1)-c1ccccc1-2",      "2-methylfluoranthene"),
    ("Cc1ccc2c3c(cccc13)-c1ccccc1-2",      "3-methylfluoranthene"),
    ("Cc1cccc2c1-c1cccc3cccc-2c13",        "7-methylfluoranthene"),
    ("Cc1ccc2c(c1)-c1cccc3cccc-2c13",      "8-methylfluoranthene"),
])
def test_phase568_fluoranthene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
