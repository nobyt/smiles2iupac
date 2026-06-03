"""Phase 230: styrene retained name (IUPAC 2013 P-31.1.3.4).

C=Cc1ccccc1 → "styrene" is a preferred IUPAC name (ethenylbenzene).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C=Cc1ccccc1", "styrene"),
    # regression: other vinylated benzenes use systematic names
    ("CC=Cc1ccccc1", "(prop-1-en-1-yl)benzene"),
    # regression: toluene still works
    ("Cc1ccccc1", "toluene"),
])
def test_phase230_styrene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
