"""Phase 410: β-Tetralone, isochroman-1-one, chroman-2-one, chroman-4-one.

IUPAC 2013 P-31.1.3: retained names for partially saturated fused ring
carbocycles and lactones derived from naphthalene, isochroman, and chromane.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # β-tetralone (3,4-dihydronaphthalen-1(2H)-one)
    ("O=C1CCCc2ccccc21",              "3,4-dihydronaphthalen-1(2H)-one"),
    # isochroman-1-one (3,4-dihydroisocoumarin)
    ("O=C1OCCc2ccccc21",              "isochroman-1-one"),
    # chroman-2-one (6-membered lactone with O adjacent to C8a)
    ("O=C1CCc2ccccc2O1",              "chroman-2-one"),
    # chroman-4-one (flavanone skeleton without substituents)
    ("O=C1CCOc2ccccc21",              "chroman-4-one"),
    # regression: indan-1-one (5-membered) unchanged
    ("O=C1CCc2ccccc21",               "indan-1-one"),
    # regression: chromane unchanged
    ("C1CCc2ccccc2O1",                "chromane"),
    # regression: isochromane unchanged
    ("C1COCc2ccccc21",                "isochromane"),
    # regression: benzene unchanged
    ("c1ccccc1",                       "benzene"),
])
def test_phase410_benzocyclicketones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
