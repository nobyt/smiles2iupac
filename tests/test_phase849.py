"""Phase 849: N-methyl heterocyclic lactams → -one suffix (not oxo- prefix).

IUPAC 2013: quinolinone, isoquinolinone, acridinone, phenanthridinone retain
the -one suffix even when the ring N is substituted. When RDKit aromatizes
these N-methyl keto forms, the C=O appears as sub_nm=="oxo"; Phase 849 converts
this to the -one suffix. Acridine locant map: N-10 now assigned locant 10.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Quinoline -ones: N-methyl (oxo→ -one suffix, no indicated-H)
    ("Cn1ccc(=O)c2ccccc21",          "1-methylquinolin-4-one"),
    ("Cn1c(=O)ccc2ccccc21",          "1-methylquinolin-2-one"),
    # Quinoline -ones: unsubstituted NH (indicated-H required)
    ("O=c1ccc2ccccc2[nH]1",          "quinolin-2(1H)-one"),
    ("O=c1cc[nH]c2ccccc12",          "quinolin-4(1H)-one"),
    # Isoquinoline -ones: N-methyl
    ("Cn1ccc2ccccc2c1=O",            "2-methylisoquinolin-1-one"),
    ("Cn1cc2ccccc2cc1=O",            "2-methylisoquinolin-3-one"),
    # Isoquinoline -ones: unsubstituted NH
    ("O=c1[nH]ccc2ccccc12",          "isoquinolin-1(2H)-one"),
    ("O=c1cc2ccccc2c[nH]1",          "isoquinolin-3(2H)-one"),
    # Acridine -ones: N-methyl (requires acridine locant map N10=10)
    ("Cn1c2ccccc2c(=O)c2ccccc21",    "10-methylacridin-9-one"),
    # Acridine -ones: unsubstituted NH
    ("O=c1c2ccccc2[nH]c2ccccc12",    "acridin-9(10H)-one"),
    # Phenanthridine -ones: N-methyl
    ("Cn1c(=O)c2ccccc2c2ccccc21",    "5-methylphenanthridin-6-one"),
    # Phenanthridine -ones: unsubstituted NH
    ("O=c1[nH]c2ccccc2c2ccccc12",    "phenanthridin-6(5H)-one"),
])
def test_phase849(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
