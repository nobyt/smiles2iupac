"""Phase 409: 3,4-Dihydroquinolin-2(1H)-one and 3,4-dihydroisoquinolin-1(2H)-one.

IUPAC 2013 P-31.1.3: retained names for partially saturated benzo-fused lactams.
Preferred over the 'N-oxotetrahydro' substituent-based names.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 3,4-dihydroquinolin-2(1H)-one (dihydroquinolinone)
    ("O=C1CCc2ccccc2N1",               "3,4-dihydroquinolin-2(1H)-one"),
    # N-methyl-3,4-dihydroquinolin-2(1H)-one (N1 adjacent to C8a junction)
    ("CN1C(=O)CCc2ccccc21",            "1-methyl-3,4-dihydroquinolin-2(1H)-one"),
    # N-methyl-3,4-dihydroisoquinolin-1(2H)-one (N2 between C1=O and C3)
    ("CN1CCc2ccccc2C1=O",              "2-methyl-3,4-dihydroisoquinolin-1(2H)-one"),
    # 3,4-dihydroisoquinolin-1(2H)-one
    ("O=C1NCCc2ccccc21",               "3,4-dihydroisoquinolin-1(2H)-one"),
    # regression: indolin-2-one (5-membered) unchanged
    ("O=C1Cc2ccccc2N1",                "indolin-2-one"),
    # regression: isoindolin-1-one unchanged
    ("O=C1CNc2ccccc21",                "isoindolin-1-one"),
    # regression: 1,2,3,4-tetrahydroquinoline unchanged
    ("c1ccc2c(c1)CCCN2",               "1,2,3,4-tetrahydroquinoline"),
    # regression: benzene unchanged
    ("c1ccccc1",                        "benzene"),
])
def test_phase409_dihydroquinolinones(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
