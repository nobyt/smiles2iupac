"""Phase 842: pyrimidin-4(3H)-one/thione and pyridazin-3(2H)-one/thione → IUPAC 2013 PINs.

IUPAC 2013 P-31.1.2: indicated H must be at the lowest possible locant.
Locant set comparison: {1,6} < {3,4} for pyrimidine; {1,6} < {2,3} for pyridazine.
Substituent locants are remapped accordingly (pyrimidine: C4↔C6; pyridazine: C3↔C6, C4↔C5).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # monocyclic pyrimidine lactam/thiolactam PINs
    ("O=C1NC=NC=C1",        "1H-pyrimidin-6-one"),
    ("S=C1NC=NC=C1",        "1H-pyrimidin-6-thione"),
    # monocyclic pyridazine lactam/thiolactam PINs
    ("O=C1C=CC=NN1",        "1H-pyridazin-6-one"),
    ("Sc1cccnn1",           "1H-pyridazin-6-thione"),
    # substituted: pyridazine C6-methyl flips to C3 after renumbering
    ("Cc1ccc(O)nn1",        "3-methyl-1H-pyridazin-6-one"),
    # substituted: pyrimidine C2-methyl stays at C2 (C2 doesn't flip)
    ("Oc1ccnc(C)n1",        "2-methyl-1H-pyrimidin-6-one"),
    # regression: pyrimidin-2(1H)-one unchanged (different position)
    ("O=C1NC=CC=N1",        "pyrimidin-2(1H)-one"),
    # regression: pyrazin-2(1H)-one unchanged
    ("O=C1NC=CN=C1",        "pyrazin-2(1H)-one"),
])
def test_phase842(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
