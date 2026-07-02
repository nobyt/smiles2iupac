"""Phase 774: 1,2,4-triazin-3(2H)-one/thione — missing C3 alpha position.

In 1,2,4-triazine (N1,N2,C3,N4,C5,C6), C3 is alpha to both N2 and N4.
Preferred tautomer puts H on N2 (lower locant).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # C3 position (between N2 and N4)
    ("Oc1nccnn1",       "1,2,4-triazin-3(2H)-one"),
    ("Sc1nccnn1",       "1,2,4-triazin-3(2H)-thione"),
    # Regression: other 1,2,4-triazine positions unchanged
    ("Oc1cncnn1",       "1,2,4-triazin-6(1H)-one"),
    ("Oc1cnncn1",       "1,2,4-triazin-5(4H)-one"),
    # Regression: parent ring
    ("c1cnncn1",        "1,2,4-triazine"),
])
def test_phase774_124triazine_c3_tautomer(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
