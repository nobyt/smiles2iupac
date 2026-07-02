"""Phase 773: pteridine α-ol/thiol → tautomers (IUPAC 2013).

Pteridine (N1, N3, N5, N8) has four non-junction alpha C positions:
- C2 (alpha to N1 and N3) → 2(1H)-one  [N1-H preferred]
- C4 (alpha to N3) → 4(3H)-one
- C6 (alpha to N5) → 6(5H)-one
- C7 (alpha to N8) → 7(8H)-one
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # C2
    ("Oc1ncc2nccnc2n1",     "pteridin-2(1H)-one"),
    ("Sc1ncc2nccnc2n1",     "pteridin-2(1H)-thione"),
    # C4
    ("Oc1ncnc2nccnc12",     "pteridin-4(3H)-one"),
    ("Sc1ncnc2nccnc12",     "pteridin-4(3H)-thione"),
    # C6
    ("Oc1cnc2ncncc2n1",     "pteridin-6(5H)-one"),
    ("Sc1cnc2ncncc2n1",     "pteridin-6(5H)-thione"),
    # C7
    ("Oc1cnc2cncnc2n1",     "pteridin-7(8H)-one"),
    ("Sc1cnc2cncnc2n1",     "pteridin-7(8H)-thione"),
    # Regression: parent ring unaffected
    ("c1cnc2ncncc2n1",      "pteridine"),
    # Regression: Phase 752 purine unchanged
    ("Oc1ncnc2[nH]cnc12",   "9H-purin-6(1H)-one"),
])
def test_phase773_pteridine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
