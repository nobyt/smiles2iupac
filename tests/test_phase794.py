"""Phase 794: furo[3,2-b]pyridine, [1,2,5]thiadiazolo/oxadiazolo[3,4-b]pyridine and [3,4-e]pyrazine α-ol/thiol → tautomers.

- furo[3,2-b]pyridine C5 → 5(4H)-one/thione
- [1,2,5]thiadiazolo[3,4-b]pyridine C5 → 5(4H)-one/thione
- [1,2,5]oxadiazolo[3,4-b]pyridine C5 → 5(4H)-one/thione
- [1,2,5]thiadiazolo[3,4-e]pyrazine C5 → 5(4H)-one/thione
- [1,2,5]oxadiazolo[3,4-e]pyrazine C5 → 5(4H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # furo[3,2-b]pyridine C5-OH/SH
    ("Oc1ccc2occc2n1",   "furo[3,2-b]pyridin-5(4H)-one"),
    ("Sc1ccc2occc2n1",   "furo[3,2-b]pyridin-5(4H)-thione"),
    # [1,2,5]thiadiazolo[3,4-b]pyridine C5-OH/SH
    ("Oc1ccc2nsnc2n1",   "[1,2,5]thiadiazolo[3,4-b]pyridin-5(4H)-one"),
    ("Sc1ccc2nsnc2n1",   "[1,2,5]thiadiazolo[3,4-b]pyridin-5(4H)-thione"),
    # [1,2,5]oxadiazolo[3,4-b]pyridine C5-OH/SH
    ("Oc1ccc2nonc2n1",   "[1,2,5]oxadiazolo[3,4-b]pyridin-5(4H)-one"),
    ("Sc1ccc2nonc2n1",   "[1,2,5]oxadiazolo[3,4-b]pyridin-5(4H)-thione"),
    # [1,2,5]thiadiazolo[3,4-e]pyrazine C5-OH/SH
    ("Oc1cnc2nsnc2n1",   "[1,2,5]thiadiazolo[3,4-e]pyrazin-5(4H)-one"),
    ("Sc1cnc2nsnc2n1",   "[1,2,5]thiadiazolo[3,4-e]pyrazin-5(4H)-thione"),
    # [1,2,5]oxadiazolo[3,4-e]pyrazine C5-OH/SH
    ("Oc1cnc2nonc2n1",   "[1,2,5]oxadiazolo[3,4-e]pyrazin-5(4H)-one"),
    ("Sc1cnc2nonc2n1",   "[1,2,5]oxadiazolo[3,4-e]pyrazin-5(4H)-thione"),
    # Regressions: parent rings unchanged
    ("c1ccc2occc2n1",    "furo[3,2-b]pyridine"),
    ("c1ccc2nsnc2n1",    "[1,2,5]thiadiazolo[3,4-b]pyridine"),
    ("c1ccc2nonc2n1",    "[1,2,5]oxadiazolo[3,4-b]pyridine"),
    ("c1cnc2nsnc2n1",    "[1,2,5]thiadiazolo[3,4-e]pyrazine"),
    ("c1cnc2nonc2n1",    "[1,2,5]oxadiazolo[3,4-e]pyrazine"),
])
def test_phase794_fused_bicyclic_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
