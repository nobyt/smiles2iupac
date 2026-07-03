"""Phase 798: pyrrolo[1,2-a]pyrimidine and [1,2,4]triazolo[1,5-a]pyridine α-ol/thiol → tautomers.

- pyrrolo[1,2-a]pyrimidine C2 → 2(1H)-one/thione; C4 → 4(3H)-one/thione; C6 → 6(4H)-one/thione
- [1,2,4]triazolo[1,5-a]pyridine C2 → 2(1H)-one/thione; C5 → 5(1H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrrolo[1,2-a]pyrimidine C2-OH/SH
    ("Oc1ccn2cccc2n1",   "pyrrolo[1,2-a]pyrimidin-2(1H)-one"),
    ("Sc1ccn2cccc2n1",   "pyrrolo[1,2-a]pyrimidin-2(1H)-thione"),
    # pyrrolo[1,2-a]pyrimidine C4-OH/SH
    ("Oc1ccnc2cccn12",   "pyrrolo[1,2-a]pyrimidin-4(3H)-one"),
    ("Sc1ccnc2cccn12",   "pyrrolo[1,2-a]pyrimidin-4(3H)-thione"),
    # pyrrolo[1,2-a]pyrimidine C6-OH/SH
    ("Oc1ccc2ncccn12",   "pyrrolo[1,2-a]pyrimidin-6(4H)-one"),
    ("Sc1ccc2ncccn12",   "pyrrolo[1,2-a]pyrimidin-6(4H)-thione"),
    # [1,2,4]triazolo[1,5-a]pyridine C2-OH/SH
    ("Oc1nc2ccccn2n1",   "[1,2,4]triazolo[1,5-a]pyridin-2(1H)-one"),
    ("Sc1nc2ccccn2n1",   "[1,2,4]triazolo[1,5-a]pyridin-2(1H)-thione"),
    # [1,2,4]triazolo[1,5-a]pyridine C5-OH/SH
    ("Oc1cccc2ncnn12",   "[1,2,4]triazolo[1,5-a]pyridin-5(1H)-one"),
    ("Sc1cccc2ncnn12",   "[1,2,4]triazolo[1,5-a]pyridin-5(1H)-thione"),
    # Regressions: parent rings unchanged
    ("c1ccc2ncccn12",    "pyrrolo[1,2-a]pyrimidine"),
    ("c1cccc2ncnn12",    "[1,2,4]triazolo[1,5-a]pyridine"),
])
def test_phase798_pyrrolo_triazolo_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
