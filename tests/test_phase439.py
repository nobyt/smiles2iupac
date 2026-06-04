"""Phase 439: oxazolo- and thiazolo-pyridine/pyrimidine bicyclic retained names
(IUPAC 2013 P-31.1.3).

These C5-6N2O/S bicyclics currently output 'benzene' (wrong).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # oxazolo[4,5-b]pyridine — O then N in 5-membered ring
    ("c1cnc2ocnc2c1",              "oxazolo[4,5-b]pyridine"),
    # oxazolo[5,4-b]pyridine — N then O in 5-membered ring
    ("c1cnc2ncoc2c1",              "oxazolo[5,4-b]pyridine"),
    # thiazolo[4,5-b]pyridine — S then N in 5-membered ring
    ("c1cnc2scnc2c1",              "thiazolo[4,5-b]pyridine"),
    # thiazolo[5,4-b]pyridine — N then S in 5-membered ring
    ("c1cnc2ncsc2c1",              "thiazolo[5,4-b]pyridine"),
    # oxazolo[4,5-d]pyrimidine — oxazole fused to pyrimidine
    ("c1cnc2ocnc2n1",              "oxazolo[4,5-d]pyrimidine"),
    # thiazolo[4,5-d]pyrimidine — thiazole fused to pyrimidine
    ("c1cnc2scnc2n1",              "thiazolo[4,5-d]pyrimidine"),
    # regression: 3H-imidazo[4,5-b]pyridine unchanged
    ("c1cnc2nc[nH]c2c1",           "3H-imidazo[4,5-b]pyridine"),
    # regression: furo[2,3-b]pyridine unchanged (Phase 436)
    ("c1cnc2occc2c1",              "furo[2,3-b]pyridine"),
])
def test_phase439_oxazolo_thiazolo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
