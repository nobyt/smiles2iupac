"""Phase 807: pyrrolo[1,2-b]pyridazine α-ol/thiol → tautomers.

- C2 → 2(1H)-one/thione; C7 → 7(1H)-one/thione
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrrolo[1,2-b]pyridazine C2-OH/SH
    ("Oc1ccc2cccn2n1",   "pyrrolo[1,2-b]pyridazin-2(1H)-one"),
    ("Sc1ccc2cccn2n1",   "pyrrolo[1,2-b]pyridazin-2(1H)-thione"),
    # pyrrolo[1,2-b]pyridazine C7-OH/SH
    ("Oc1ccc2cccnn12",   "pyrrolo[1,2-b]pyridazin-7(1H)-one"),
    ("Sc1ccc2cccnn12",   "pyrrolo[1,2-b]pyridazin-7(1H)-thione"),
    # Regression: parent ring unchanged
    ("c1cnn2cccc2c1",    "pyrrolo[1,2-b]pyridazine"),
])
def test_phase807_pyrrolo_pyridazine_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
