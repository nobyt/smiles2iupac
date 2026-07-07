"""Phase 752: pyrimidine-2,4-diol (enol) → pyrimidine-2,4(1H,3H)-dione (IUPAC 2013).

The enol (di-OH) form of uracil/thymine prefers the di-lactam (dione) tautomer.
Extends the mono-OH → mono-one rule (Phase 744) to the pyrimidine-2,4 diol case.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # uracil enol (di-OH) → dione
    ("Oc1ccnc(O)n1",              "pyrimidine-2,4(1H,3H)-dione"),
    ("Oc1nc(O)ccn1",              "pyrimidine-2,4(1H,3H)-dione"),
    # thymine enol
    ("Oc1c(C)cnc(O)n1",           "5-methylpyrimidine-2,4(1H,3H)-dione"),
    # Regression: keto (dione) SMILES unchanged
    ("O=c1cc[nH]c(=O)[nH]1",      "pyrimidine-2,4(1H,3H)-dione"),
    ("Cc1c[nH]c(=O)[nH]c1=O",     "5-methylpyrimidine-2,4(1H,3H)-dione"),
    # Regression: mono-OH tautomers unchanged
    ("Oc1ccncn1",                  "1H-pyrimidin-6-one"),
    ("Oc1ncccn1",                  "1H-pyrimidin-2-one"),
    # Regression: plain pyrimidine unchanged
    ("c1ccncn1",                   "pyrimidine"),
])
def test_phase752_pyrimidine_diol_tautomer(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
