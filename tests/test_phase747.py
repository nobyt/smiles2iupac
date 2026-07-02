"""Phase 747: 2-hydroxypyrazine, 3-hydroxypyridazine, 4-hydroxypyrimidine →
their preferred lactam tautomers (IUPAC 2013 PINs).

All α-hydroxy positions adjacent to ring N prefer the N-H keto form.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 2-hydroxypyrazine → pyrazin-2(1H)-one
    ("Oc1cnccn1",            "pyrazin-2(1H)-one"),
    # 3-hydroxypyridazine → pyridazin-3(2H)-one
    ("Oc1cccnn1",            "pyridazin-3(2H)-one"),
    # 4-hydroxypyrimidine → pyrimidin-4(3H)-one
    ("Oc1ccncn1",            "pyrimidin-4(3H)-one"),
    # Substituted cases
    ("Cc1cnc(O)cn1",         "5-methylpyrazin-2(1H)-one"),
    ("Cc1ccc(O)nn1",         "6-methylpyridazin-3(2H)-one"),
    ("Oc1ccnc(C)n1",         "2-methylpyrimidin-4(3H)-one"),
    # Regression: keto SMILES already correct
    ("O=C1NC=CN=C1",         "pyrazin-2(1H)-one"),
    ("O=C1C=CC=NN1",         "pyridazin-3(2H)-one"),
    ("O=C1NC=NC=C1",         "pyrimidin-4(3H)-one"),
    # Regression: Phase 744 pyrimidine-2-ol unaffected
    ("Oc1ncccn1",            "1H-pyrimidin-2-one"),
    # Regression: parent rings unaffected
    ("c1cnccn1",             "pyrazine"),
    ("c1ccnnc1",             "pyridazine"),
])
def test_phase747_azine_2ol_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
