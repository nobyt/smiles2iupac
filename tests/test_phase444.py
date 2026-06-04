"""Phase 444: 1H-pyrrolo[2,3-b]quinoline, thieno/furo-quinoline isomers,
and benzo[g]quinoxaline retained names (IUPAC 2013 P-31.1.3).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-pyrrolo[2,3-b]quinoline
    ("c1ccc2nc3[nH]ccc3cc2c1",  "1H-pyrrolo[2,3-b]quinoline"),
    # thieno[2,3-b]quinoline
    ("c1ccc2nc3sccc3cc2c1",     "thieno[2,3-b]quinoline"),
    # thieno[3,4-b]quinoline
    ("c1ccc2nc3cscc3cc2c1",     "thieno[3,4-b]quinoline"),
    # thieno[3,2-b]quinoline
    ("c1ccc2nc3ccsc3cc2c1",     "thieno[3,2-b]quinoline"),
    # furo[2,3-b]quinoline
    ("c1ccc2nc3occc3cc2c1",     "furo[2,3-b]quinoline"),
    # benzo[g]quinoxaline
    ("c1ccc2cc3nccnc3cc2c1",    "benzo[g]quinoxaline"),
    # regressions
    ("c1ccc2nc3[nH]ccc3cc2c1",  "1H-pyrrolo[2,3-b]quinoline"),
    ("c1cnc2[nH]ccc2c1",        "1H-pyrrolo[2,3-b]pyridine"),
    ("c1cnc2cc[nH]c2c1",        "1H-pyrrolo[3,2-b]pyridine"),
    ("c1ccc2nc3cc[nH]c3cc2c1",  "1H-pyrrolo[3,2-b]quinoline"),
])
def test_phase444_quinoline_hetero(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
