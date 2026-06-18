"""Phase 587: Substituted pyrido[3,4-g]quinoxaline, pyrido[2,3-g]quinoxaline,
and pyrazino[2,3-b]quinoline naming.
Each has 3 N atoms; substitutable C positions: 7 per compound.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrido[3,4-g]quinoxaline (N at 1,4,7; sub C: 2,3,5,6,8,9,10)
    ("c1cc2cc3nccnc3cc2cn1",         "pyrido[3,4-g]quinoxaline"),
    ("Cc1cnc2cc3cnccc3cc2n1",        "2-methylpyrido[3,4-g]quinoxaline"),
    ("Cc1cnc2cc3ccncc3cc2n1",        "3-methylpyrido[3,4-g]quinoxaline"),
    ("Cc1c2cnccc2cc2nccnc12",        "5-methylpyrido[3,4-g]quinoxaline"),
    ("Cc1nccc2cc3nccnc3cc12",        "6-methylpyrido[3,4-g]quinoxaline"),
    ("Cc1cc2cc3nccnc3cc2cn1",        "8-methylpyrido[3,4-g]quinoxaline"),
    ("Cc1cncc2cc3nccnc3cc12",        "9-methylpyrido[3,4-g]quinoxaline"),
    ("Cc1c2ccncc2cc2nccnc12",        "10-methylpyrido[3,4-g]quinoxaline"),
    # pyrido[2,3-g]quinoxaline (N at 1,4,6; sub C: 2,3,5,7,8,9,10)
    ("c1cnc2cc3nccnc3cc2c1",         "pyrido[2,3-g]quinoxaline"),
    ("Cc1cnc2cc3ncccc3cc2n1",        "2-methylpyrido[2,3-g]quinoxaline"),
    ("Cc1cnc2cc3cccnc3cc2n1",        "3-methylpyrido[2,3-g]quinoxaline"),
    ("Cc1c2ncccc2cc2nccnc12",        "5-methylpyrido[2,3-g]quinoxaline"),
    ("Cc1ccc2cc3nccnc3cc2n1",        "7-methylpyrido[2,3-g]quinoxaline"),
    ("Cc1cnc2cc3nccnc3cc2c1",        "8-methylpyrido[2,3-g]quinoxaline"),
    ("Cc1ccnc2cc3nccnc3cc12",        "9-methylpyrido[2,3-g]quinoxaline"),
    ("Cc1c2cccnc2cc2nccnc12",        "10-methylpyrido[2,3-g]quinoxaline"),
    # pyrazino[2,3-b]quinoline (N at 1,4,5; sub C: 2,3,6,7,8,9,10)
    ("c1ccc2nc3nccnc3cc2c1",         "pyrazino[2,3-b]quinoline"),
    ("Cc1cnc2nc3ccccc3cc2n1",        "2-methylpyrazino[2,3-b]quinoline"),
    ("Cc1cnc2cc3ccccc3nc2n1",        "3-methylpyrazino[2,3-b]quinoline"),
    ("Cc1cccc2cc3nccnc3nc12",        "6-methylpyrazino[2,3-b]quinoline"),
    ("Cc1ccc2cc3nccnc3nc2c1",        "7-methylpyrazino[2,3-b]quinoline"),
    ("Cc1ccc2nc3nccnc3cc2c1",        "8-methylpyrazino[2,3-b]quinoline"),
    ("Cc1cccc2nc3nccnc3cc12",        "9-methylpyrazino[2,3-b]quinoline"),
    ("Cc1c2ccccc2nc2nccnc12",        "10-methylpyrazino[2,3-b]quinoline"),
])
def test_phase587_pyrido_pyrazino_quinoxaline(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
