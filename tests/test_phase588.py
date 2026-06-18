"""Phase 588: Substituted pyrido[3,2-g]quinazoline, pyrido[2,3-g]quinazoline,
pyrimido[4,5-b]quinoline, and pyrimido[5,4-b]quinoline naming.
Each has 3 N atoms; substitutable C positions: 7 per compound.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrido[3,2-g]quinazoline (N at 1,3,9; sub C: 2,4-8,10)
    ("c1cnc2cc3ncncc3cc2c1",         "pyrido[3,2-g]quinazoline"),
    ("Cc1ncc2cc3cccnc3cc2n1",        "2-methylpyrido[3,2-g]quinazoline"),
    ("Cc1ncnc2cc3ncccc3cc12",        "4-methylpyrido[3,2-g]quinazoline"),
    ("Cc1c2cccnc2cc2ncncc12",        "5-methylpyrido[3,2-g]quinazoline"),
    ("Cc1ccnc2cc3ncncc3cc12",        "6-methylpyrido[3,2-g]quinazoline"),
    ("Cc1cnc2cc3ncncc3cc2c1",        "7-methylpyrido[3,2-g]quinazoline"),
    ("Cc1ccc2cc3cncnc3cc2n1",        "8-methylpyrido[3,2-g]quinazoline"),
    ("Cc1c2ncccc2cc2cncnc12",        "10-methylpyrido[3,2-g]quinazoline"),
    # pyrido[2,3-g]quinazoline (N at 1,3,6; sub C: 2,4,5,7-10)
    ("c1cnc2cc3cncnc3cc2c1",         "pyrido[2,3-g]quinazoline"),
    ("Cc1ncc2cc3ncccc3cc2n1",        "2-methylpyrido[2,3-g]quinazoline"),
    ("Cc1ncnc2cc3cccnc3cc12",        "4-methylpyrido[2,3-g]quinazoline"),
    ("Cc1c2cncnc2cc2cccnc12",        "5-methylpyrido[2,3-g]quinazoline"),
    ("Cc1ccc2cc3ncncc3cc2n1",        "7-methylpyrido[2,3-g]quinazoline"),
    ("Cc1cnc2cc3cncnc3cc2c1",        "8-methylpyrido[2,3-g]quinazoline"),
    ("Cc1ccnc2cc3cncnc3cc12",        "9-methylpyrido[2,3-g]quinazoline"),
    ("Cc1c2cccnc2cc2cncnc12",        "10-methylpyrido[2,3-g]quinazoline"),
    # pyrimido[4,5-b]quinoline (N at 1,3,10; sub C: 2,4-9)
    ("c1ccc2nc3ncncc3cc2c1",         "pyrimido[4,5-b]quinoline"),
    ("Cc1ncc2cc3ccccc3nc2n1",        "2-methylpyrimido[4,5-b]quinoline"),
    ("Cc1ncnc2nc3ccccc3cc12",        "4-methylpyrimido[4,5-b]quinoline"),
    ("Cc1c2ccccc2nc2ncncc12",        "5-methylpyrimido[4,5-b]quinoline"),
    ("Cc1cccc2nc3ncncc3cc12",        "6-methylpyrimido[4,5-b]quinoline"),
    ("Cc1ccc2nc3ncncc3cc2c1",        "7-methylpyrimido[4,5-b]quinoline"),
    ("Cc1ccc2cc3cncnc3nc2c1",        "8-methylpyrimido[4,5-b]quinoline"),
    ("Cc1cccc2cc3cncnc3nc12",        "9-methylpyrimido[4,5-b]quinoline"),
    # pyrimido[5,4-b]quinoline (N at 1,3,5; sub C: 2,4,6-10)
    ("c1ccc2nc3cncnc3cc2c1",         "pyrimido[5,4-b]quinoline"),
    ("Cc1ncc2nc3ccccc3cc2n1",        "2-methylpyrimido[5,4-b]quinoline"),
    ("Cc1ncnc2cc3ccccc3nc12",        "4-methylpyrimido[5,4-b]quinoline"),
    ("Cc1cccc2cc3ncncc3nc12",        "6-methylpyrimido[5,4-b]quinoline"),
    ("Cc1ccc2cc3ncncc3nc2c1",        "7-methylpyrimido[5,4-b]quinoline"),
    ("Cc1ccc2nc3cncnc3cc2c1",        "8-methylpyrimido[5,4-b]quinoline"),
    ("Cc1cccc2nc3cncnc3cc12",        "9-methylpyrimido[5,4-b]quinoline"),
    ("Cc1c2ccccc2nc2cncnc12",        "10-methylpyrimido[5,4-b]quinoline"),
])
def test_phase588_pyrido_pyrimido_quinazoline(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
