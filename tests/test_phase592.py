"""Phase 592: Substituted pyrimido[4,5-b]quinoxaline naming.
4 N atoms (N at 1,3,5,10); 6 substitutable C: 2,4,6,7,8,9.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrimido[4,5-b]quinoxaline (N at 1,3,5,10; sub C: 2,4,6,7,8,9)
    ("c1ccc2nc3ncncc3nc2c1",                  "pyrimido[4,5-b]quinoxaline"),
    ("Cc1ncc2nc3ccccc3nc2n1",                 "2-methylpyrimido[4,5-b]quinoxaline"),
    ("Cc1ncnc2nc3ccccc3nc12",                 "4-methylpyrimido[4,5-b]quinoxaline"),
    ("Cc1cccc2nc3ncncc3nc12",                 "6-methylpyrimido[4,5-b]quinoxaline"),
    ("Cc1ccc2nc3ncncc3nc2c1",                 "7-methylpyrimido[4,5-b]quinoxaline"),
    ("Cc1ccc2nc3cncnc3nc2c1",                 "8-methylpyrimido[4,5-b]quinoxaline"),
    ("Cc1cccc2nc3cncnc3nc12",                 "9-methylpyrimido[4,5-b]quinoxaline"),
])
def test_phase592_pyrimido_quinoxaline(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
