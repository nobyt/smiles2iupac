"""Phase 589: Substituted pyrimido[5,4-g][1,5]naphthyridine,
pyrimido[5,4-g][1,8]naphthyridine, and pyrido[2,3-g][1,5]naphthyridine naming.
4 N atoms in first two compounds (6 sub C each); 3 N atoms in third (7 sub C).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrimido[5,4-g][1,5]naphthyridine (N at 1,3,5,9; sub C: 2,4,6,7,8,10)
    ("c1cnc2cc3ncncc3nc2c1",                  "pyrimido[5,4-g][1,5]naphthyridine"),
    ("Cc1ncc2nc3cccnc3cc2n1",                 "2-methylpyrimido[5,4-g][1,5]naphthyridine"),
    ("Cc1ncnc2cc3ncccc3nc12",                 "4-methylpyrimido[5,4-g][1,5]naphthyridine"),
    ("Cc1ccnc2cc3ncncc3nc12",                 "6-methylpyrimido[5,4-g][1,5]naphthyridine"),
    ("Cc1cnc2cc3ncncc3nc2c1",                 "7-methylpyrimido[5,4-g][1,5]naphthyridine"),
    ("Cc1ccc2nc3cncnc3cc2n1",                 "8-methylpyrimido[5,4-g][1,5]naphthyridine"),
    ("Cc1c2ncccc2nc2cncnc12",                 "10-methylpyrimido[5,4-g][1,5]naphthyridine"),
    # pyrimido[5,4-g][1,8]naphthyridine (N at 1,3,9,10; sub C: 2,4-8)
    ("c1cnc2nc3ncncc3cc2c1",                  "pyrimido[5,4-g][1,8]naphthyridine"),
    ("Cc1ncc2cc3cccnc3nc2n1",                 "2-methylpyrimido[5,4-g][1,8]naphthyridine"),
    ("Cc1ncnc2nc3ncccc3cc12",                 "4-methylpyrimido[5,4-g][1,8]naphthyridine"),
    ("Cc1c2cccnc2nc2ncncc12",                 "5-methylpyrimido[5,4-g][1,8]naphthyridine"),
    ("Cc1ccnc2nc3ncncc3cc12",                 "6-methylpyrimido[5,4-g][1,8]naphthyridine"),
    ("Cc1cnc2nc3ncncc3cc2c1",                 "7-methylpyrimido[5,4-g][1,8]naphthyridine"),
    ("Cc1ccc2cc3cncnc3nc2n1",                 "8-methylpyrimido[5,4-g][1,8]naphthyridine"),
    # pyrido[2,3-g][1,5]naphthyridine (N at 1,5,6; sub C: 2-4,7-10)
    ("c1cnc2nc3cccnc3cc2c1",                  "pyrido[2,3-g][1,5]naphthyridine"),
    ("Cc1ccc2nc3ncccc3cc2n1",                 "2-methylpyrido[2,3-g][1,5]naphthyridine"),
    ("Cc1cnc2cc3cccnc3nc2c1",                 "3-methylpyrido[2,3-g][1,5]naphthyridine"),
    ("Cc1ccnc2cc3cccnc3nc12",                 "4-methylpyrido[2,3-g][1,5]naphthyridine"),
    ("Cc1ccc2cc3ncccc3nc2n1",                 "7-methylpyrido[2,3-g][1,5]naphthyridine"),
    ("Cc1cnc2nc3cccnc3cc2c1",                 "8-methylpyrido[2,3-g][1,5]naphthyridine"),
    ("Cc1ccnc2nc3cccnc3cc12",                 "9-methylpyrido[2,3-g][1,5]naphthyridine"),
    ("Cc1c2cccnc2nc2cccnc12",                 "10-methylpyrido[2,3-g][1,5]naphthyridine"),
])
def test_phase589_pyrimido_pyrido_naphthyridine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
