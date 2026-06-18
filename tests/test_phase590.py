"""Phase 590: Substituted pyrazino[2,3-g][1,8]naphthyridine,
pyrazino[2,3-b][1,5]naphthyridine, pyrazino[2,3-b]quinoxaline,
and pyrazino[2,3-g]quinoxaline naming.
Highly symmetric compounds: pyrazino[2,3-b]quinoxaline has 3 unique sub C,
pyrazino[2,3-g]quinoxaline has 2 unique sub C.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrazino[2,3-g][1,8]naphthyridine (N at 1,4,5,6; sub C: 2,3,7,8,9,10)
    ("c1cnc2nc3nccnc3cc2c1",                  "pyrazino[2,3-g][1,8]naphthyridine"),
    ("Cc1cnc2nc3ncccc3cc2n1",                 "2-methylpyrazino[2,3-g][1,8]naphthyridine"),
    ("Cc1cnc2cc3cccnc3nc2n1",                 "3-methylpyrazino[2,3-g][1,8]naphthyridine"),
    ("Cc1ccc2cc3nccnc3nc2n1",                 "7-methylpyrazino[2,3-g][1,8]naphthyridine"),
    ("Cc1cnc2nc3nccnc3cc2c1",                 "8-methylpyrazino[2,3-g][1,8]naphthyridine"),
    ("Cc1ccnc2nc3nccnc3cc12",                 "9-methylpyrazino[2,3-g][1,8]naphthyridine"),
    ("Cc1c2cccnc2nc2nccnc12",                 "10-methylpyrazino[2,3-g][1,8]naphthyridine"),
    # pyrazino[2,3-b][1,5]naphthyridine (N at 1,4,5,9; sub C: 2,3,6,7,8,10)
    ("c1cnc2cc3nccnc3nc2c1",                  "pyrazino[2,3-b][1,5]naphthyridine"),
    ("Cc1cnc2nc3cccnc3cc2n1",                 "2-methylpyrazino[2,3-b][1,5]naphthyridine"),
    ("Cc1cnc2cc3ncccc3nc2n1",                 "3-methylpyrazino[2,3-b][1,5]naphthyridine"),
    ("Cc1ccnc2cc3nccnc3nc12",                 "6-methylpyrazino[2,3-b][1,5]naphthyridine"),
    ("Cc1cnc2cc3nccnc3nc2c1",                 "7-methylpyrazino[2,3-b][1,5]naphthyridine"),
    ("Cc1ccc2nc3nccnc3cc2n1",                 "8-methylpyrazino[2,3-b][1,5]naphthyridine"),
    ("Cc1c2ncccc2nc2nccnc12",                 "10-methylpyrazino[2,3-b][1,5]naphthyridine"),
    # pyrazino[2,3-b]quinoxaline (N at 1,4,5,10; 3 unique sub C: 2,6,7)
    ("c1ccc2nc3nccnc3nc2c1",                  "pyrazino[2,3-b]quinoxaline"),
    ("Cc1cnc2nc3ccccc3nc2n1",                 "2-methylpyrazino[2,3-b]quinoxaline"),
    ("Cc1cccc2nc3nccnc3nc12",                 "6-methylpyrazino[2,3-b]quinoxaline"),
    ("Cc1ccc2nc3nccnc3nc2c1",                 "7-methylpyrazino[2,3-b]quinoxaline"),
    # pyrazino[2,3-g]quinoxaline (N at 1,4,6,9; 2 unique sub C: 2,5)
    ("c1cnc2cc3nccnc3cc2n1",                  "pyrazino[2,3-g]quinoxaline"),
    ("Cc1cnc2cc3nccnc3cc2n1",                 "2-methylpyrazino[2,3-g]quinoxaline"),
    ("Cc1c2nccnc2cc2nccnc12",                 "5-methylpyrazino[2,3-g]quinoxaline"),
])
def test_phase590_pyrazino_fused(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
