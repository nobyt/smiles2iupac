"""Phase 591: Substituted pyrido[2,3-b]quinoxaline, pyrido[2,3-b][1,7]naphthyridine,
pyrido[2,3-b][1,6]naphthyridine, pyrido[2,3-b][1,8]naphthyridine,
and pyrazino[2,3-h][1,6]naphthyridine naming.
pyrido[2,3-b][1,8]naphthyridine is symmetric (4 unique sub C); others have 6-7.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrido[2,3-b]quinoxaline (N at 1,5,10; sub C: 2,3,4,6,7,8,9)
    ("c1ccc2nc3ncccc3nc2c1",                  "pyrido[2,3-b]quinoxaline"),
    ("Cc1ccc2nc3ccccc3nc2n1",                 "2-methylpyrido[2,3-b]quinoxaline"),
    ("Cc1cnc2nc3ccccc3nc2c1",                 "3-methylpyrido[2,3-b]quinoxaline"),
    ("Cc1ccnc2nc3ccccc3nc12",                 "4-methylpyrido[2,3-b]quinoxaline"),
    ("Cc1cccc2nc3ncccc3nc12",                 "6-methylpyrido[2,3-b]quinoxaline"),
    ("Cc1ccc2nc3ncccc3nc2c1",                 "7-methylpyrido[2,3-b]quinoxaline"),
    ("Cc1ccc2nc3cccnc3nc2c1",                 "8-methylpyrido[2,3-b]quinoxaline"),
    ("Cc1cccc2nc3cccnc3nc12",                 "9-methylpyrido[2,3-b]quinoxaline"),
    # pyrido[2,3-b][1,7]naphthyridine (N at 1,8,10; sub C: 2,3,4,5,6,7,9)
    ("c1cnc2nc3cnccc3cc2c1",                  "pyrido[2,3-b][1,7]naphthyridine"),
    ("Cc1ccc2cc3ccncc3nc2n1",                 "2-methylpyrido[2,3-b][1,7]naphthyridine"),
    ("Cc1cnc2nc3cnccc3cc2c1",                 "3-methylpyrido[2,3-b][1,7]naphthyridine"),
    ("Cc1ccnc2nc3cnccc3cc12",                 "4-methylpyrido[2,3-b][1,7]naphthyridine"),
    ("Cc1c2ccncc2nc2ncccc12",                 "5-methylpyrido[2,3-b][1,7]naphthyridine"),
    ("Cc1cncc2nc3ncccc3cc12",                 "6-methylpyrido[2,3-b][1,7]naphthyridine"),
    ("Cc1cc2cc3cccnc3nc2cn1",                 "7-methylpyrido[2,3-b][1,7]naphthyridine"),
    ("Cc1nccc2cc3cccnc3nc12",                 "9-methylpyrido[2,3-b][1,7]naphthyridine"),
    # pyrido[2,3-b][1,6]naphthyridine (N at 1,7,10; sub C: 2,3,4,5,6,8,9)
    ("c1cnc2nc3ccncc3cc2c1",                  "pyrido[2,3-b][1,6]naphthyridine"),
    ("Cc1ccc2cc3cnccc3nc2n1",                 "2-methylpyrido[2,3-b][1,6]naphthyridine"),
    ("Cc1cnc2nc3ccncc3cc2c1",                 "3-methylpyrido[2,3-b][1,6]naphthyridine"),
    ("Cc1ccnc2nc3ccncc3cc12",                 "4-methylpyrido[2,3-b][1,6]naphthyridine"),
    ("Cc1c2cnccc2nc2ncccc12",                 "5-methylpyrido[2,3-b][1,6]naphthyridine"),
    ("Cc1nccc2nc3ncccc3cc12",                 "6-methylpyrido[2,3-b][1,6]naphthyridine"),
    ("Cc1cc2nc3ncccc3cc2cn1",                 "8-methylpyrido[2,3-b][1,6]naphthyridine"),
    ("Cc1cncc2cc3cccnc3nc12",                 "9-methylpyrido[2,3-b][1,6]naphthyridine"),
    # pyrido[2,3-b][1,8]naphthyridine (N at 1,9,10; symm; 4 unique sub C: 2,3,4,5)
    ("c1cnc2nc3ncccc3cc2c1",                  "pyrido[2,3-b][1,8]naphthyridine"),
    ("Cc1ccc2cc3cccnc3nc2n1",                 "2-methylpyrido[2,3-b][1,8]naphthyridine"),
    ("Cc1cnc2nc3ncccc3cc2c1",                 "3-methylpyrido[2,3-b][1,8]naphthyridine"),
    ("Cc1ccnc2nc3ncccc3cc12",                 "4-methylpyrido[2,3-b][1,8]naphthyridine"),
    ("Cc1c2cccnc2nc2ncccc12",                 "5-methylpyrido[2,3-b][1,8]naphthyridine"),
    # pyrazino[2,3-h][1,6]naphthyridine (N at 1,4,5,10; sub C: 2,3,6,7,8,9)
    ("c1cnc2c(c1)cnc1nccnc12",               "pyrazino[2,3-h][1,6]naphthyridine"),
    ("Cc1cnc2ncc3cccnc3c2n1",                "2-methylpyrazino[2,3-h][1,6]naphthyridine"),
    ("Cc1cnc2c(ncc3cccnc32)n1",              "3-methylpyrazino[2,3-h][1,6]naphthyridine"),
    ("Cc1nc2nccnc2c2ncccc12",                "6-methylpyrazino[2,3-h][1,6]naphthyridine"),
    ("Cc1ccnc2c1cnc1nccnc12",               "7-methylpyrazino[2,3-h][1,6]naphthyridine"),
    ("Cc1cnc2c(cnc3nccnc32)c1",              "8-methylpyrazino[2,3-h][1,6]naphthyridine"),
    ("Cc1ccc2cnc3nccnc3c2n1",               "9-methylpyrazino[2,3-h][1,6]naphthyridine"),
])
def test_phase591_pyrido_pyrazino_fused(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
