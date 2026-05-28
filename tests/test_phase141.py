"""Phase 141: 単環式ヘテロ芳香族 (IUPAC 2013 P-31.1.3)

isoxazole, oxazole, thiazole, isothiazole,
oxadiazoles (1,2,3 / 1,2,4 / 1,2,5 / 1,3,4),
thiadiazoles (1,2,3 / 1,2,4 / 1,2,5 / 1,3,4),
1H-pyrazole, 1H-1,2,3-triazole, 1H-1,2,4-triazole, 1H-tetrazole,
1,2,3-triazine, 1,2,4-triazine, 1,3,5-triazine, 1,2,4,5-tetrazine
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 5-membered with two different heteroatoms
    ("c1cnoc1",   "isoxazole"),
    ("c1cocn1",   "oxazole"),
    ("c1cscn1",   "thiazole"),
    ("c1cnsc1",   "isothiazole"),
    # oxadiazoles
    ("c1conn1",   "1,2,3-oxadiazole"),
    ("c1ncon1",   "1,2,4-oxadiazole"),
    ("c1cnon1",   "1,2,5-oxadiazole"),
    ("c1nnco1",   "1,3,4-oxadiazole"),
    # thiadiazoles
    ("c1csnn1",   "1,2,3-thiadiazole"),
    ("c1ncsn1",   "1,2,4-thiadiazole"),
    ("c1cnsn1",   "1,2,5-thiadiazole"),
    ("c1nncs1",   "1,3,4-thiadiazole"),
    # pyrazole and triazoles
    ("c1cn[nH]c1",  "1H-pyrazole"),
    ("c1cn[nH]n1",  "1H-1,2,3-triazole"),
    ("c1nc[nH]n1",  "1H-1,2,4-triazole"),
    ("c1nn[nH]n1",  "1H-tetrazole"),
    # triazines and tetrazine
    ("c1cnnnc1",   "1,2,3-triazine"),
    ("c1cnncn1",   "1,2,4-triazine"),
    ("c1ncncn1",   "1,3,5-triazine"),
    ("c1nncnn1",   "1,2,4,5-tetrazine"),
    # 回帰: single-N/O/S rings unchanged
    ("c1ccncc1",  "pyridine"),
    ("c1ccoc1",   "furan"),
    ("c1ccsc1",   "thiophene"),
    ("c1cc[nH]c1", "1H-pyrrole"),
    ("c1cnc[nH]1", "1H-imidazole"),
    # 回帰: fused heteroaromatics unchanged
    ("c1ccc2ncccc2c1",  "quinoline"),
    ("c1ccc2[nH]ccc2c1", "1H-indole"),
])
def test_phase141_heteroaromatics(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
