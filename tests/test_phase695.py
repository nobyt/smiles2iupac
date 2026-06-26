"""Phase 695: methyl derivatives of c-fused 9-atom bicyclic heterocycles
(1H-imidazo[4,5-c]pyridine, 1H-pyrazolo[3,4-c]pyridine,
1H-pyrazolo[4,5-c]pyridine, 2H-[1,2,3]triazolo[4,5-c]pyridine,
1H-[1,2,3]triazolo[5,4-c]pyridine).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-imidazo[4,5-c]pyridine (CH at 2,4,6,7)
    ("c1cc2nc[nH]c2cn1",  "1H-imidazo[4,5-c]pyridine"),
    ("Cc1nc2ccncc2[nH]1", "2-methyl-1H-imidazo[4,5-c]pyridine"),
    ("Cc1nccc2nc[nH]c12", "4-methyl-1H-imidazo[4,5-c]pyridine"),
    ("Cc1cc2nc[nH]c2cn1", "6-methyl-1H-imidazo[4,5-c]pyridine"),
    ("Cc1cncc2[nH]cnc12", "7-methyl-1H-imidazo[4,5-c]pyridine"),
    # 1H-pyrazolo[3,4-c]pyridine (CH at 3,4,5,7)
    ("c1cc2c[nH]nc2cn1",  "1H-pyrazolo[3,4-c]pyridine"),
    ("Cc1[nH]nc2cnccc12", "3-methyl-1H-pyrazolo[3,4-c]pyridine"),
    ("Cc1cncc2n[nH]cc12", "4-methyl-1H-pyrazolo[3,4-c]pyridine"),
    ("Cc1cc2c[nH]nc2cn1", "5-methyl-1H-pyrazolo[3,4-c]pyridine"),
    ("Cc1nccc2c[nH]nc12", "7-methyl-1H-pyrazolo[3,4-c]pyridine"),
    # 1H-pyrazolo[4,5-c]pyridine (CH at 3,4,6,7)
    ("c1cc2[nH]ncc2cn1",  "1H-pyrazolo[4,5-c]pyridine"),
    ("Cc1n[nH]c2ccncc12", "3-methyl-1H-pyrazolo[4,5-c]pyridine"),
    ("Cc1nccc2[nH]ncc12", "4-methyl-1H-pyrazolo[4,5-c]pyridine"),
    ("Cc1cc2[nH]ncc2cn1", "6-methyl-1H-pyrazolo[4,5-c]pyridine"),
    ("Cc1cncc2cn[nH]c12", "7-methyl-1H-pyrazolo[4,5-c]pyridine"),
    # 2H-[1,2,3]triazolo[4,5-c]pyridine (CH at 4,6,7)
    ("c1cc2n[nH]nc2cn1",  "2H-[1,2,3]triazolo[4,5-c]pyridine"),
    ("Cc1nccc2n[nH]nc12", "4-methyl-2H-[1,2,3]triazolo[4,5-c]pyridine"),
    ("Cc1cc2n[nH]nc2cn1", "6-methyl-2H-[1,2,3]triazolo[4,5-c]pyridine"),
    ("Cc1cncc2n[nH]nc12", "7-methyl-2H-[1,2,3]triazolo[4,5-c]pyridine"),
    # 1H-[1,2,3]triazolo[5,4-c]pyridine (CH at 4,6,7)
    ("c1cc2nn[nH]c2cn1",  "1H-[1,2,3]triazolo[5,4-c]pyridine"),
    ("Cc1nccc2nn[nH]c12", "4-methyl-1H-[1,2,3]triazolo[5,4-c]pyridine"),
    ("Cc1cc2nn[nH]c2cn1", "6-methyl-1H-[1,2,3]triazolo[5,4-c]pyridine"),
    ("Cc1cncc2[nH]nnc12", "7-methyl-1H-[1,2,3]triazolo[5,4-c]pyridine"),
])
def test_phase695(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
