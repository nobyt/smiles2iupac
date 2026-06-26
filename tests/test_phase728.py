"""Phase 728: methyl derivatives of pyrazolo-fused pyridines/pyridazines and
pyrrolo-fused pyridines/pyridazines:
1H-pyrazolo[3,4-b]pyridine, 1H-pyrazolo[4,3-b]pyridine,
1H-pyrazolo[4,3-c]pyridazine, 1H-pyrazolo[4,5-b]pyridine,
1H-pyrazolo[5,4-c]pyridazine, 1H-pyrazolo[5,4-c]pyridine,
1H-pyrrolo[2,3-b]pyridine, 1H-pyrrolo[3,2-b]pyridine,
1H-pyrrolo[3,4-b]pyridine, 1H-pyrrolo[2,3-c]pyridine,
1H-pyrrolo[3,2-c]pyridine, 1H-pyrrolo[2,3-c]pyridazine,
1H-pyrrolo[3,2-c]pyridazine, 1H-pyrrolo[3,4-c]pyridazine.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-pyrazolo[3,4-b]pyridine (CH at 3,4,5,6; two tautomeric SMILES)
    ("c1cnc2[nH]ncc2c1",        "1H-pyrazolo[3,4-b]pyridine"),
    ("Cc1n[nH]c2ncccc12",       "3-methyl-1H-pyrazolo[3,4-b]pyridine"),
    ("Cc1ccnc2[nH]ncc12",       "4-methyl-1H-pyrazolo[3,4-b]pyridine"),
    ("Cc1cnc2[nH]ncc2c1",       "5-methyl-1H-pyrazolo[3,4-b]pyridine"),
    ("Cc1ccc2cn[nH]c2n1",       "6-methyl-1H-pyrazolo[3,4-b]pyridine"),
    ("c1cnc2n[nH]cc2c1",        "1H-pyrazolo[3,4-b]pyridine"),
    ("Cc1[nH]nc2ncccc12",       "3-methyl-1H-pyrazolo[3,4-b]pyridine"),
    ("Cc1ccnc2n[nH]cc12",       "4-methyl-1H-pyrazolo[3,4-b]pyridine"),
    ("Cc1cnc2n[nH]cc2c1",       "5-methyl-1H-pyrazolo[3,4-b]pyridine"),
    ("Cc1ccc2c[nH]nc2n1",       "6-methyl-1H-pyrazolo[3,4-b]pyridine"),
    # 1H-pyrazolo[4,3-b]pyridine (CH at 3,5,6,7)
    ("c1cnc2c[nH]nc2c1",        "1H-pyrazolo[4,3-b]pyridine"),
    ("Cc1[nH]nc2cccnc12",       "3-methyl-1H-pyrazolo[4,3-b]pyridine"),
    ("Cc1ccc2n[nH]cc2n1",       "5-methyl-1H-pyrazolo[4,3-b]pyridine"),
    ("Cc1cnc2c[nH]nc2c1",       "6-methyl-1H-pyrazolo[4,3-b]pyridine"),
    ("Cc1ccnc2c[nH]nc12",       "7-methyl-1H-pyrazolo[4,3-b]pyridine"),
    # 1H-pyrazolo[4,3-c]pyridazine (CH at 3,6,7)
    ("c1cc2n[nH]cc2nn1",        "1H-pyrazolo[4,3-c]pyridazine"),
    ("Cc1[nH]nc2ccnnc12",       "3-methyl-1H-pyrazolo[4,3-c]pyridazine"),
    ("Cc1cc2n[nH]cc2nn1",       "6-methyl-1H-pyrazolo[4,3-c]pyridazine"),
    ("Cc1cnnc2c[nH]nc12",       "7-methyl-1H-pyrazolo[4,3-c]pyridazine"),
    # 1H-pyrazolo[4,5-b]pyridine (CH at 3,5,6,7)
    ("c1cnc2cn[nH]c2c1",        "1H-pyrazolo[4,5-b]pyridine"),
    ("Cc1n[nH]c2cccnc12",       "3-methyl-1H-pyrazolo[4,5-b]pyridine"),
    ("Cc1ccc2[nH]ncc2n1",       "5-methyl-1H-pyrazolo[4,5-b]pyridine"),
    ("Cc1cnc2cn[nH]c2c1",       "6-methyl-1H-pyrazolo[4,5-b]pyridine"),
    ("Cc1ccnc2cn[nH]c12",       "7-methyl-1H-pyrazolo[4,5-b]pyridine"),
    # 1H-pyrazolo[5,4-c]pyridazine (CH at 3,4,5)
    ("c1cc2cn[nH]c2nn1",        "1H-pyrazolo[5,4-c]pyridazine"),
    ("Cc1n[nH]c2nnccc12",       "3-methyl-1H-pyrazolo[5,4-c]pyridazine"),
    ("Cc1cnnc2[nH]ncc12",       "4-methyl-1H-pyrazolo[5,4-c]pyridazine"),
    ("Cc1cc2cn[nH]c2nn1",       "5-methyl-1H-pyrazolo[5,4-c]pyridazine"),
    # 1H-pyrazolo[5,4-c]pyridine (CH at 3,4,5,7)
    ("c1cc2cn[nH]c2cn1",        "1H-pyrazolo[5,4-c]pyridine"),
    ("Cc1n[nH]c2cnccc12",       "3-methyl-1H-pyrazolo[5,4-c]pyridine"),
    ("Cc1cncc2[nH]ncc12",       "4-methyl-1H-pyrazolo[5,4-c]pyridine"),
    ("Cc1cc2cn[nH]c2cn1",       "5-methyl-1H-pyrazolo[5,4-c]pyridine"),
    ("Cc1nccc2cn[nH]c12",       "7-methyl-1H-pyrazolo[5,4-c]pyridine"),
    # 1H-pyrrolo[2,3-b]pyridine (CH at 2,3,4,5,6)
    ("c1cnc2[nH]ccc2c1",        "1H-pyrrolo[2,3-b]pyridine"),
    ("Cc1cc2cccnc2[nH]1",       "2-methyl-1H-pyrrolo[2,3-b]pyridine"),
    ("Cc1c[nH]c2ncccc12",       "3-methyl-1H-pyrrolo[2,3-b]pyridine"),
    ("Cc1ccnc2[nH]ccc12",       "4-methyl-1H-pyrrolo[2,3-b]pyridine"),
    ("Cc1cnc2[nH]ccc2c1",       "5-methyl-1H-pyrrolo[2,3-b]pyridine"),
    ("Cc1ccc2cc[nH]c2n1",       "6-methyl-1H-pyrrolo[2,3-b]pyridine"),
    # 1H-pyrrolo[3,2-b]pyridine (CH at 2,3,5,6,7)
    ("c1cnc2cc[nH]c2c1",        "1H-pyrrolo[3,2-b]pyridine"),
    ("Cc1cc2ncccc2[nH]1",       "2-methyl-1H-pyrrolo[3,2-b]pyridine"),
    ("Cc1c[nH]c2cccnc12",       "3-methyl-1H-pyrrolo[3,2-b]pyridine"),
    ("Cc1ccc2[nH]ccc2n1",       "5-methyl-1H-pyrrolo[3,2-b]pyridine"),
    ("Cc1cnc2cc[nH]c2c1",       "6-methyl-1H-pyrrolo[3,2-b]pyridine"),
    ("Cc1ccnc2cc[nH]c12",       "7-methyl-1H-pyrrolo[3,2-b]pyridine"),
    # 1H-pyrrolo[3,4-b]pyridine (CH at 2,3,4,5,7)
    ("c1c[nH]c2cncc-2c1",       "1H-pyrrolo[3,4-b]pyridine"),
    ("Cc1ccc2cncc-2[nH]1",      "2-methyl-1H-pyrrolo[3,4-b]pyridine"),
    ("Cc1c[nH]c2cncc-2c1",      "3-methyl-1H-pyrrolo[3,4-b]pyridine"),
    ("Cc1cc[nH]c2cncc1-2",      "4-methyl-1H-pyrrolo[3,4-b]pyridine"),
    ("Cc1ncc2[nH]cccc1-2",      "5-methyl-1H-pyrrolo[3,4-b]pyridine"),
    ("Cc1ncc2ccc[nH]c1-2",      "7-methyl-1H-pyrrolo[3,4-b]pyridine"),
    # 1H-pyrrolo[2,3-c]pyridine (CH at 2,3,4,5,7)
    ("c1cc2cc[nH]c2cn1",        "1H-pyrrolo[2,3-c]pyridine"),
    ("Cc1cc2ccncc2[nH]1",       "2-methyl-1H-pyrrolo[2,3-c]pyridine"),
    ("Cc1c[nH]c2cnccc12",       "3-methyl-1H-pyrrolo[2,3-c]pyridine"),
    ("Cc1cncc2[nH]ccc12",       "4-methyl-1H-pyrrolo[2,3-c]pyridine"),
    ("Cc1cc2cc[nH]c2cn1",       "5-methyl-1H-pyrrolo[2,3-c]pyridine"),
    ("Cc1nccc2cc[nH]c12",       "7-methyl-1H-pyrrolo[2,3-c]pyridine"),
    # 1H-pyrrolo[3,2-c]pyridine (CH at 2,3,4,6,7)
    ("c1cc2[nH]ccc2cn1",        "1H-pyrrolo[3,2-c]pyridine"),
    ("Cc1cc2cnccc2[nH]1",       "2-methyl-1H-pyrrolo[3,2-c]pyridine"),
    ("Cc1c[nH]c2ccncc12",       "3-methyl-1H-pyrrolo[3,2-c]pyridine"),
    ("Cc1nccc2[nH]ccc12",       "4-methyl-1H-pyrrolo[3,2-c]pyridine"),
    ("Cc1cc2[nH]ccc2cn1",       "6-methyl-1H-pyrrolo[3,2-c]pyridine"),
    ("Cc1cncc2cc[nH]c12",       "7-methyl-1H-pyrrolo[3,2-c]pyridine"),
    # 1H-pyrrolo[2,3-c]pyridazine (CH at 3,4,5,6)
    ("c1cc2cc[nH]c2nn1",        "1H-pyrrolo[2,3-c]pyridazine"),
    ("Cc1cc2cc[nH]c2nn1",       "3-methyl-1H-pyrrolo[2,3-c]pyridazine"),
    ("Cc1cnnc2[nH]ccc12",       "4-methyl-1H-pyrrolo[2,3-c]pyridazine"),
    ("Cc1c[nH]c2nnccc12",       "5-methyl-1H-pyrrolo[2,3-c]pyridazine"),
    ("Cc1cc2ccnnc2[nH]1",       "6-methyl-1H-pyrrolo[2,3-c]pyridazine"),
    # 1H-pyrrolo[3,2-c]pyridazine (CH at 3,4,6,7)
    ("c1cc2nccc-2[nH]n1",       "1H-pyrrolo[3,2-c]pyridazine"),
    ("Cc1cc2nccc-2[nH]n1",      "3-methyl-1H-pyrrolo[3,2-c]pyridazine"),
    ("Cc1cn[nH]c2ccnc1-2",      "4-methyl-1H-pyrrolo[3,2-c]pyridazine"),
    ("Cc1cc2[nH]nccc-2n1",      "6-methyl-1H-pyrrolo[3,2-c]pyridazine"),
    ("Cc1cnc2ccn[nH]c1-2",      "7-methyl-1H-pyrrolo[3,2-c]pyridazine"),
    # 1H-pyrrolo[3,4-c]pyridazine (CH at 3,4,5,7)
    ("c1cc2cncc-2[nH]n1",       "1H-pyrrolo[3,4-c]pyridazine"),
    ("Cc1cc2cncc-2[nH]n1",      "3-methyl-1H-pyrrolo[3,4-c]pyridazine"),
    ("Cc1cn[nH]c2cncc1-2",      "4-methyl-1H-pyrrolo[3,4-c]pyridazine"),
    ("Cc1ncc2[nH]nccc1-2",      "5-methyl-1H-pyrrolo[3,4-c]pyridazine"),
    ("Cc1ncc2ccn[nH]c1-2",      "7-methyl-1H-pyrrolo[3,4-c]pyridazine"),
])
def test_phase728(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
