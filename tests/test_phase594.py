"""Phase 594: Substituted pyrido[2,3-e]pyridazine, pyrido[2,3-d]pyridazine,
pyrido[3,4-c]pyridazine, pyrido[3,4-e]pyridazine, pyrido[2,3-c]pyridazine,
1,2,4-benzotriazine, and pyrido[3,4-e]pyrazine naming.
10-atom 6+6 bicyclics with 2-3 N atoms.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrido[2,3-e]pyridazine (N at 1,2,8; sub C: 3,4,6,7,8)
    ("c1cnc2ccnnc2c1",         "pyrido[2,3-e]pyridazine"),
    ("Cc1cc2ncccc2nn1",        "3-methylpyrido[2,3-e]pyridazine"),
    ("Cc1cnnc2cccnc12",        "4-methylpyrido[2,3-e]pyridazine"),
    ("Cc1ccc2nnccc2n1",        "6-methylpyrido[2,3-e]pyridazine"),
    ("Cc1cnc2ccnnc2c1",        "7-methylpyrido[2,3-e]pyridazine"),
    ("Cc1ccnc2ccnnc12",        "8-methylpyrido[2,3-e]pyridazine"),
    # pyrido[2,3-d]pyridazine (N at 1,5,6; sub C: 2,3,4,5,8)
    ("c1cnc2cnncc2c1",         "pyrido[2,3-d]pyridazine"),
    ("Cc1ccc2cnncc2n1",        "2-methylpyrido[2,3-d]pyridazine"),
    ("Cc1cnc2cnncc2c1",        "3-methylpyrido[2,3-d]pyridazine"),
    ("Cc1ccnc2cnncc12",        "4-methylpyrido[2,3-d]pyridazine"),
    ("Cc1nncc2ncccc12",        "5-methylpyrido[2,3-d]pyridazine"),
    ("Cc1nncc2cccnc12",        "8-methylpyrido[2,3-d]pyridazine"),
    # pyrido[3,4-c]pyridazine (N at 1,2,9; sub C: 3,4,5,6,8)
    ("c1cc2ccnnc2cn1",         "pyrido[3,4-c]pyridazine"),
    ("Cc1cc2ccncc2nn1",        "3-methylpyrido[3,4-c]pyridazine"),
    ("Cc1cnnc2cnccc12",        "4-methylpyrido[3,4-c]pyridazine"),
    ("Cc1cncc2nnccc12",        "5-methylpyrido[3,4-c]pyridazine"),
    ("Cc1cc2ccnnc2cn1",        "6-methylpyrido[3,4-c]pyridazine"),
    ("Cc1nccc2ccnnc12",        "8-methylpyrido[3,4-c]pyridazine"),
    # pyrido[3,4-e]pyridazine (N at 1,2,6; sub C: 3,4,5,7,8)
    ("c1cc2nnccc2cn1",         "pyrido[3,4-e]pyridazine"),
    ("Cc1cc2cnccc2nn1",        "3-methylpyrido[3,4-e]pyridazine"),
    ("Cc1cnnc2ccncc12",        "4-methylpyrido[3,4-e]pyridazine"),
    ("Cc1nccc2nnccc12",        "5-methylpyrido[3,4-e]pyridazine"),
    ("Cc1cc2nnccc2cn1",        "7-methylpyrido[3,4-e]pyridazine"),
    ("Cc1cncc2ccnnc12",        "8-methylpyrido[3,4-e]pyridazine"),
    # pyrido[2,3-c]pyridazine (N at 1,2,9; sub C: 3,4,5,6,7)
    ("c1cnc2nnccc2c1",         "pyrido[2,3-c]pyridazine"),
    ("Cc1cc2cccnc2nn1",        "3-methylpyrido[2,3-c]pyridazine"),
    ("Cc1cnnc2ncccc12",        "4-methylpyrido[2,3-c]pyridazine"),
    ("Cc1ccnc2nnccc12",        "5-methylpyrido[2,3-c]pyridazine"),
    ("Cc1cnc2nnccc2c1",        "6-methylpyrido[2,3-c]pyridazine"),
    ("Cc1ccc2ccnnc2n1",        "7-methylpyrido[2,3-c]pyridazine"),
    # 1,2,4-benzotriazine (N at 1,2,4; sub C: 3,5,6,7,8)
    ("c1ccc2nncnc2c1",         "1,2,4-benzotriazine"),
    ("Cc1nnc2ccccc2n1",        "3-methyl-1,2,4-benzotriazine"),
    ("Cc1cccc2nncnc12",        "5-methyl-1,2,4-benzotriazine"),
    ("Cc1ccc2nncnc2c1",        "6-methyl-1,2,4-benzotriazine"),
    ("Cc1ccc2ncnnc2c1",        "7-methyl-1,2,4-benzotriazine"),
    ("Cc1cccc2ncnnc12",        "8-methyl-1,2,4-benzotriazine"),
    # pyrido[3,4-e]pyrazine (N at 1,4,9; sub C: 2,3,5,7,8)
    ("c1cc2nccnc2cn1",         "pyrido[3,4-e]pyrazine"),
    ("Cc1cnc2cnccc2n1",        "2-methylpyrido[3,4-e]pyrazine"),
    ("Cc1cnc2ccncc2n1",        "3-methylpyrido[3,4-e]pyrazine"),
    ("Cc1nccc2nccnc12",        "5-methylpyrido[3,4-e]pyrazine"),
    ("Cc1cc2nccnc2cn1",        "7-methylpyrido[3,4-e]pyrazine"),
    ("Cc1cncc2nccnc12",        "8-methylpyrido[3,4-e]pyrazine"),
])
def test_phase594_pyrido_pyridazine_pyrazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
