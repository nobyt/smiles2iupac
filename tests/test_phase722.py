"""Phase 722: methyl derivatives of pyrido[x,y-c/d/e]pyrimidine,
pyrido[x,y-c/d/e]pyridazine, and pyrido[x,y-e]pyrazine series.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrido[2,3-d]pyrimidine (CH at 2,4,5,6,7)
    ("c1cnc2ncncc2c1",             "pyrido[2,3-d]pyrimidine"),
    ("Cc1ncc2cccnc2n1",            "2-methylpyrido[2,3-d]pyrimidine"),
    ("Cc1ncnc2ncccc12",            "4-methylpyrido[2,3-d]pyrimidine"),
    ("Cc1ccnc2ncncc12",            "5-methylpyrido[2,3-d]pyrimidine"),
    ("Cc1cnc2ncncc2c1",            "6-methylpyrido[2,3-d]pyrimidine"),
    ("Cc1ccc2cncnc2n1",            "7-methylpyrido[2,3-d]pyrimidine"),
    # pyrido[3,4-d]pyrimidine (CH at 2,4,5,6,8)
    ("c1cc2cncnc2cn1",             "pyrido[3,4-d]pyrimidine"),
    ("Cc1ncc2ccncc2n1",            "2-methylpyrido[3,4-d]pyrimidine"),
    ("Cc1ncnc2cnccc12",            "4-methylpyrido[3,4-d]pyrimidine"),
    ("Cc1cncc2ncncc12",            "5-methylpyrido[3,4-d]pyrimidine"),
    ("Cc1cc2cncnc2cn1",            "6-methylpyrido[3,4-d]pyrimidine"),
    ("Cc1nccc2cncnc12",            "8-methylpyrido[3,4-d]pyrimidine"),
    # pyrido[2,3-e]pyrimidine (CH at 2,4,6,7,8)
    ("c1cnc2cncnc2c1",             "pyrido[2,3-e]pyrimidine"),
    ("Cc1ncc2ncccc2n1",            "2-methylpyrido[2,3-e]pyrimidine"),
    ("Cc1ncnc2cccnc12",            "4-methylpyrido[2,3-e]pyrimidine"),
    ("Cc1ccc2ncncc2n1",            "6-methylpyrido[2,3-e]pyrimidine"),
    ("Cc1cnc2cncnc2c1",            "7-methylpyrido[2,3-e]pyrimidine"),
    ("Cc1ccnc2cncnc12",            "8-methylpyrido[2,3-e]pyrimidine"),
    # pyrido[3,4-e]pyrimidine (CH at 2,4,5,7,8)
    ("c1cc2ncncc2cn1",             "pyrido[3,4-e]pyrimidine"),
    ("Cc1ncc2cnccc2n1",            "2-methylpyrido[3,4-e]pyrimidine"),
    ("Cc1ncnc2ccncc12",            "4-methylpyrido[3,4-e]pyrimidine"),
    ("Cc1nccc2ncncc12",            "5-methylpyrido[3,4-e]pyrimidine"),
    ("Cc1cc2ncncc2cn1",            "7-methylpyrido[3,4-e]pyrimidine"),
    ("Cc1cncc2cncnc12",            "8-methylpyrido[3,4-e]pyrimidine"),
    # pyrido[2,3-e]pyridazine (CH at 3,4,6,7,8)
    ("c1cnc2ccnnc2c1",             "pyrido[2,3-e]pyridazine"),
    ("Cc1cc2ncccc2nn1",            "3-methylpyrido[2,3-e]pyridazine"),
    ("Cc1cnnc2cccnc12",            "4-methylpyrido[2,3-e]pyridazine"),
    ("Cc1ccc2nnccc2n1",            "6-methylpyrido[2,3-e]pyridazine"),
    ("Cc1cnc2ccnnc2c1",            "7-methylpyrido[2,3-e]pyridazine"),
    ("Cc1ccnc2ccnnc12",            "8-methylpyrido[2,3-e]pyridazine"),
    # pyrido[2,3-d]pyridazine (CH at 2,3,4,5,8)
    ("c1cnc2cnncc2c1",             "pyrido[2,3-d]pyridazine"),
    ("Cc1ccc2cnncc2n1",            "2-methylpyrido[2,3-d]pyridazine"),
    ("Cc1cnc2cnncc2c1",            "3-methylpyrido[2,3-d]pyridazine"),
    ("Cc1ccnc2cnncc12",            "4-methylpyrido[2,3-d]pyridazine"),
    ("Cc1nncc2ncccc12",            "5-methylpyrido[2,3-d]pyridazine"),
    ("Cc1nncc2cccnc12",            "8-methylpyrido[2,3-d]pyridazine"),
    # pyrido[3,4-c]pyridazine (CH at 3,4,5,6,8)
    ("c1cc2ccnnc2cn1",             "pyrido[3,4-c]pyridazine"),
    ("Cc1cc2ccncc2nn1",            "3-methylpyrido[3,4-c]pyridazine"),
    ("Cc1cnnc2cnccc12",            "4-methylpyrido[3,4-c]pyridazine"),
    ("Cc1cncc2nnccc12",            "5-methylpyrido[3,4-c]pyridazine"),
    ("Cc1cc2ccnnc2cn1",            "6-methylpyrido[3,4-c]pyridazine"),
    ("Cc1nccc2ccnnc12",            "8-methylpyrido[3,4-c]pyridazine"),
    # pyrido[3,4-e]pyridazine (CH at 3,4,5,7,8)
    ("c1cc2nnccc2cn1",             "pyrido[3,4-e]pyridazine"),
    ("Cc1cc2cnccc2nn1",            "3-methylpyrido[3,4-e]pyridazine"),
    ("Cc1cnnc2ccncc12",            "4-methylpyrido[3,4-e]pyridazine"),
    ("Cc1nccc2nnccc12",            "5-methylpyrido[3,4-e]pyridazine"),
    ("Cc1cc2nnccc2cn1",            "7-methylpyrido[3,4-e]pyridazine"),
    ("Cc1cncc2ccnnc12",            "8-methylpyrido[3,4-e]pyridazine"),
    # pyrido[2,3-c]pyridazine (CH at 3,4,5,6,7)
    ("c1cnc2nnccc2c1",             "pyrido[2,3-c]pyridazine"),
    ("Cc1cc2cccnc2nn1",            "3-methylpyrido[2,3-c]pyridazine"),
    ("Cc1cnnc2ncccc12",            "4-methylpyrido[2,3-c]pyridazine"),
    ("Cc1ccnc2nnccc12",            "5-methylpyrido[2,3-c]pyridazine"),
    ("Cc1cnc2nnccc2c1",            "6-methylpyrido[2,3-c]pyridazine"),
    ("Cc1ccc2ccnnc2n1",            "7-methylpyrido[2,3-c]pyridazine"),
    # pyrido[3,4-e]pyrazine (CH at 2,3,5,7,8)
    ("c1cc2nccnc2cn1",             "pyrido[3,4-e]pyrazine"),
    ("Cc1cnc2cnccc2n1",            "2-methylpyrido[3,4-e]pyrazine"),
    ("Cc1cnc2ccncc2n1",            "3-methylpyrido[3,4-e]pyrazine"),
    ("Cc1nccc2nccnc12",            "5-methylpyrido[3,4-e]pyrazine"),
    ("Cc1cc2nccnc2cn1",            "7-methylpyrido[3,4-e]pyrazine"),
    ("Cc1cncc2nccnc12",            "8-methylpyrido[3,4-e]pyrazine"),
])
def test_phase722(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
