"""Phase 708: methyl derivatives of isoxazolo[x,y-b/c/d/e] series
(28 parent compounds: all b/c/d/e ring-fusion isomers with pyridine,
pyridazine, pyrimidine, pyrazine, and triazine partners).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isoxazolo[3,4-b]pyridine (CH at 3,4,5,6)
    ("c1cnc2nocc2c1",       "isoxazolo[3,4-b]pyridine"),
    ("Cc1onc2ncccc12",      "3-methylisoxazolo[3,4-b]pyridine"),
    ("Cc1ccnc2nocc12",      "4-methylisoxazolo[3,4-b]pyridine"),
    ("Cc1cnc2nocc2c1",      "5-methylisoxazolo[3,4-b]pyridine"),
    ("Cc1ccc2conc2n1",      "6-methylisoxazolo[3,4-b]pyridine"),
    # isoxazolo[3,4-c]pyridazine (CH at 3,4,5)
    ("c1cc2conc2nn1",       "isoxazolo[3,4-c]pyridazine"),
    ("Cc1onc2nnccc12",      "3-methylisoxazolo[3,4-c]pyridazine"),
    ("Cc1cnnc2nocc12",      "4-methylisoxazolo[3,4-c]pyridazine"),
    ("Cc1cc2conc2nn1",      "5-methylisoxazolo[3,4-c]pyridazine"),
    # isoxazolo[3,4-c]pyridine (CH at 3,4,5,7)
    ("c1cc2conc2cn1",       "isoxazolo[3,4-c]pyridine"),
    ("Cc1onc2cnccc12",      "3-methylisoxazolo[3,4-c]pyridine"),
    ("Cc1cncc2nocc12",      "4-methylisoxazolo[3,4-c]pyridine"),
    ("Cc1cc2conc2cn1",      "5-methylisoxazolo[3,4-c]pyridine"),
    ("Cc1nccc2conc12",      "7-methylisoxazolo[3,4-c]pyridine"),
    # isoxazolo[3,4-d][1,2,3]triazine (CH at 4,5)
    ("c1nnnc2nocc12",       "isoxazolo[3,4-d][1,2,3]triazine"),
    ("Cc1nnnc2nocc12",      "4-methylisoxazolo[3,4-d][1,2,3]triazine"),
    ("Cc1onc2nnncc12",      "5-methylisoxazolo[3,4-d][1,2,3]triazine"),
    # isoxazolo[3,4-d]pyridazine (CH at 3,4,7)
    ("c1nncc2nocc12",       "isoxazolo[3,4-d]pyridazine"),
    ("Cc1onc2cnncc12",      "3-methylisoxazolo[3,4-d]pyridazine"),
    ("Cc1nncc2nocc12",      "4-methylisoxazolo[3,4-d]pyridazine"),
    ("Cc1nncc2conc12",      "7-methylisoxazolo[3,4-d]pyridazine"),
    # isoxazolo[3,4-d]pyrimidine (CH at 3,4,6)
    ("c1ncc2conc2n1",       "isoxazolo[3,4-d]pyrimidine"),
    ("Cc1onc2ncncc12",      "3-methylisoxazolo[3,4-d]pyrimidine"),
    ("Cc1ncnc2nocc12",      "4-methylisoxazolo[3,4-d]pyrimidine"),
    ("Cc1ncc2conc2n1",      "6-methylisoxazolo[3,4-d]pyrimidine"),
    # isoxazolo[3,4-e][1,2,4]triazine (CH at 3,7)
    ("c1nnc2conc2n1",       "isoxazolo[3,4-e][1,2,4]triazine"),
    ("Cc1nnc2conc2n1",      "3-methylisoxazolo[3,4-e][1,2,4]triazine"),
    ("Cc1onc2ncnnc12",      "7-methylisoxazolo[3,4-e][1,2,4]triazine"),
    # isoxazolo[3,4-e]pyrazine (CH at 3,5,6)
    ("c1cnc2nocc2n1",       "isoxazolo[3,4-e]pyrazine"),
    ("Cc1onc2nccnc12",      "3-methylisoxazolo[3,4-e]pyrazine"),
    ("Cc1cnc2nocc2n1",      "5-methylisoxazolo[3,4-e]pyrazine"),
    ("Cc1cnc2conc2n1",      "6-methylisoxazolo[3,4-e]pyrazine"),
    # isoxazolo[4,3-b]pyridine (CH at 3,5,6,7)
    ("c1cnc2conc2c1",       "isoxazolo[4,3-b]pyridine"),
    ("Cc1onc2cccnc12",      "3-methylisoxazolo[4,3-b]pyridine"),
    ("Cc1ccc2nocc2n1",      "5-methylisoxazolo[4,3-b]pyridine"),
    ("Cc1cnc2conc2c1",      "6-methylisoxazolo[4,3-b]pyridine"),
    ("Cc1ccnc2conc12",      "7-methylisoxazolo[4,3-b]pyridine"),
    # isoxazolo[4,3-c]pyridazine (CH at 3,6,7)
    ("c1cc2nocc2nn1",       "isoxazolo[4,3-c]pyridazine"),
    ("Cc1onc2ccnnc12",      "3-methylisoxazolo[4,3-c]pyridazine"),
    ("Cc1cc2nocc2nn1",      "6-methylisoxazolo[4,3-c]pyridazine"),
    ("Cc1cnnc2conc12",      "7-methylisoxazolo[4,3-c]pyridazine"),
    # isoxazolo[4,3-c]pyridine (CH at 3,4,6,7)
    ("c1cc2nocc2cn1",       "isoxazolo[4,3-c]pyridine"),
    ("Cc1onc2ccncc12",      "3-methylisoxazolo[4,3-c]pyridine"),
    ("Cc1nccc2nocc12",      "4-methylisoxazolo[4,3-c]pyridine"),
    ("Cc1cc2nocc2cn1",      "6-methylisoxazolo[4,3-c]pyridine"),
    ("Cc1cncc2conc12",      "7-methylisoxazolo[4,3-c]pyridine"),
    # isoxazolo[4,3-d][1,2,3]triazine (CH at 4,7)
    ("c1onc2cnnnc12",       "isoxazolo[4,3-d][1,2,3]triazine"),
    ("Cc1nnnc2conc12",      "4-methylisoxazolo[4,3-d][1,2,3]triazine"),
    ("Cc1onc2cnnnc12",      "7-methylisoxazolo[4,3-d][1,2,3]triazine"),
    # isoxazolo[4,3-d]pyrimidine (CH at 3,5,7)
    ("c1ncc2nocc2n1",       "isoxazolo[4,3-d]pyrimidine"),
    ("Cc1onc2cncnc12",      "3-methylisoxazolo[4,3-d]pyrimidine"),
    ("Cc1ncc2nocc2n1",      "5-methylisoxazolo[4,3-d]pyrimidine"),
    ("Cc1ncnc2conc12",      "7-methylisoxazolo[4,3-d]pyrimidine"),
    # isoxazolo[4,3-e][1,2,4]triazine (CH at 3,5)
    ("c1nnc2nocc2n1",       "isoxazolo[4,3-e][1,2,4]triazine"),
    ("Cc1onc2nncnc12",      "3-methylisoxazolo[4,3-e][1,2,4]triazine"),
    ("Cc1nnc2nocc2n1",      "5-methylisoxazolo[4,3-e][1,2,4]triazine"),
    # isoxazolo[4,5-b]pyridine (CH at 3,5,6,7)
    ("c1cnc2cnoc2c1",       "isoxazolo[4,5-b]pyridine"),
    ("Cc1noc2cccnc12",      "3-methylisoxazolo[4,5-b]pyridine"),
    ("Cc1ccc2oncc2n1",      "5-methylisoxazolo[4,5-b]pyridine"),
    ("Cc1cnc2cnoc2c1",      "6-methylisoxazolo[4,5-b]pyridine"),
    ("Cc1ccnc2cnoc12",      "7-methylisoxazolo[4,5-b]pyridine"),
    # isoxazolo[4,5-c]pyridazine (CH at 3,6,7)
    ("c1cc2oncc2nn1",       "isoxazolo[4,5-c]pyridazine"),
    ("Cc1noc2ccnnc12",      "3-methylisoxazolo[4,5-c]pyridazine"),
    ("Cc1cc2oncc2nn1",      "6-methylisoxazolo[4,5-c]pyridazine"),
    ("Cc1cnnc2cnoc12",      "7-methylisoxazolo[4,5-c]pyridazine"),
    # isoxazolo[4,5-c]pyridine (CH at 3,4,6,7)
    ("c1cc2oncc2cn1",       "isoxazolo[4,5-c]pyridine"),
    ("Cc1noc2ccncc12",      "3-methylisoxazolo[4,5-c]pyridine"),
    ("Cc1nccc2oncc12",      "4-methylisoxazolo[4,5-c]pyridine"),
    ("Cc1cc2oncc2cn1",      "6-methylisoxazolo[4,5-c]pyridine"),
    ("Cc1cncc2cnoc12",      "7-methylisoxazolo[4,5-c]pyridine"),
    # isoxazolo[4,5-d][1,2,3]triazine (CH at 4,7)
    ("c1noc2cnnnc12",       "isoxazolo[4,5-d][1,2,3]triazine"),
    ("Cc1nnnc2cnoc12",      "4-methylisoxazolo[4,5-d][1,2,3]triazine"),
    ("Cc1noc2cnnnc12",      "7-methylisoxazolo[4,5-d][1,2,3]triazine"),
    # isoxazolo[4,5-d]pyrimidine (CH at 3,5,7)
    ("c1ncc2oncc2n1",       "isoxazolo[4,5-d]pyrimidine"),
    ("Cc1noc2cncnc12",      "3-methylisoxazolo[4,5-d]pyrimidine"),
    ("Cc1ncc2oncc2n1",      "5-methylisoxazolo[4,5-d]pyrimidine"),
    ("Cc1ncnc2cnoc12",      "7-methylisoxazolo[4,5-d]pyrimidine"),
    # isoxazolo[4,5-e][1,2,4]triazine (CH at 3,5)
    ("c1nnc2oncc2n1",       "isoxazolo[4,5-e][1,2,4]triazine"),
    ("Cc1noc2nncnc12",      "3-methylisoxazolo[4,5-e][1,2,4]triazine"),
    ("Cc1nnc2oncc2n1",      "5-methylisoxazolo[4,5-e][1,2,4]triazine"),
    # isoxazolo[4,5-e]pyrazine (CH at 3,5,6)
    ("c1cnc2oncc2n1",       "isoxazolo[4,5-e]pyrazine"),
    ("Cc1noc2nccnc12",      "3-methylisoxazolo[4,5-e]pyrazine"),
    ("Cc1cnc2oncc2n1",      "5-methylisoxazolo[4,5-e]pyrazine"),
    ("Cc1cnc2cnoc2n1",      "6-methylisoxazolo[4,5-e]pyrazine"),
    # isoxazolo[5,4-b]pyridine (CH at 3,4,5,6)
    ("c1cnc2oncc2c1",       "isoxazolo[5,4-b]pyridine"),
    ("Cc1noc2ncccc12",      "3-methylisoxazolo[5,4-b]pyridine"),
    ("Cc1ccnc2oncc12",      "4-methylisoxazolo[5,4-b]pyridine"),
    ("Cc1cnc2oncc2c1",      "5-methylisoxazolo[5,4-b]pyridine"),
    ("Cc1ccc2cnoc2n1",      "6-methylisoxazolo[5,4-b]pyridine"),
    # isoxazolo[5,4-c]pyridazine (CH at 3,4,5)
    ("c1cc2cnoc2nn1",       "isoxazolo[5,4-c]pyridazine"),
    ("Cc1noc2nnccc12",      "3-methylisoxazolo[5,4-c]pyridazine"),
    ("Cc1cnnc2oncc12",      "4-methylisoxazolo[5,4-c]pyridazine"),
    ("Cc1cc2cnoc2nn1",      "5-methylisoxazolo[5,4-c]pyridazine"),
    # isoxazolo[5,4-c]pyridine (CH at 3,4,5,7)
    ("c1cc2cnoc2cn1",       "isoxazolo[5,4-c]pyridine"),
    ("Cc1noc2cnccc12",      "3-methylisoxazolo[5,4-c]pyridine"),
    ("Cc1cncc2oncc12",      "4-methylisoxazolo[5,4-c]pyridine"),
    ("Cc1cc2cnoc2cn1",      "5-methylisoxazolo[5,4-c]pyridine"),
    ("Cc1nccc2cnoc12",      "7-methylisoxazolo[5,4-c]pyridine"),
    # isoxazolo[5,4-d][1,2,3]triazine (CH at 4,5)
    ("c1nnnc2oncc12",       "isoxazolo[5,4-d][1,2,3]triazine"),
    ("Cc1nnnc2oncc12",      "4-methylisoxazolo[5,4-d][1,2,3]triazine"),
    ("Cc1noc2nnncc12",      "5-methylisoxazolo[5,4-d][1,2,3]triazine"),
    # isoxazolo[5,4-d]pyridazine (CH at 3,4,7)
    ("c1nncc2oncc12",       "isoxazolo[5,4-d]pyridazine"),
    ("Cc1noc2cnncc12",      "3-methylisoxazolo[5,4-d]pyridazine"),
    ("Cc1nncc2oncc12",      "4-methylisoxazolo[5,4-d]pyridazine"),
    ("Cc1nncc2cnoc12",      "7-methylisoxazolo[5,4-d]pyridazine"),
    # isoxazolo[5,4-d]pyrimidine (CH at 3,4,6)
    ("c1ncc2cnoc2n1",       "isoxazolo[5,4-d]pyrimidine"),
    ("Cc1noc2ncncc12",      "3-methylisoxazolo[5,4-d]pyrimidine"),
    ("Cc1ncnc2oncc12",      "4-methylisoxazolo[5,4-d]pyrimidine"),
    ("Cc1ncc2cnoc2n1",      "6-methylisoxazolo[5,4-d]pyrimidine"),
    # isoxazolo[5,4-e][1,2,4]triazine (CH at 3,7)
    ("c1nnc2cnoc2n1",       "isoxazolo[5,4-e][1,2,4]triazine"),
    ("Cc1nnc2cnoc2n1",      "3-methylisoxazolo[5,4-e][1,2,4]triazine"),
    ("Cc1noc2ncnnc12",      "7-methylisoxazolo[5,4-e][1,2,4]triazine"),
])
def test_phase708(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
