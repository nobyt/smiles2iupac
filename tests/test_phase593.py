"""Phase 593: Substituted pteridine, pyrido[2,3-d]pyrimidine,
pyrido[3,4-d]pyrimidine, pyrido[2,3-e]pyrimidine, and
pyrido[3,4-e]pyrimidine naming.
10-atom 6+6 bicyclics with 3-4 N atoms.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pteridine (N at 1,3,5,8; sub C: 2,4,6,7)
    ("c1cnc2ncncc2n1",         "pteridine"),
    ("Cc1ncc2nccnc2n1",        "2-methylpteridine"),
    ("Cc1ncnc2nccnc12",        "4-methylpteridine"),
    ("Cc1cnc2ncncc2n1",        "6-methylpteridine"),
    ("Cc1cnc2cncnc2n1",        "7-methylpteridine"),
    # pyrido[2,3-d]pyrimidine (N at 1,3,8; sub C: 2,4,5,6,7)
    ("c1cnc2ncncc2c1",         "pyrido[2,3-d]pyrimidine"),
    ("Cc1ncc2cccnc2n1",        "2-methylpyrido[2,3-d]pyrimidine"),
    ("Cc1ncnc2ncccc12",        "4-methylpyrido[2,3-d]pyrimidine"),
    ("Cc1ccnc2ncncc12",        "5-methylpyrido[2,3-d]pyrimidine"),
    ("Cc1cnc2ncncc2c1",        "6-methylpyrido[2,3-d]pyrimidine"),
    ("Cc1ccc2cncnc2n1",        "7-methylpyrido[2,3-d]pyrimidine"),
    # pyrido[3,4-d]pyrimidine (N at 1,3,7; sub C: 2,4,5,6,8)
    ("c1cc2cncnc2cn1",         "pyrido[3,4-d]pyrimidine"),
    ("Cc1ncc2ccncc2n1",        "2-methylpyrido[3,4-d]pyrimidine"),
    ("Cc1ncnc2cnccc12",        "4-methylpyrido[3,4-d]pyrimidine"),
    ("Cc1cncc2ncncc12",        "5-methylpyrido[3,4-d]pyrimidine"),
    ("Cc1cc2cncnc2cn1",        "6-methylpyrido[3,4-d]pyrimidine"),
    ("Cc1nccc2cncnc12",        "8-methylpyrido[3,4-d]pyrimidine"),
    # pyrido[2,3-e]pyrimidine (N at 1,3,6; sub C: 2,4,6,7,8)
    ("c1cnc2cncnc2c1",         "pyrido[2,3-e]pyrimidine"),
    ("Cc1ncc2ncccc2n1",        "2-methylpyrido[2,3-e]pyrimidine"),
    ("Cc1ncnc2cccnc12",        "4-methylpyrido[2,3-e]pyrimidine"),
    ("Cc1ccc2ncncc2n1",        "6-methylpyrido[2,3-e]pyrimidine"),
    ("Cc1cnc2cncnc2c1",        "7-methylpyrido[2,3-e]pyrimidine"),
    ("Cc1ccnc2cncnc12",        "8-methylpyrido[2,3-e]pyrimidine"),
    # pyrido[3,4-e]pyrimidine (N at 1,3,9; sub C: 2,4,5,7,8)
    ("c1cc2ncncc2cn1",         "pyrido[3,4-e]pyrimidine"),
    ("Cc1ncc2cnccc2n1",        "2-methylpyrido[3,4-e]pyrimidine"),
    ("Cc1ncnc2ccncc12",        "4-methylpyrido[3,4-e]pyrimidine"),
    ("Cc1nccc2ncncc12",        "5-methylpyrido[3,4-e]pyrimidine"),
    ("Cc1cc2ncncc2cn1",        "7-methylpyrido[3,4-e]pyrimidine"),
    ("Cc1cncc2cncnc12",        "8-methylpyrido[3,4-e]pyrimidine"),
])
def test_phase593_pteridine_pyrido_pyrimidine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
