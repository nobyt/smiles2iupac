"""Phase 586: Substituted pyrido[3,4-g]-, pyrido[3,2-g]-, pyrido[2,3-g]-,
and pyrido[4,3-g]quinoline naming.
pyrido[3,2-g]- and pyrido[2,3-g]quinoline are symmetric (equivalent positions
share the same canonical SMILES and map to the lower locant).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # pyrido[3,4-g]quinoline (N at 1,7; sub C: 2-6,8-10)
    ("c1cnc2cc3ccncc3cc2c1",         "pyrido[3,4-g]quinoline"),
    ("Cc1ccc2cc3cnccc3cc2n1",        "2-methylpyrido[3,4-g]quinoline"),
    ("Cc1cnc2cc3ccncc3cc2c1",        "3-methylpyrido[3,4-g]quinoline"),
    ("Cc1ccnc2cc3ccncc3cc12",        "4-methylpyrido[3,4-g]quinoline"),
    ("Cc1c2cnccc2cc2ncccc12",        "5-methylpyrido[3,4-g]quinoline"),
    ("Cc1nccc2cc3ncccc3cc12",        "6-methylpyrido[3,4-g]quinoline"),
    ("Cc1cc2cc3ncccc3cc2cn1",        "8-methylpyrido[3,4-g]quinoline"),
    ("Cc1cncc2cc3cccnc3cc12",        "9-methylpyrido[3,4-g]quinoline"),
    ("Cc1c2ccncc2cc2cccnc12",        "10-methylpyrido[3,4-g]quinoline"),
    # pyrido[3,2-g]quinoline (N at 1,9; symmetric 2=8, 3=7, 4=6; sub C: 2-5,10)
    ("c1cnc2cc3ncccc3cc2c1",         "pyrido[3,2-g]quinoline"),
    ("Cc1ccc2cc3cccnc3cc2n1",        "2-methylpyrido[3,2-g]quinoline"),
    ("Cc1cnc2cc3ncccc3cc2c1",        "3-methylpyrido[3,2-g]quinoline"),
    ("Cc1ccnc2cc3ncccc3cc12",        "4-methylpyrido[3,2-g]quinoline"),
    ("Cc1c2cccnc2cc2ncccc12",        "5-methylpyrido[3,2-g]quinoline"),
    ("Cc1c2ncccc2cc2cccnc12",        "10-methylpyrido[3,2-g]quinoline"),
    # pyrido[2,3-g]quinoline (N at 1,6; symmetric 2=7, 3=8, 4=9, 5=10; sub C: 2-5)
    ("c1cnc2cc3cccnc3cc2c1",         "pyrido[2,3-g]quinoline"),
    ("Cc1ccc2cc3ncccc3cc2n1",        "2-methylpyrido[2,3-g]quinoline"),
    ("Cc1cnc2cc3cccnc3cc2c1",        "3-methylpyrido[2,3-g]quinoline"),
    ("Cc1ccnc2cc3cccnc3cc12",        "4-methylpyrido[2,3-g]quinoline"),
    ("Cc1c2cccnc2cc2cccnc12",        "5-methylpyrido[2,3-g]quinoline"),
    # pyrido[4,3-g]quinoline (N at 1,8; sub C: 2-7,9-10)
    ("c1cnc2cc3cnccc3cc2c1",         "pyrido[4,3-g]quinoline"),
    ("Cc1ccc2cc3ccncc3cc2n1",        "2-methylpyrido[4,3-g]quinoline"),
    ("Cc1cnc2cc3cnccc3cc2c1",        "3-methylpyrido[4,3-g]quinoline"),
    ("Cc1ccnc2cc3cnccc3cc12",        "4-methylpyrido[4,3-g]quinoline"),
    ("Cc1c2ccncc2cc2ncccc12",        "5-methylpyrido[4,3-g]quinoline"),
    ("Cc1cncc2cc3ncccc3cc12",        "6-methylpyrido[4,3-g]quinoline"),
    ("Cc1cc2cc3cccnc3cc2cn1",        "7-methylpyrido[4,3-g]quinoline"),
    ("Cc1nccc2cc3cccnc3cc12",        "9-methylpyrido[4,3-g]quinoline"),
    ("Cc1c2cnccc2cc2cccnc12",        "10-methylpyrido[4,3-g]quinoline"),
])
def test_phase586_pyrido_quinoline(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
