"""Phase 571: Substituted acenaphthylene and acenaphthene naming.
Acenaphthylene: C2v-symmetric (12 atoms, C=C bridge), 4 unique pairs:
{1,2}, {3,8}, {4,7}, {5,6}. IUPAC uses minimum locant.
Acenaphthene: C2v-symmetric (12 atoms, CH2-CH2 bridge), same 4 unique pairs.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # acenaphthylene
    ("C1=Cc2cccc3cccc1c23",              "acenaphthylene"),
    ("CC1=Cc2cccc3cccc1c23",             "1-methylacenaphthylene"),
    ("Cc1ccc2cccc3c2c1C=C3",             "3-methylacenaphthylene"),
    ("Cc1cc2c3c(cccc3c1)C=C2",           "4-methylacenaphthylene"),
    ("Cc1ccc2c3c(cccc13)C=C2",           "5-methylacenaphthylene"),
    # acenaphthene
    ("c1cc2c3c(cccc3c1)CC2",             "acenaphthene"),
    ("CC1Cc2cccc3cccc1c23",              "1-methylacenaphthene"),
    ("Cc1ccc2cccc3c2c1CC3",              "3-methylacenaphthene"),
    ("Cc1cc2c3c(cccc3c1)CC2",            "4-methylacenaphthene"),
    ("Cc1ccc2c3c(cccc13)CC2",            "5-methylacenaphthene"),
])
def test_phase571_acenaphtho(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
