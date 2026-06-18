"""Phase 626: 2,3-dihydrobenzofuran, 2,3-dihydrobenzothiophene, and their 1,3-dihydro isomers
with substituents at all positions (IUPAC 2013).

Also fixes the previously wrong retained names:
  c1ccc2c(c1)CCO2 was '1,3-dihydro-2-benzofuran' → now correctly '2,3-dihydrobenzofuran'
  c1ccc2c(c1)CCS2 was '1,3-dihydro-2-benzothiophene' → now correctly '2,3-dihydrobenzothiophene'
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("c1ccc2c(c1)CCO2", "2,3-dihydrobenzofuran"),
    ("c1ccc2c(c1)COC2", "1,3-dihydro-2-benzofuran"),
    ("c1ccc2c(c1)CCS2", "2,3-dihydrobenzothiophene"),
    ("c1ccc2c(c1)CSC2", "1,3-dihydro-2-benzothiophene"),
    # 2,3-dihydrobenzofuran substituted (O at 1, CH2 at 2,3, benzene at 4-7)
    ("Cc1cccc2c1CCO2",  "4-methyl-2,3-dihydrobenzofuran"),
    ("Cc1ccc2c(c1)CCO2", "5-methyl-2,3-dihydrobenzofuran"),
    ("Cc1ccc2c(c1)OCC2", "6-methyl-2,3-dihydrobenzofuran"),
    ("Cc1cccc2c1OCC2",  "7-methyl-2,3-dihydrobenzofuran"),
    ("CC1Cc2ccccc2O1",  "2-methyl-2,3-dihydrobenzofuran"),
    ("CC1COc2ccccc21",  "3-methyl-2,3-dihydrobenzofuran"),
    # 2,3-dihydrobenzothiophene substituted (S at 1, CH2 at 2,3, benzene at 4-7)
    ("Cc1cccc2c1CCS2",  "4-methyl-2,3-dihydrobenzothiophene"),
    ("Cc1ccc2c(c1)CCS2", "5-methyl-2,3-dihydrobenzothiophene"),
    ("Cc1ccc2c(c1)SCC2", "6-methyl-2,3-dihydrobenzothiophene"),
    ("Cc1cccc2c1SCC2",  "7-methyl-2,3-dihydrobenzothiophene"),
    # 1,3-dihydro-2-benzofuran substituted (C at 1,3 symmetric; O at 2)
    ("Cc1cccc2c1COC2",  "4-methyl-1,3-dihydro-2-benzofuran"),
    ("Cc1ccc2c(c1)COC2", "5-methyl-1,3-dihydro-2-benzofuran"),
    ("CC1OCc2ccccc21",  "1-methyl-1,3-dihydro-2-benzofuran"),
])
def test_phase626_dihydrobenzo_heterocycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
