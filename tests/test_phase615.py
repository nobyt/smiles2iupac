"""Phase 615: azulene (5-7 fused PAH) and chrysene (4-ring PAH) with methyl substituents."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # azulene (C2v; 1=3, 4=8, 5=7 by symmetry → lowest locant wins)
    ('c1ccc2cccc-2cc1',    'azulene'),
    ('Cc1ccc2cccccc1-2',   '1-methylazulene'),
    ('Cc1cc2cccccc-2c1',   '2-methylazulene'),
    ('Cc1ccccc2cccc1-2',   '4-methylazulene'),
    ('Cc1cccc2cccc-2c1',   '5-methylazulene'),
    ('Cc1ccc2cccc-2cc1',   '6-methylazulene'),
    # chrysene (D2h; 1=7, 2=8, 3=9, 4=10, 5=11, 6=12 → lowest locant wins)
    ('c1ccc2c(c1)ccc1c3ccccc3ccc21',   'chrysene'),
    ('Cc1cccc2c1ccc1c3ccccc3ccc21',    '1-methylchrysene'),
    ('Cc1ccc2c(ccc3c4ccccc4ccc23)c1',  '2-methylchrysene'),
    ('Cc1ccc2ccc3c4ccccc4ccc3c2c1',    '3-methylchrysene'),
    ('Cc1cccc2ccc3c4ccccc4ccc3c12',    '4-methylchrysene'),
    ('Cc1cc2ccccc2c2ccc3ccccc3c12',    '5-methylchrysene'),
    ('Cc1cc2c3ccccc3ccc2c2ccccc12',    '6-methylchrysene'),
])
def test_phase615_pah(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
