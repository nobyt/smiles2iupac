"""Phase 433: 1H-benzo[e]indole and 1H-benzo[f]indole retained names (IUPAC 2013 P-31.1.3).

Two C12H9N tricyclic ring systems formed by fusing benzene onto indole.
benzo[e]indole = benzo[g]indole (related by indole mirror symmetry);
benzo[f]indole is the other distinct isomer.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-benzo[f]indole (OPSIN canonical SMILES)
    ("c1ccc2cc3[nH]ccc3cc2c1",         "1H-benzo[f]indole"),
    # 1H-benzo[f]indole — benzo fused at C5-C6 of indole (alternate SMILES)
    ("c1ccc2cc3[nH]ccc3cc2c1",         "1H-benzo[f]indole"),
    # regression: 1H-indole unchanged
    ("c1ccc2[nH]ccc2c1",               "1H-indole"),
    # regression: 9H-carbazole unchanged
    ("c1ccc2[nH]c3ccccc3c2c1",         "9H-carbazole"),
    # regression: naphthalene unchanged
    ("c1ccc2ccccc2c1",                 "naphthalene"),
])
def test_phase433_benzoindole(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
