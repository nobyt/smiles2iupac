"""Phase 434: naphtho[2,1-b]furan and naphtho[2,1-b]thiophene retained names (IUPAC 2013 P-31.1.3).

Two C12H8X tricyclic linear-fused ring systems (furan or thiophene fused at the
2,1-bond of naphthalene) that currently output 'naphthalene' (wrong).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # naphtho[2,1-b]furan — linear tricyclic, O in terminal 5-membered ring
    ("c1ccc2cc3occc3cc2c1",          "naphtho[2,1-b]furan"),
    # naphtho[2,1-b]thiophene — linear tricyclic, S in terminal 5-membered ring
    ("c1ccc2cc3sccc3cc2c1",          "naphtho[2,1-b]thiophene"),
    # regression: benzofuran unchanged
    ("c1ccc2occc2c1",                "benzofuran"),
    # regression: benzo[b]thiophene unchanged
    ("c1ccc2sccc2c1",                "benzo[b]thiophene"),
    # regression: naphthalene unchanged
    ("c1ccc2ccccc2c1",               "naphthalene"),
])
def test_phase434_naphtho_furan_thiophene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
