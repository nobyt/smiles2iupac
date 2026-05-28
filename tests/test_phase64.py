"""Phase 64: N-置換イミン命名 (C=N-R)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 一級イミン (C=N-H)
    ("CC=N", "ethan-1-imine"),
    ("CCC=N", "propan-1-imine"),
    ("C=N", "methan-1-imine"),
    # N-単置換イミン (C=N-R)
    ("CC=NC", "N-methylethan-1-imine"),
    ("CCC=NC", "N-methylpropan-1-imine"),
    ("CC=NCC", "N-ethylethan-1-imine"),
    ("C=NC", "N-methylmethan-1-imine"),
    # N-二置換イミン (C=N-R₂)
    ("CC=N(CC)CC", None),  # skip tertiary (quaternary N, invalid)
])
def test_phase64_n_substituted_imine(smiles, expected):
    if expected is None:
        pytest.skip("complex case not yet supported")
    assert smiles_to_iupac(smiles) == expected
