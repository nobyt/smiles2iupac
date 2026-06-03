"""Phase 239: carbon tetrachloride retained name + trione naming.

IUPAC 2013 P-31.1.3.4: carbon tetrachloride = CCl4 (retained name).
Ketone triones: pentane-2,3,4-trione etc.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # carbon tetrachloride retained name
    ("ClC(Cl)(Cl)Cl",     "carbon tetrachloride"),
    # triones
    ("CC(=O)C(=O)C(=O)C", "pentane-2,3,4-trione"),
    ("CC(=O)C(=O)C(=O)CC","hexane-2,3,4-trione"),
    # regression: dione unchanged
    ("CC(=O)CC(=O)C",     "pentane-2,4-dione"),
    ("CC(=O)C(=O)C",      "butane-2,3-dione"),
    # regression: other halomethanes unchanged
    ("ClCCl",             "dichloromethane"),
    ("ClC(Cl)Cl",         "trichloromethane"),
])
def test_phase239_ccl4_trione(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
