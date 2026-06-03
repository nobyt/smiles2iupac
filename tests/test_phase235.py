"""Phase 235: dicarboxylic acid ester retained names (IUPAC 2013 P-65.1.2.4).

Saturated diacid esters use retained names: malonate, succinate, glutarate,
adipate, etc.  Unsaturated esters still use systematic names.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("COC(=O)C(=O)OC",      "dimethyl oxalate"),
    ("COC(=O)CC(=O)OC",     "dimethyl malonate"),
    ("COC(=O)CCC(=O)OC",    "dimethyl succinate"),
    ("COC(=O)CCCC(=O)OC",   "dimethyl glutarate"),
    ("COC(=O)CCCCC(=O)OC",  "dimethyl adipate"),
    ("CCOC(=O)CC(=O)OCC",   "diethyl malonate"),
    # unsaturated → systematic
    ("COC(=O)C=CC(=O)OC",   "dimethyl but-2-enedioate"),
    # regression: simple monoester unchanged
    ("CCOC(=O)C",            "ethyl acetate"),
])
def test_phase235_diester_retained(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
