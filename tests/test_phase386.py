"""Phase 386: Sulfonate anion naming (IUPAC 2013 P-65.3.1.3).

C-S(=O)₂-[O⁻] is a sulfonate anion named as [alkyl/aryl]sulfonate,
the deprotonated form of the corresponding sulfonic acid.

Bug: CS(=O)(=O)[O-] was giving 'sulfonylmethane' instead of 'methanesulfonate'.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # simple alkyl sulfonates
    ("CS(=O)(=O)[O-]",             "methanesulfonate"),
    ("CCS(=O)(=O)[O-]",            "ethanesulfonate"),
    ("CCCS(=O)(=O)[O-]",           "propane-1-sulfonate"),
    ("CCCCS(=O)(=O)[O-]",          "butane-1-sulfonate"),
    # aryl sulfonate
    ("c1ccc(cc1)S(=O)(=O)[O-]",    "benzenesulfonate"),
    # sodium salt with sulfonate anion
    ("[Na+].CS(=O)(=O)[O-]",       "sodium methanesulfonate"),
    ("[Na+].CCS(=O)(=O)[O-]",      "sodium ethanesulfonate"),
    # regression: neutral sulfonic acid still works
    ("CS(=O)(=O)O",                "methanesulfonic acid"),
    ("CCS(=O)(=O)O",               "ethanesulfonic acid"),
    # regression: sulfonate ester still works
    ("CS(=O)(=O)OCC",              "ethyl methanesulfonate"),
])
def test_phase386_sulfonate_anion(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
