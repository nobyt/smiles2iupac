"""Phase 225: sulfite ester naming (IUPAC 2013 P-67.2).

COS(=O)OC → dimethyl sulfite; COS(=O)O → methyl hydrogen sulfite.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # dialkyl sulfite
    ("COS(=O)OC",   "dimethyl sulfite"),
    ("CCOS(=O)OCC", "diethyl sulfite"),
    ("COS(=O)OCC",  "ethyl methyl sulfite"),
    # monoalkyl (methyl hydrogen sulfite)
    ("COS(=O)O",    "methyl hydrogen sulfite"),
    # regression: dimethyl sulfate still works
    ("COS(=O)(=O)OC", "dimethyl sulfate"),
])
def test_phase225_sulfite_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
