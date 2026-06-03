"""Phase 211: sulfate and sulfamate ester naming (IUPAC 2013 P-67.2)

  COS(=O)(=O)O   → methyl hydrogen sulfate
  CCOS(=O)(=O)O  → ethyl hydrogen sulfate
  COS(=O)(=O)OC  → dimethyl sulfate
  NS(=O)(=O)OC   → methyl sulfamate

Sulfate monoesters: alkyl hydrogen sulfate
Sulfate diesters: dialkyl sulfate
Sulfamate esters: alkyl sulfamate
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Sulfate monoesters
    ("COS(=O)(=O)O",    "methyl hydrogen sulfate"),
    ("CCOS(=O)(=O)O",   "ethyl hydrogen sulfate"),
    # Sulfate diesters
    ("COS(=O)(=O)OC",   "dimethyl sulfate"),
    ("CCOS(=O)(=O)OCC", "diethyl sulfate"),
    # Sulfamate ester
    ("NS(=O)(=O)OC",    "methyl sulfamate"),
    # regression: sulfonic acid unaffected
    ("CS(=O)(=O)O",     "methanesulfonic acid"),
    # regression: sulfonate ester unaffected
    ("CS(=O)(=O)OC",    "methyl methanesulfonate"),
])
def test_phase211_sulfate_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
