"""Phase 217: phenylacetate esters and phenylacetic acid derivatives (IUPAC 2013)

  CCOC(=O)Cc1ccccc1   → ethyl 2-phenylacetate
  OC(=O)Cc1ccccc1      → 2-phenylacetic acid
  COC(=O)Cc1ccccc1     → methyl 2-phenylacetate

  collect_acid_chain must stop at aromatic ring boundary so
  benzene ring carbons are not included in the acid chain.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Phenylacetate esters
    ("CCOC(=O)Cc1ccccc1",  "ethyl 2-phenylacetate"),
    ("COC(=O)Cc1ccccc1",   "methyl 2-phenylacetate"),
    # Phenylacetic acid (IUPAC 2013 retained name)
    ("OC(=O)Cc1ccccc1",    "phenylacetic acid"),
    # Regression: simple aliphatic esters still work
    ("CCOC(=O)C",          "ethyl acetate"),
    ("CCOC(=O)CC",         "ethyl propanoate"),
])
def test_phase217_phenylacetate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
