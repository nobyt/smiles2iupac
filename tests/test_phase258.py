"""Phase 258: hydrazinecarboxylate and dithiocarbamic acid (IUPAC 2013).

  NNC(=O)OC  → methyl hydrazinecarboxylate  (N-N detected in carbamate)
  NNC(=O)OCC → ethyl hydrazinecarboxylate
  NC(=S)S    → carbamodithioic acid          (retained name, IUPAC 2013)
  OC(=S)S    → carbonodithioic O-acid
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # hydrazinecarboxylates (carbazates)
    ("NNC(=O)OC",   "methyl hydrazinecarboxylate"),
    ("NNC(=O)OCC",  "ethyl hydrazinecarboxylate"),
    ("NNC(=O)OCCC", "propyl hydrazinecarboxylate"),
    # dithiocarbamic acid
    ("NC(=S)S",     "carbamodithioic acid"),
    ("OC(=S)S",     "carbonodithioic O-acid"),
    # regression: carbamates unchanged
    ("NC(=O)OC",    "methyl carbamate"),
    ("NC(=O)OCC",   "ethyl carbamate"),
    ("CNC(=O)OCC",  "ethyl N-methylcarbamate"),
    # regression: thiourea/urea unchanged
    ("NC(=S)N",     "thiourea"),
    ("NC(=O)N",     "urea"),
])
def test_phase258_hydrazinecarboxylate_dithiocarbamic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
