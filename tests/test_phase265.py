"""Phase 265: benzoate ion naming; semicarbazone false-positive fix (IUPAC 2013).

  [Na+].[O-]C(=O)c1ccccc1   → sodium benzoate   (was 'sodium formate')
  CC=NNC(=O)c1ccccc1         → no longer 'ethanal semicarbazone'
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzoate ion in metal salts (was 'formate' due to aryl-chain bug)
    ("[Na+].[O-]C(=O)c1ccccc1",                 "sodium benzoate"),
    ("[Mg+2].[O-]C(=O)c1ccccc1.[O-]C(=O)c1ccccc1", "magnesium dibenzoate"),
    ("[K+].[O-]C(=O)c1ccccc1",                  "potassium benzoate"),
    # regression: alkyl carboxylates unchanged
    ("[Na+].CC(=O)[O-]",   "sodium acetate"),
    ("[Ca+2].CC(=O)[O-].CC(=O)[O-]", "calcium diacetate"),
    # regression: true semicarbazone still detected correctly
    ("CC=NNC(=O)N",    "ethanal semicarbazone"),
    ("CC(=NNC(=O)N)C", "propan-2-one semicarbazone"),
])
def test_phase265_benzoate_semicarbazone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
