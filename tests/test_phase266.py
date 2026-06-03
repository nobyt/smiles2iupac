"""Phase 266: benzoyloxy substituent and benzenecarbothioic acid naming (IUPAC 2013).

  OC(=O)c1ccc(OC(=O)c2ccccc2)cc1 → 4-(benzoyloxy)benzoic acid   (was 'formyloxy')
  SC(=O)c1ccccc1                   → benzenecarbothioic S-acid    (was 'methanethioic')
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzoyloxy substituent on benzene ring (was wrongly 'formyloxy')
    ("OC(=O)c1ccc(OC(=O)c2ccccc2)cc1",  "4-(benzoyloxy)benzoic acid"),
    # benzoyloxy on aliphatic chain
    ("OC(=O)CCCOC(=O)c1ccccc1",          "4-(benzoyloxy)butanoic acid"),
    # thioic acid with phenyl directly on carbonyl
    ("SC(=O)c1ccccc1",                    "benzenecarbothioic S-acid"),
    # regression: formyloxy still works (HC(=O)O- substituent)
    ("OC(=O)COC=O",                       "2-(formyloxy)acetic acid"),
    # regression: benzoate salt still works
    ("[Na+].[O-]C(=O)c1ccccc1",           "sodium benzoate"),
    # regression: methyl benzoate and ethyl benzoate
    ("COC(=O)c1ccccc1",                   "methyl benzoate"),
    ("CCOC(=O)c1ccccc1",                  "ethyl benzoate"),
    # regression: acetyl substituent (ketone retained as acetoacetic acid)
    ("OC(=O)CC(=O)C",                     "acetoacetic acid"),
])
def test_phase266_benzoyloxy_thioic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
