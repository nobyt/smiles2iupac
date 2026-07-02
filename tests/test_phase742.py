"""Phase 742: benzo[a/b/c]fluorene → 11H/7H-benzo[x]fluorene (IUPAC 2013 PINs).

Retained names benzo[a/b]fluorene → 11H-benzo[a/b]fluorene;
benzo[c]fluorene → 7H-benzo[c]fluorene.
benzo[a] also requires locant renumbering (reversal of positions 1–10).
benzo[c] requires a full locant remap.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Parent names
    ("c1ccc2c(c1)Cc1c-2ccc2ccccc12",   "11H-benzo[a]fluorene"),
    ("c1ccc2c(c1)Cc1cc3ccccc3cc1-2",   "11H-benzo[b]fluorene"),
    ("c1ccc2c(c1)Cc1ccc3ccccc3c1-2",   "7H-benzo[c]fluorene"),
    # benzo[a]: locant reversal (1↔10, 2↔9, …, 5↔6; 11 unchanged)
    ("Cc1cccc2ccc3c(c12)Cc1ccccc1-3",  "1-methyl-11H-benzo[a]fluorene"),
    ("Cc1cccc2c1Cc1c-2ccc2ccccc12",    "10-methyl-11H-benzo[a]fluorene"),
    ("CC1c2ccccc2-c2ccc3ccccc3c21",    "11-methyl-11H-benzo[a]fluorene"),
    # benzo[b]: locants unchanged, add 11H prefix
    ("Cc1cccc2c1Cc1cc3ccccc3cc1-2",    "1-methyl-11H-benzo[b]fluorene"),
    ("Cc1c2c(cc3ccccc13)-c1ccccc1C2",  "10-methyl-11H-benzo[b]fluorene"),
    ("CC1c2ccccc2-c2cc3ccccc3cc21",    "11-methyl-11H-benzo[b]fluorene"),
    # benzo[c]: complex remap; CH2 stays at 7
    ("Cc1cccc2ccc3c(c12)-c1ccccc1C3",  "1-methyl-7H-benzo[c]fluorene"),
    ("Cc1ccc2ccc3c(c2c1)-c1ccccc1C3",  "2-methyl-7H-benzo[c]fluorene"),
    ("Cc1ccc2c3c(ccc2c1)Cc1ccccc1-3",  "3-methyl-7H-benzo[c]fluorene"),
    ("CC1c2ccccc2-c2c1ccc1ccccc21",    "7-methyl-7H-benzo[c]fluorene"),
    # Regression: 9H-fluorene unchanged
    ("c1ccc2c(c1)Cc1ccccc1-2",         "9H-fluorene"),
    ("Cc1cccc2c1Cc1ccccc1-2",          "1-methyl-9H-fluorene"),
])
def test_phase742_benzofluorene_pins(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
