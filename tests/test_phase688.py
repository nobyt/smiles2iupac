"""Phase 688: benzo[a/b/c]fluorene methyl derivatives (positions 1–11)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[a]fluorene (11H): positions 1–10 aromatic, 11 = CH2
    ("c1ccc2c(c1)Cc1c-2ccc2ccccc12",   "benzo[a]fluorene"),
    ("Cc1cccc2c1Cc1c-2ccc2ccccc12",    "1-methylbenzo[a]fluorene"),
    ("Cc1ccc2c(c1)Cc1c-2ccc2ccccc12",  "2-methylbenzo[a]fluorene"),
    ("Cc1ccc2c(c1)-c1ccc3ccccc3c1C2",  "3-methylbenzo[a]fluorene"),
    ("Cc1cccc2c1-c1ccc3ccccc3c1C2",    "4-methylbenzo[a]fluorene"),
    ("Cc1cc2ccccc2c2c1-c1ccccc1C2",    "5-methylbenzo[a]fluorene"),
    ("Cc1cc2c(c3ccccc13)Cc1ccccc1-2",  "6-methylbenzo[a]fluorene"),
    ("Cc1cccc2c3c(ccc12)-c1ccccc1C3",  "7-methylbenzo[a]fluorene"),
    ("Cc1ccc2c3c(ccc2c1)-c1ccccc1C3",  "8-methylbenzo[a]fluorene"),
    ("Cc1ccc2ccc3c(c2c1)Cc1ccccc1-3",  "9-methylbenzo[a]fluorene"),
    ("Cc1cccc2ccc3c(c12)Cc1ccccc1-3",  "10-methylbenzo[a]fluorene"),
    ("CC1c2ccccc2-c2ccc3ccccc3c21",    "11-methylbenzo[a]fluorene"),
    # benzo[b]fluorene (11H): positions 1–10 aromatic, 11 = CH2
    ("c1ccc2c(c1)Cc1cc3ccccc3cc1-2",       "benzo[b]fluorene"),
    ("Cc1cccc2c1Cc1cc3ccccc3cc1-2",        "1-methylbenzo[b]fluorene"),
    ("Cc1ccc2c(c1)Cc1cc3ccccc3cc1-2",      "2-methylbenzo[b]fluorene"),
    ("Cc1ccc2c(c1)-c1cc3ccccc3cc1C2",      "3-methylbenzo[b]fluorene"),
    ("Cc1cccc2c1-c1cc3ccccc3cc1C2",        "4-methylbenzo[b]fluorene"),
    ("Cc1c2c(cc3ccccc13)Cc1ccccc1-2",      "5-methylbenzo[b]fluorene"),
    ("Cc1cccc2cc3c(cc12)-c1ccccc1C3",      "6-methylbenzo[b]fluorene"),
    ("Cc1ccc2cc3c(cc2c1)-c1ccccc1C3",      "7-methylbenzo[b]fluorene"),
    ("Cc1ccc2cc3c(cc2c1)Cc1ccccc1-3",      "8-methylbenzo[b]fluorene"),
    ("Cc1cccc2cc3c(cc12)Cc1ccccc1-3",      "9-methylbenzo[b]fluorene"),
    ("Cc1c2c(cc3ccccc13)-c1ccccc1C2",      "10-methylbenzo[b]fluorene"),
    ("CC1c2ccccc2-c2cc3ccccc3cc21",        "11-methylbenzo[b]fluorene"),
    # benzo[c]fluorene (7H): positions 1–6,8–11 aromatic, 7 = CH2
    ("c1ccc2c(c1)Cc1ccc3ccccc3c1-2",       "benzo[c]fluorene"),
    ("Cc1ccc2ccc3c(c2c1)-c1ccccc1C3",      "1-methylbenzo[c]fluorene"),
    ("Cc1cccc2ccc3c(c12)-c1ccccc1C3",      "2-methylbenzo[c]fluorene"),
    ("Cc1cccc2c1-c1c(ccc3ccccc13)C2",      "3-methylbenzo[c]fluorene"),
    ("Cc1ccc2c(c1)-c1c(ccc3ccccc13)C2",    "4-methylbenzo[c]fluorene"),
    ("Cc1ccc2c(c1)Cc1ccc3ccccc3c1-2",      "5-methylbenzo[c]fluorene"),
    ("Cc1cccc2c1Cc1ccc3ccccc3c1-2",        "6-methylbenzo[c]fluorene"),
    ("CC1c2ccccc2-c2c1ccc1ccccc21",        "7-methylbenzo[c]fluorene"),
    ("Cc1cc2ccccc2c2c1Cc1ccccc1-2",        "8-methylbenzo[c]fluorene"),
    ("Cc1cc2c(c3ccccc13)-c1ccccc1C2",      "9-methylbenzo[c]fluorene"),
    ("Cc1cccc2c3c(ccc12)Cc1ccccc1-3",      "10-methylbenzo[c]fluorene"),
    ("Cc1ccc2c3c(ccc2c1)Cc1ccccc1-3",      "11-methylbenzo[c]fluorene"),
])
def test_phase688(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
