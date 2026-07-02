"""Phase 688: 11H-benzo[a/b]fluorene and 7H-benzo[c]fluorene methyl derivatives (IUPAC 2013 PINs)."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 11H-benzo[a]fluorene: positions 1–10 aromatic, 11 = CH2
    ("c1ccc2c(c1)Cc1c-2ccc2ccccc12",   "11H-benzo[a]fluorene"),
    ("Cc1cccc2c1Cc1c-2ccc2ccccc12",    "10-methyl-11H-benzo[a]fluorene"),
    ("Cc1ccc2c(c1)Cc1c-2ccc2ccccc12",  "9-methyl-11H-benzo[a]fluorene"),
    ("Cc1ccc2c(c1)-c1ccc3ccccc3c1C2",  "8-methyl-11H-benzo[a]fluorene"),
    ("Cc1cccc2c1-c1ccc3ccccc3c1C2",    "7-methyl-11H-benzo[a]fluorene"),
    ("Cc1cc2ccccc2c2c1-c1ccccc1C2",    "6-methyl-11H-benzo[a]fluorene"),
    ("Cc1cc2c(c3ccccc13)Cc1ccccc1-2",  "5-methyl-11H-benzo[a]fluorene"),
    ("Cc1cccc2c3c(ccc12)-c1ccccc1C3",  "4-methyl-11H-benzo[a]fluorene"),
    ("Cc1ccc2c3c(ccc2c1)-c1ccccc1C3",  "3-methyl-11H-benzo[a]fluorene"),
    ("Cc1ccc2ccc3c(c2c1)Cc1ccccc1-3",  "2-methyl-11H-benzo[a]fluorene"),
    ("Cc1cccc2ccc3c(c12)Cc1ccccc1-3",  "1-methyl-11H-benzo[a]fluorene"),
    ("CC1c2ccccc2-c2ccc3ccccc3c21",    "11-methyl-11H-benzo[a]fluorene"),
    # 11H-benzo[b]fluorene: positions 1–10 aromatic, 11 = CH2
    ("c1ccc2c(c1)Cc1cc3ccccc3cc1-2",       "11H-benzo[b]fluorene"),
    ("Cc1cccc2c1Cc1cc3ccccc3cc1-2",        "1-methyl-11H-benzo[b]fluorene"),
    ("Cc1ccc2c(c1)Cc1cc3ccccc3cc1-2",      "2-methyl-11H-benzo[b]fluorene"),
    ("Cc1ccc2c(c1)-c1cc3ccccc3cc1C2",      "3-methyl-11H-benzo[b]fluorene"),
    ("Cc1cccc2c1-c1cc3ccccc3cc1C2",        "4-methyl-11H-benzo[b]fluorene"),
    ("Cc1c2c(cc3ccccc13)Cc1ccccc1-2",      "5-methyl-11H-benzo[b]fluorene"),
    ("Cc1cccc2cc3c(cc12)-c1ccccc1C3",      "6-methyl-11H-benzo[b]fluorene"),
    ("Cc1ccc2cc3c(cc2c1)-c1ccccc1C3",      "7-methyl-11H-benzo[b]fluorene"),
    ("Cc1ccc2cc3c(cc2c1)Cc1ccccc1-3",      "8-methyl-11H-benzo[b]fluorene"),
    ("Cc1cccc2cc3c(cc12)Cc1ccccc1-3",      "9-methyl-11H-benzo[b]fluorene"),
    ("Cc1c2c(cc3ccccc13)-c1ccccc1C2",      "10-methyl-11H-benzo[b]fluorene"),
    ("CC1c2ccccc2-c2cc3ccccc3cc21",        "11-methyl-11H-benzo[b]fluorene"),
    # 7H-benzo[c]fluorene: positions 1–6,8–11 aromatic, 7 = CH2
    ("c1ccc2c(c1)Cc1ccc3ccccc3c1-2",       "7H-benzo[c]fluorene"),
    ("Cc1ccc2ccc3c(c2c1)-c1ccccc1C3",      "2-methyl-7H-benzo[c]fluorene"),
    ("Cc1cccc2ccc3c(c12)-c1ccccc1C3",      "1-methyl-7H-benzo[c]fluorene"),
    ("Cc1cccc2c1-c1c(ccc3ccccc13)C2",      "11-methyl-7H-benzo[c]fluorene"),
    ("Cc1ccc2c(c1)-c1c(ccc3ccccc13)C2",    "10-methyl-7H-benzo[c]fluorene"),
    ("Cc1ccc2c(c1)Cc1ccc3ccccc3c1-2",      "9-methyl-7H-benzo[c]fluorene"),
    ("Cc1cccc2c1Cc1ccc3ccccc3c1-2",        "8-methyl-7H-benzo[c]fluorene"),
    ("CC1c2ccccc2-c2c1ccc1ccccc21",        "7-methyl-7H-benzo[c]fluorene"),
    ("Cc1cc2ccccc2c2c1Cc1ccccc1-2",        "6-methyl-7H-benzo[c]fluorene"),
    ("Cc1cc2c(c3ccccc13)-c1ccccc1C2",      "5-methyl-7H-benzo[c]fluorene"),
    ("Cc1cccc2c3c(ccc12)Cc1ccccc1-3",      "4-methyl-7H-benzo[c]fluorene"),
    ("Cc1ccc2c3c(ccc2c1)Cc1ccccc1-3",      "3-methyl-7H-benzo[c]fluorene"),
])
def test_phase688(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
