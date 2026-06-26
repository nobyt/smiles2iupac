"""Phase 690: methyl derivatives of acridine, phenanthridine, phenazine, 9H-carbazole, phenothiazine, phenoxazine."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # acridine (C2v symmetry: 1=8, 2=7, 3=6, 4=5, 9 unique → 5 unique methyls)
    ("c1ccc2nc3ccccc3cc2c1",     "acridine"),
    ("Cc1cccc2cc3ccccc3nc12",    "1-methylacridine"),
    ("Cc1ccc2cc3ccccc3nc2c1",    "2-methylacridine"),
    ("Cc1ccc2nc3ccccc3cc2c1",    "3-methylacridine"),
    ("Cc1cccc2nc3ccccc3cc12",    "4-methylacridine"),
    ("Cc1c2ccccc2nc2ccccc12",    "9-methylacridine"),
    # phenanthridine (no symmetry; 10 unique positions 1-4, 6-10; pos 5 = N)
    ("c1ccc2c(c1)cnc1ccccc12",         "phenanthridine"),
    ("Cc1cccc2c1cnc1ccccc12",          "4-methylphenanthridine"),
    ("Cc1ccc2cnc3ccccc3c2c1",          "2-methylphenanthridine"),
    ("Cc1ccc2c(cnc3ccccc32)c1",        "3-methylphenanthridine"),
    ("Cc1cccc2cnc3ccccc3c12",          "1-methylphenanthridine"),
    ("Cc1nc2ccccc2c2ccccc12",          "6-methylphenanthridine"),
    ("Cc1cccc2c1ncc1ccccc12",          "7-methylphenanthridine"),
    ("Cc1ccc2c(c1)ncc1ccccc12",        "8-methylphenanthridine"),
    ("Cc1ccc2ncc3ccccc3c2c1",          "9-methylphenanthridine"),
    ("Cc1cccc2ncc3ccccc3c12",          "10-methylphenanthridine"),
    # phenazine (C2h symmetry: 1=6, 2=7, 3=8, 4=9 → 2 unique methyls; N at 5,10)
    ("c1ccc2nc3ccccc3nc2c1",    "phenazine"),
    ("Cc1cccc2nc3ccccc3nc12",   "1-methylphenazine"),
    ("Cc1ccc2nc3ccccc3nc2c1",   "2-methylphenazine"),
    # 9H-carbazole (C2v symmetry: 1=8, 2=7, 3=6, 4=5, 9H unique)
    ("c1ccc2c(c1)[nH]c1ccccc12",  "9H-carbazole"),
    ("Cc1cccc2c1[nH]c1ccccc12",   "1-methyl-9H-carbazole"),
    ("Cc1ccc2c(c1)[nH]c1ccccc12", "2-methyl-9H-carbazole"),
    ("Cc1ccc2[nH]c3ccccc3c2c1",   "3-methyl-9H-carbazole"),
    ("Cc1cccc2[nH]c3ccccc3c12",   "4-methyl-9H-carbazole"),
    ("Cn1c2ccccc2c2ccccc21",       "9-methyl-9H-carbazole"),
    # phenothiazine (C2v symmetry; N at 10, S bridging; 1=9, 2=8, 3=7, 4=6)
    ("c1ccc2c(c1)Nc1ccccc1S2",    "phenothiazine"),
    ("Cc1cccc2c1Nc1ccccc1S2",     "1-methylphenothiazine"),
    ("Cc1ccc2c(c1)Nc1ccccc1S2",   "2-methylphenothiazine"),
    ("Cc1ccc2c(c1)Sc1ccccc1N2",   "3-methylphenothiazine"),
    ("Cc1cccc2c1Sc1ccccc1N2",     "4-methylphenothiazine"),
    ("CN1c2ccccc2Sc2ccccc21",     "10-methylphenothiazine"),
    # phenoxazine (C2v symmetry; N at 10, O bridging; 1=9, 2=8, 3=7, 4=6)
    ("c1ccc2c(c1)Nc1ccccc1O2",    "phenoxazine"),
    ("Cc1cccc2c1Nc1ccccc1O2",     "1-methylphenoxazine"),
    ("Cc1ccc2c(c1)Nc1ccccc1O2",   "2-methylphenoxazine"),
    ("Cc1ccc2c(c1)Oc1ccccc1N2",   "3-methylphenoxazine"),
    ("Cc1cccc2c1Oc1ccccc1N2",     "4-methylphenoxazine"),
    ("CN1c2ccccc2Oc2ccccc21",     "10-methylphenoxazine"),
])
def test_phase690(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
