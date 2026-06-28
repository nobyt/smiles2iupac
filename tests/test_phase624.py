"""Phase 624: methyl-substituted tricyclic compounds (IUPAC 2013).

Xanthene, thioxanthene, 9H-carbazole, phenoxazine, phenothiazine,
xanthen-9-one, dibenzofuran, dibenzothiophene — all carbon positions.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # xanthene (C9 sp3)
    ("CC1c2ccccc2Oc2ccccc21",     "9-methyl-9H-xanthene"),
    ("Cc1cccc2c1Cc1ccccc1O2",     "1-methyl-9H-xanthene"),
    ("Cc1ccc2c(c1)Cc1ccccc1O2",   "2-methyl-9H-xanthene"),
    ("Cc1ccc2c(c1)Oc1ccccc1C2",   "3-methyl-9H-xanthene"),
    ("Cc1cccc2c1Oc1ccccc1C2",     "4-methyl-9H-xanthene"),
    # thioxanthene (C9 sp3)
    ("CC1c2ccccc2Sc2ccccc21",     "9-methyl-9H-thioxanthene"),
    ("Cc1cccc2c1Cc1ccccc1S2",     "1-methyl-9H-thioxanthene"),
    ("Cc1ccc2c(c1)Cc1ccccc1S2",   "2-methyl-9H-thioxanthene"),
    ("Cc1ccc2c(c1)Sc1ccccc1C2",   "3-methyl-9H-thioxanthene"),
    ("Cc1cccc2c1Sc1ccccc1C2",     "4-methyl-9H-thioxanthene"),
    # 9H-carbazole (NH at 9)
    ("Cc1ccc2c(c1)[nH]c1ccccc12", "2-methyl-9H-carbazole"),
    ("Cc1ccc2[nH]c3ccccc3c2c1",   "3-methyl-9H-carbazole"),
    ("Cc1cccc2[nH]c3ccccc3c12",   "4-methyl-9H-carbazole"),
    ("Cc1cccc2c1[nH]c1ccccc12",   "1-methyl-9H-carbazole"),
    # phenoxazine (O and NH)
    ("Cc1ccc2c(c1)Nc1ccccc1O2",   "2-methyl-10H-phenoxazine"),
    ("Cc1ccc2c(c1)Oc1ccccc1N2",   "3-methyl-10H-phenoxazine"),
    ("Cc1cccc2c1Oc1ccccc1N2",     "4-methyl-10H-phenoxazine"),
    ("Cc1cccc2c1Nc1ccccc1O2",     "1-methyl-10H-phenoxazine"),
    # phenothiazine (S and NH)
    ("Cc1ccc2c(c1)Nc1ccccc1S2",   "2-methyl-10H-phenothiazine"),
    ("Cc1ccc2c(c1)Sc1ccccc1N2",   "3-methyl-10H-phenothiazine"),
    ("Cc1cccc2c1Sc1ccccc1N2",     "4-methyl-10H-phenothiazine"),
    ("Cc1cccc2c1Nc1ccccc1S2",     "1-methyl-10H-phenothiazine"),
    # xanthen-9-one (C=O at 9)
    ("Cc1cccc2oc3ccccc3c(=O)c12", "1-methylxanthen-9-one"),
    ("Cc1ccc2oc3ccccc3c(=O)c2c1", "2-methylxanthen-9-one"),
    ("Cc1ccc2c(=O)c3ccccc3oc2c1", "3-methylxanthen-9-one"),
    ("Cc1cccc2c(=O)c3ccccc3oc12", "4-methylxanthen-9-one"),
    # dibenzofuran
    ("Cc1ccc2c(c1)oc1ccccc12",    "3-methyldibenzofuran"),
    ("Cc1ccc2oc3ccccc3c2c1",      "2-methyldibenzofuran"),
    ("Cc1cccc2oc3ccccc3c12",      "1-methyldibenzofuran"),
    ("Cc1cccc2c1oc1ccccc12",      "4-methyldibenzofuran"),
    # dibenzothiophene
    ("Cc1ccc2c(c1)sc1ccccc12",    "3-methyldibenzothiophene"),
    ("Cc1ccc2sc3ccccc3c2c1",      "2-methyldibenzothiophene"),
    ("Cc1cccc2sc3ccccc3c12",      "1-methyldibenzothiophene"),
    ("Cc1cccc2c1sc1ccccc12",      "4-methyldibenzothiophene"),
])
def test_phase624_methyl_tricyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
