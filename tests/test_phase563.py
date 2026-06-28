"""Phase 563: Substituted dibenzofuran, dibenzothiophene, fluorene, xanthene,
thioxanthene, phenoxazine, phenothiazine, phenoxathiin, and thianthrene naming.
Dibenzofuran/dibenzothiophene: C2-symmetric, positions 1-4 and 6-9 (pairs sum to 10).
Fluorene: C2-symmetric, positions 1-8 + C9 (pairs sum to 9).
Xanthene/thioxanthene: C2-symmetric, positions 1-9 (pairs sum to 9).
Phenoxazine/phenothiazine: C2-symmetric, positions 1-4,6-9 + N10 (pairs sum to 10).
Phenoxathiin: C2-symmetric, positions 1-4 and 6-9 (pairs sum to 10).
Thianthrene: D2h-symmetric, only 2 unique C environments (locants 1 and 2).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # dibenzofuran (C2-symmetric; 4 unique C positions: 1,2,3,4)
    ("c1ccc2c(c1)oc1ccccc12",   "dibenzofuran"),
    ("Cc1cccc2oc3ccccc3c12",    "1-methyldibenzofuran"),
    ("Cc1ccc2oc3ccccc3c2c1",    "2-methyldibenzofuran"),
    ("Cc1ccc2c(c1)oc1ccccc12",  "3-methyldibenzofuran"),
    ("Cc1cccc2c1oc1ccccc12",    "4-methyldibenzofuran"),
    # dibenzothiophene (C2-symmetric; 4 unique C positions: 1,2,3,4)
    ("c1ccc2c(c1)sc1ccccc12",   "dibenzothiophene"),
    ("Cc1cccc2sc3ccccc3c12",    "1-methyldibenzothiophene"),
    ("Cc1ccc2sc3ccccc3c2c1",    "2-methyldibenzothiophene"),
    ("Cc1ccc2c(c1)sc1ccccc12",  "3-methyldibenzothiophene"),
    ("Cc1cccc2c1sc1ccccc12",    "4-methyldibenzothiophene"),
    # fluorene (C2-symmetric; 5 unique positions: 1,2,3,4,9)
    ("c1ccc2c(c1)Cc1ccccc1-2",      "9H-fluorene"),
    ("Cc1cccc2c1Cc1ccccc1-2",        "1-methyl-9H-fluorene"),
    ("Cc1ccc2c(c1)Cc1ccccc1-2",      "2-methyl-9H-fluorene"),
    ("Cc1ccc2c(c1)-c1ccccc1C2",      "3-methyl-9H-fluorene"),
    ("Cc1cccc2c1-c1ccccc1C2",        "4-methyl-9H-fluorene"),
    ("CC1c2ccccc2-c2ccccc21",        "9-methyl-9H-fluorene"),
    # xanthene (C2-symmetric; 5 unique positions: 1,2,3,4,9)
    ("c1ccc2c(c1)Cc1ccccc1O2",       "9H-xanthene"),
    ("Cc1cccc2c1Cc1ccccc1O2",        "1-methyl-9H-xanthene"),
    ("Cc1ccc2c(c1)Cc1ccccc1O2",      "2-methyl-9H-xanthene"),
    ("Cc1ccc2c(c1)Oc1ccccc1C2",      "3-methyl-9H-xanthene"),
    ("Cc1cccc2c1Oc1ccccc1C2",        "4-methyl-9H-xanthene"),
    ("CC1c2ccccc2Oc2ccccc21",        "9-methyl-9H-xanthene"),
    # thioxanthene (C2-symmetric; 5 unique positions: 1,2,3,4,9)
    ("c1ccc2c(c1)Cc1ccccc1S2",       "9H-thioxanthene"),
    ("Cc1cccc2c1Cc1ccccc1S2",        "1-methyl-9H-thioxanthene"),
    ("Cc1ccc2c(c1)Cc1ccccc1S2",      "2-methyl-9H-thioxanthene"),
    ("Cc1ccc2c(c1)Sc1ccccc1C2",      "3-methyl-9H-thioxanthene"),
    ("Cc1cccc2c1Sc1ccccc1C2",        "4-methyl-9H-thioxanthene"),
    ("CC1c2ccccc2Sc2ccccc21",        "9-methyl-9H-thioxanthene"),
    # phenoxazine (C2-symmetric; 5 unique positions: 1,2,3,4,10)
    ("c1ccc2c(c1)Nc1ccccc1O2",       "10H-phenoxazine"),
    ("Cc1cccc2c1Nc1ccccc1O2",        "1-methyl-10H-phenoxazine"),
    ("Cc1ccc2c(c1)Nc1ccccc1O2",      "2-methyl-10H-phenoxazine"),
    ("Cc1ccc2c(c1)Oc1ccccc1N2",      "3-methyl-10H-phenoxazine"),
    ("Cc1cccc2c1Oc1ccccc1N2",        "4-methyl-10H-phenoxazine"),
    ("CN1c2ccccc2Oc2ccccc21",        "10-methyl-10H-phenoxazine"),
    # phenothiazine (C2-symmetric; 5 unique positions: 1,2,3,4,10)
    ("c1ccc2c(c1)Nc1ccccc1S2",       "10H-phenothiazine"),
    ("Cc1cccc2c1Nc1ccccc1S2",        "1-methyl-10H-phenothiazine"),
    ("Cc1ccc2c(c1)Nc1ccccc1S2",      "2-methyl-10H-phenothiazine"),
    ("Cc1ccc2c(c1)Sc1ccccc1N2",      "3-methyl-10H-phenothiazine"),
    ("Cc1cccc2c1Sc1ccccc1N2",        "4-methyl-10H-phenothiazine"),
    ("CN1c2ccccc2Sc2ccccc21",        "10-methyl-10H-phenothiazine"),
    # phenoxathiin (C2-symmetric; 4 unique positions: 1,2,3,4)
    ("c1ccc2c(c1)Oc1ccccc1S2",       "phenoxathiin"),
    ("Cc1cccc2c1Sc1ccccc1O2",        "1-methylphenoxathiin"),
    ("Cc1ccc2c(c1)Sc1ccccc1O2",      "2-methylphenoxathiin"),
    ("Cc1ccc2c(c1)Oc1ccccc1S2",      "3-methylphenoxathiin"),
    ("Cc1cccc2c1Oc1ccccc1S2",        "4-methylphenoxathiin"),
    # thianthrene (D2h-symmetric; 2 unique C environments: locants 1 and 2)
    ("c1ccc2c(c1)Sc1ccccc1S2",       "thianthrene"),
    ("Cc1cccc2c1Sc1ccccc1S2",        "1-methylthianthrene"),
    ("Cc1ccc2c(c1)Sc1ccccc1S2",      "2-methylthianthrene"),
])
def test_phase563_fused_hetero(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
