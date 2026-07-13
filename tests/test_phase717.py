"""Phase 717: methyl derivatives of acridine, phenanthridine, phenazine,
phenothiazine, phenoxazine, phenoxathiin, dibenzofuran, dibenzothiophene,
perimidine, pteridine, 7H-purine, 9H-purine, and 1,2,4-benzotriazine.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # acridine (CH at 1,2,3,4,9; 5=4, 6=3, 7=2, 8=1 by symmetry)
    # Phase 854: corrected 1<->4, 2<->3 (locant map ring-swap fix)
    ("c1ccc2nc3ccccc3cc2c1",             "acridine"),
    ("Cc1cccc2cc3ccccc3nc12",            "4-methylacridine"),
    ("Cc1ccc2cc3ccccc3nc2c1",            "3-methylacridine"),
    ("Cc1ccc2nc3ccccc3cc2c1",            "2-methylacridine"),
    ("Cc1cccc2nc3ccccc3cc12",            "1-methylacridine"),
    ("Cc1c2ccccc2nc2ccccc12",            "9-methylacridine"),
    # phenanthridine (CH at 1,2,3,4,6,7,8,9,10)
    # Phase 854: corrected 1<->10, 2<->9, 3<->8, 4<->7 (locant map ring-swap fix)
    ("c1ccc2c(c1)cnc1ccccc12",           "phenanthridine"),
    ("Cc1cccc2cnc3ccccc3c12",            "10-methylphenanthridine"),
    ("Cc1ccc2cnc3ccccc3c2c1",            "9-methylphenanthridine"),
    ("Cc1ccc2c(cnc3ccccc32)c1",          "8-methylphenanthridine"),
    ("Cc1cccc2c1cnc1ccccc12",            "7-methylphenanthridine"),
    ("Cc1nc2ccccc2c2ccccc12",            "6-methylphenanthridine"),
    ("Cc1cccc2c1ncc1ccccc12",            "4-methylphenanthridine"),
    ("Cc1ccc2c(c1)ncc1ccccc12",          "3-methylphenanthridine"),
    ("Cc1ccc2ncc3ccccc3c2c1",            "2-methylphenanthridine"),
    ("Cc1cccc2ncc3ccccc3c12",            "1-methylphenanthridine"),
    # phenazine (CH at 1,2; 3=2, 4=1 by symmetry)
    ("c1ccc2nc3ccccc3nc2c1",             "phenazine"),
    ("Cc1cccc2nc3ccccc3nc12",            "1-methylphenazine"),
    ("Cc1ccc2nc3ccccc3nc2c1",            "2-methylphenazine"),
    # phenothiazine (CH at 1,2,3,4; 6=4, 7=3, 8=2, 9=1 by symmetry)
    ("c1ccc2c(c1)Nc1ccccc1S2",           "10H-phenothiazine"),
    ("Cc1cccc2c1Nc1ccccc1S2",            "1-methyl-10H-phenothiazine"),
    ("Cc1ccc2c(c1)Nc1ccccc1S2",          "2-methyl-10H-phenothiazine"),
    ("Cc1ccc2c(c1)Sc1ccccc1N2",          "3-methyl-10H-phenothiazine"),
    ("Cc1cccc2c1Sc1ccccc1N2",            "4-methyl-10H-phenothiazine"),
    # phenoxazine (CH at 1,2,3,4; 6=4, 7=3, 8=2, 9=1 by symmetry)
    ("c1ccc2c(c1)Nc1ccccc1O2",           "10H-phenoxazine"),
    ("Cc1cccc2c1Nc1ccccc1O2",            "1-methyl-10H-phenoxazine"),
    ("Cc1ccc2c(c1)Nc1ccccc1O2",          "2-methyl-10H-phenoxazine"),
    ("Cc1ccc2c(c1)Oc1ccccc1N2",          "3-methyl-10H-phenoxazine"),
    ("Cc1cccc2c1Oc1ccccc1N2",            "4-methyl-10H-phenoxazine"),
    # phenoxathiin (CH at 1,2,3,4; no symmetry — O and S differ)
    ("c1ccc2c(c1)Sc1ccccc1O2",           "phenoxathiin"),
    ("Cc1cccc2c1Sc1ccccc1O2",            "1-methylphenoxathiin"),
    ("Cc1ccc2c(c1)Sc1ccccc1O2",          "2-methylphenoxathiin"),
    ("Cc1ccc2c(c1)Oc1ccccc1S2",          "3-methylphenoxathiin"),
    ("Cc1cccc2c1Oc1ccccc1S2",            "4-methylphenoxathiin"),
    # dibenzofuran (CH at 1,2,3,4; 6=4, 7=3, 8=2, 9=1 by symmetry)
    ("c1ccc2c(c1)oc1ccccc12",            "dibenzofuran"),
    ("Cc1cccc2oc3ccccc3c12",             "1-methyldibenzofuran"),
    ("Cc1ccc2oc3ccccc3c2c1",             "2-methyldibenzofuran"),
    ("Cc1ccc2c(c1)oc1ccccc12",           "3-methyldibenzofuran"),
    ("Cc1cccc2c1oc1ccccc12",             "4-methyldibenzofuran"),
    # dibenzothiophene (CH at 1,2,3,4; 6=4, 7=3, 8=2, 9=1 by symmetry)
    ("c1ccc2c(c1)sc1ccccc12",            "dibenzothiophene"),
    ("Cc1cccc2sc3ccccc3c12",             "1-methyldibenzothiophene"),
    ("Cc1ccc2sc3ccccc3c2c1",             "2-methyldibenzothiophene"),
    ("Cc1ccc2c(c1)sc1ccccc12",           "3-methyldibenzothiophene"),
    ("Cc1cccc2c1sc1ccccc12",             "4-methyldibenzothiophene"),
    # perimidine (CH at 2,4,5,6; with tautomers at 4/5/6)
    ("C1=Nc2cccc3cccc(c23)N1",           "perimidine"),
    ("CC1=Nc2cccc3cccc(c23)N1",          "2-methylperimidine"),
    ("Cc1ccc2cccc3c2c1N=CN3",            "4-methylperimidine"),
    ("Cc1ccc2cccc3c2c1NC=N3",            "4-methylperimidine"),
    ("Cc1cc2c3c(cccc3c1)NC=N2",          "5-methylperimidine"),
    ("Cc1cc2c3c(cccc3c1)N=CN2",          "5-methylperimidine"),
    ("Cc1ccc2c3c(cccc13)NC=N2",          "6-methylperimidine"),
    ("Cc1ccc2c3c(cccc13)N=CN2",          "6-methylperimidine"),
    # pteridine (CH at 2,4,6,7)
    ("c1cnc2ncncc2n1",                   "pteridine"),
    ("Cc1ncc2nccnc2n1",                  "2-methylpteridine"),
    ("Cc1ncnc2nccnc12",                  "4-methylpteridine"),
    ("Cc1cnc2ncncc2n1",                  "6-methylpteridine"),
    ("Cc1cnc2cncnc2n1",                  "7-methylpteridine"),
    # 7H-purine (CH at 2,6,8)
    ("c1ncc2[nH]cnc2n1",                 "7H-purine"),
    ("Cc1ncc2[nH]cnc2n1",               "2-methyl-7H-purine"),
    ("Cc1ncnc2nc[nH]c12",               "6-methyl-7H-purine"),
    ("Cc1nc2ncncc2[nH]1",               "8-methyl-7H-purine"),
    # 9H-purine (CH at 2,6,8)
    ("c1ncc2nc[nH]c2n1",                 "9H-purine"),
    ("Cc1ncc2nc[nH]c2n1",               "2-methyl-9H-purine"),
    ("Cc1ncnc2[nH]cnc12",               "6-methyl-9H-purine"),
    ("Cc1nc2cncnc2[nH]1",               "8-methyl-9H-purine"),
    # 1,2,4-benzotriazine (CH at 3,5,6,7,8)
    ("c1ccc2nncnc2c1",                   "1,2,4-benzotriazine"),
    ("Cc1nnc2ccccc2n1",                  "3-methyl-1,2,4-benzotriazine"),
    ("Cc1cccc2nncnc12",                  "5-methyl-1,2,4-benzotriazine"),
    ("Cc1ccc2nncnc2c1",                  "6-methyl-1,2,4-benzotriazine"),
    ("Cc1ccc2ncnnc2c1",                  "7-methyl-1,2,4-benzotriazine"),
    ("Cc1cccc2ncnnc12",                  "8-methyl-1,2,4-benzotriazine"),
])
def test_phase717(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
