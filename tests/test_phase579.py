"""Phase 579: Substituted 1H-naphtho[2,3-d]imidazole, 1H-naphtho[2,3-d]pyrazole,
1H-naphtho[2,1-d]imidazole, 1H-naphtho[2,1-d]pyrazole, naphtho[1,2-d]oxazole,
and naphtho[1,2-d]thiazole naming.
Imidazoles: N-H at pos 1, N at pos 3; sub C: 2, 4-9.
Pyrazoles: N-H at pos 1, N at pos 2; sub C: 3-9.
Oxazole/thiazole: O/S at pos 1 (→ None), N at pos 3; sub C: 2, 4-9.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-naphtho[2,3-d]imidazole
    ("c1ccc2cc3[nH]cnc3cc2c1",          "1H-naphtho[2,3-d]imidazole"),
    ("Cc1nc2cc3ccccc3cc2[nH]1",         "2-methyl-1H-naphtho[2,3-d]imidazole"),
    ("Cc1c2ccccc2cc2[nH]cnc12",         "4-methyl-1H-naphtho[2,3-d]imidazole"),
    ("Cc1cccc2cc3[nH]cnc3cc12",         "5-methyl-1H-naphtho[2,3-d]imidazole"),
    ("Cc1ccc2cc3[nH]cnc3cc2c1",         "6-methyl-1H-naphtho[2,3-d]imidazole"),
    ("Cc1ccc2cc3nc[nH]c3cc2c1",         "7-methyl-1H-naphtho[2,3-d]imidazole"),
    ("Cc1cccc2cc3nc[nH]c3cc12",         "8-methyl-1H-naphtho[2,3-d]imidazole"),
    ("Cc1c2ccccc2cc2nc[nH]c12",         "9-methyl-1H-naphtho[2,3-d]imidazole"),
    # 1H-naphtho[2,3-d]pyrazole
    ("c1ccc2cc3[nH]ncc3cc2c1",          "1H-naphtho[2,3-d]pyrazole"),
    ("Cc1n[nH]c2cc3ccccc3cc12",         "3-methyl-1H-naphtho[2,3-d]pyrazole"),
    ("Cc1c2ccccc2cc2[nH]ncc12",         "4-methyl-1H-naphtho[2,3-d]pyrazole"),
    ("Cc1cccc2cc3[nH]ncc3cc12",         "5-methyl-1H-naphtho[2,3-d]pyrazole"),
    ("Cc1ccc2cc3[nH]ncc3cc2c1",         "6-methyl-1H-naphtho[2,3-d]pyrazole"),
    ("Cc1ccc2cc3cn[nH]c3cc2c1",         "7-methyl-1H-naphtho[2,3-d]pyrazole"),
    ("Cc1cccc2cc3cn[nH]c3cc12",         "8-methyl-1H-naphtho[2,3-d]pyrazole"),
    ("Cc1c2ccccc2cc2cn[nH]c12",         "9-methyl-1H-naphtho[2,3-d]pyrazole"),
    # 1H-naphtho[2,1-d]imidazole
    ("c1ccc2c(c1)ccc1[nH]cnc12",        "1H-naphtho[2,1-d]imidazole"),
    ("Cc1nc2ccc3ccccc3c2[nH]1",         "2-methyl-1H-naphtho[2,1-d]imidazole"),
    ("Cc1cc2ccccc2c2[nH]cnc12",         "4-methyl-1H-naphtho[2,1-d]imidazole"),
    ("Cc1cc2nc[nH]c2c2ccccc12",         "5-methyl-1H-naphtho[2,1-d]imidazole"),
    ("Cc1cccc2c1ccc1nc[nH]c12",         "6-methyl-1H-naphtho[2,1-d]imidazole"),
    ("Cc1ccc2c(ccc3nc[nH]c32)c1",       "7-methyl-1H-naphtho[2,1-d]imidazole"),
    ("Cc1ccc2ccc3nc[nH]c3c2c1",         "8-methyl-1H-naphtho[2,1-d]imidazole"),
    ("Cc1cccc2ccc3nc[nH]c3c12",         "9-methyl-1H-naphtho[2,1-d]imidazole"),
    # 1H-naphtho[2,1-d]pyrazole
    ("c1ccc2c(c1)ccc1cn[nH]c12",        "1H-naphtho[2,1-d]pyrazole"),
    ("Cc1n[nH]c2c1ccc1ccccc12",         "3-methyl-1H-naphtho[2,1-d]pyrazole"),
    ("Cc1cc2ccccc2c2[nH]ncc12",         "4-methyl-1H-naphtho[2,1-d]pyrazole"),
    ("Cc1cc2cn[nH]c2c2ccccc12",         "5-methyl-1H-naphtho[2,1-d]pyrazole"),
    ("Cc1cccc2c1ccc1cn[nH]c12",         "6-methyl-1H-naphtho[2,1-d]pyrazole"),
    ("Cc1ccc2c(ccc3cn[nH]c32)c1",       "7-methyl-1H-naphtho[2,1-d]pyrazole"),
    ("Cc1ccc2ccc3cn[nH]c3c2c1",         "8-methyl-1H-naphtho[2,1-d]pyrazole"),
    ("Cc1cccc2ccc3cn[nH]c3c12",         "9-methyl-1H-naphtho[2,1-d]pyrazole"),
    # naphtho[1,2-d]oxazole
    ("c1ccc2c(c1)ccc1ocnc12",           "naphtho[1,2-d]oxazole"),
    ("Cc1nc2c(ccc3ccccc32)o1",          "2-methylnaphtho[1,2-d]oxazole"),
    ("Cc1cc2ccccc2c2ncoc12",            "4-methylnaphtho[1,2-d]oxazole"),
    ("Cc1cc2ocnc2c2ccccc12",            "5-methylnaphtho[1,2-d]oxazole"),
    ("Cc1cccc2c1ccc1ocnc12",            "6-methylnaphtho[1,2-d]oxazole"),
    ("Cc1ccc2c(ccc3ocnc32)c1",          "7-methylnaphtho[1,2-d]oxazole"),
    ("Cc1ccc2ccc3ocnc3c2c1",            "8-methylnaphtho[1,2-d]oxazole"),
    ("Cc1cccc2ccc3ocnc3c12",            "9-methylnaphtho[1,2-d]oxazole"),
    # naphtho[1,2-d]thiazole
    ("c1ccc2c(c1)ccc1scnc12",           "naphtho[1,2-d]thiazole"),
    ("Cc1nc2c(ccc3ccccc32)s1",          "2-methylnaphtho[1,2-d]thiazole"),
    ("Cc1cc2ccccc2c2ncsc12",            "4-methylnaphtho[1,2-d]thiazole"),
    ("Cc1cc2scnc2c2ccccc12",            "5-methylnaphtho[1,2-d]thiazole"),
    ("Cc1cccc2c1ccc1scnc12",            "6-methylnaphtho[1,2-d]thiazole"),
    ("Cc1ccc2c(ccc3scnc32)c1",          "7-methylnaphtho[1,2-d]thiazole"),
    ("Cc1ccc2ccc3scnc3c2c1",            "8-methylnaphtho[1,2-d]thiazole"),
    ("Cc1cccc2ccc3scnc3c12",            "9-methylnaphtho[1,2-d]thiazole"),
])
def test_phase579_naphtho_n_heterocycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
