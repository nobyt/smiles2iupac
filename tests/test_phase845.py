"""Phase 845: drop nH- prefix indicated-H when the N at that locant is substituted.

IUPAC 2013 P-14.7.1: when the atom bearing the indicated hydrogen is substituted
by a substituent group, the indicated hydrogen is omitted.  This extends Phase 844
(which covered inline (nH) tokens) to the prefix form used by fused heteroaromatics:
1H-indole, 1H-benzimidazole, 1H-benzotriazole, 1H-indazole, etc.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-indole: N1-H → keep 1H-
    ("c1ccc2[nH]ccc2c1",     "1H-indole"),
    # N1-methyl → drop 1H-
    ("Cn1ccc2ccccc21",        "1-methylindole"),
    # C5-methyl: N1 still has H → keep 1H-
    ("Cc1ccc2[nH]ccc2c1",    "5-methyl-1H-indole"),
    # 1H-benzimidazole: N1-H → keep 1H-
    ("c1ccc2[nH]cnc2c1",     "1H-benzimidazole"),
    # N1-methyl → drop 1H-
    ("Cn1cnc2ccccc21",        "1-methylbenzimidazole"),
    # C5-methyl: N1 still has H → keep 1H-
    ("Cc1ccc2[nH]cnc2c1",    "5-methyl-1H-benzimidazole"),
    # 1H-benzotriazole: N1-H → keep 1H-
    ("c1ccc2[nH]nnc2c1",     "1H-benzotriazole"),
    # N1-methyl → drop 1H-
    ("Cn1nnc2ccccc21",        "1-methylbenzotriazole"),
    # C5-methyl: N1 still has H → keep 1H-
    ("Cc1ccc2[nH]nnc2c1",    "5-methyl-1H-benzotriazole"),
    # 1H-indazole: N1-H → keep 1H-
    ("c1ccc2[nH]ncc2c1",     "1H-indazole"),
    # N1-methyl → drop 1H-
    ("Cn1ncc2ccccc21",        "1-methylindazole"),
    # C5-methyl: N1 still has H → keep 1H-
    ("Cc1ccc2[nH]ncc2c1",    "5-methyl-1H-indazole"),
    # isatin (1H-indole-2,3-dione): N1-methyl → drop 1H-
    ("O=C1N(C)c2ccccc2C1=O", "1-methylindole-2,3-dione"),
    # isatin: C5-methyl → keep 1H-
    ("O=C1Nc2cc(C)ccc2C1=O", "5-methyl-1H-indole-2,3-dione"),
])
def test_phase845(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
