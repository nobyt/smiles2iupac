"""Phase 137: 核酸塩基 保留名 (IUPAC 2013 P-14.5)

adenine, guanine, cytosine, uracil, thymine, hypoxanthine, xanthine, 9H-purine
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # adenine (6-aminopurine)
    ("Nc1ncnc2[nH]cnc12", "adenine"),
    # guanine (2-amino-1H-purin-6(9H)-one)
    ("Nc1nc2[nH]cnc2c(=O)[nH]1", "guanine"),
    # cytosine (4-aminopyrimidin-2(1H)-one)
    ("Nc1ccnc(=O)[nH]1", "cytosine"),
    # uracil: IUPAC 2013 preferred = pyrimidine-2,4(1H,3H)-dione (Phase 401)
    ("O=c1cc[nH]c(=O)[nH]1", "pyrimidine-2,4(1H,3H)-dione"),
    # thymine: IUPAC 2013 preferred = 5-methylpyrimidine-2,4(1H,3H)-dione (Phase 401)
    ("Cc1c[nH]c(=O)[nH]c1=O", "5-methylpyrimidine-2,4(1H,3H)-dione"),
    # hypoxanthine (3,9-dihydro-1H-purin-6(2H)-one)
    ("O=c1[nH]cnc2[nH]cnc12", "hypoxanthine"),
    # xanthine (3,7-dihydro-1H-purine-2,6-dione)
    ("O=c1[nH]c(=O)c2[nH]cnc2[nH]1", "xanthine"),
    # 7H-purine (NH at N7, adjacent to C5)
    ("c1ncc2[nH]cnc2n1", "7H-purine"),
    # 回帰: fused heteroaromatics unchanged (Phase 131-134)
    ("c1ccc2c(c1)[nH]c1ccccc12", "9H-carbazole"),
    ("c1cnc2ncccc2c1",           "1,8-naphthyridine"),
    # 回帰: amino acids unchanged (Phase 135)
    ("NCC(=O)O", "glycine"),
    # 回帰: benzene unchanged
    ("c1ccccc1", "benzene"),
])
def test_phase137_nucleobases(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
