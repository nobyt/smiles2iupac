"""Phase 131: 三環縮合ヘテロ芳香族保留名 (IUPAC 2013 P-31.1.3)

carbazole, dibenzothiophene, dibenzofuran
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 9H-carbazole (dibenzo[b,d]pyrrole)
    ("c1ccc2c(c1)[nH]c1ccccc12", "9H-carbazole"),
    # dibenzothiophene
    ("c1ccc2c(c1)sc1ccccc12", "dibenzothiophene"),
    # dibenzofuran
    ("c1ccc2c(c1)oc1ccccc12", "dibenzofuran"),
    # 回帰: Phase 130 fused compounds unchanged
    ("c1ccc2ocnc2c1", "1,3-benzoxazole"),
    ("c1ccc2scnc2c1", "1,3-benzothiazole"),
    ("c1ccc2nc3ccccc3nc2c1", "phenazine"),
    ("c1cnc2ncncc2n1", "pteridine"),
    # 回帰: Phase 17 fused compounds unchanged
    ("c1ccc2[nH]ccc2c1", "1H-indole"),
    ("c1ccc2ncccc2c1", "quinoline"),
])
def test_phase131_tricyclic_and_naphthyridine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
