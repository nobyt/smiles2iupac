"""Phase 134: 縮合環保留名追加 (IUPAC 2013 P-31.1.3)

fluorene, xanthene, thioxanthene, chromen-2-one (coumarin), isocoumarin,
2H-chromene, 4H-chromene, phenoxazine, phenothiazine, xanthen-9-one, perimidine
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # fluorene (dibenzo[a,d][7]annulene with CH2)
    ("C1c2ccccc2-c2ccccc21", "fluorene"),
    # xanthene (9H-xanthene)
    ("C1c2ccccc2Oc2ccccc21", "xanthene"),
    # thioxanthene
    ("C1c2ccccc2Sc2ccccc21", "thioxanthene"),
    # coumarin (retained PIN, IUPAC 2013 P-31.1.3.4)
    ("O=c1ccc2ccccc2o1", "coumarin"),
    # isocoumarin (retained PIN, IUPAC 2013 P-31.1.3.4)
    ("O=c1occc2ccccc12", "isocoumarin"),
    # 2H-chromene
    ("C1=COc2ccccc2C1", "2H-chromene"),
    # 4H-chromene
    ("C1=Cc2ccccc2OC1", "4H-chromene"),
    # phenoxazine
    ("c1ccc2c(c1)Nc1ccccc1O2", "phenoxazine"),
    # phenothiazine
    ("c1ccc2c(c1)Nc1ccccc1S2", "phenothiazine"),
    # xanthen-9-one (xanthone)
    ("O=C1c2ccccc2Oc2ccccc21", "xanthen-9-one"),
    # perimidine (2,3-dihydro-1H-perimidine parent skeleton)
    ("C1=Nc2cccc3cccc1c23", "perimidine"),
    # 回帰: Phase 133 partially saturated fused compounds unchanged
    ("C1CCc2ccccc21", "indane"),
    ("C1COc2ccccc2C1", "chromane"),
    # 回帰: Phase 131–132 fused heteroaromatics unchanged
    ("c1ccc2c(c1)[nH]c1ccccc12", "9H-carbazole"),
    ("c1cnc2ncccc2c1", "1,8-naphthyridine"),
])
def test_phase134_more_fused_retained(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
