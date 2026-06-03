"""Phase 263: coumarin / isocoumarin retained names; oxazolidine / thiazolidine
without locant prefix (IUPAC 2013 P-31.1.3.4).

  O=c1ccc2ccccc2o1  → coumarin      (was 'chromen-2-one')
  O=c1occc2ccccc12  → isocoumarin   (was '1H-2-benzopyran-1-one')
  C1CNCO1           → oxazolidine   (was '1,3-oxazolidine')
  C1CSCN1           → thiazolidine  (was '1,3-thiazolidine')
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # coumarin / isocoumarin (retained PINs, IUPAC 2013 P-31.1.3.4)
    ("O=c1ccc2ccccc2o1",  "coumarin"),
    ("O=c1occc2ccccc12",  "isocoumarin"),
    # oxazolidine / isoxazolidine (retained names without locant prefix)
    ("C1CNCO1",           "oxazolidine"),
    ("C1CCNO1",           "isoxazolidine"),
    # thiazolidine / isothiazolidine
    ("C1CSCN1",           "thiazolidine"),
    # lactam/lactone derivatives preserve systematic locant-2-one suffix
    ("O=C1OCCN1",         "oxazolidin-2-one"),
    ("O=C1SCCN1",         "thiazolidin-2-one"),
    # regression: unrelated fused retained names unchanged
    ("C1=COc2ccccc2C1",   "2H-chromene"),
])
def test_phase263_coumarin_oxazolidine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
