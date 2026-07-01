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
    ("O=c1occc2ccccc12",  "isochromen-1-one"),
    # 1,3-oxazolidine / 1,2-oxazolidine PINs (Phase 737)
    ("C1CNCO1",           "1,3-oxazolidine"),
    ("C1CCNO1",           "1,2-oxazolidine"),
    # 1,3-thiazolidine PIN (Phase 737)
    ("C1CSCN1",           "1,3-thiazolidine"),
    # lactam derivatives with updated locants (Phase 737)
    ("O=C1OCCN1",         "1,3-oxazolidin-2-one"),
    ("O=C1SCCN1",         "1,3-thiazolidin-2-one"),
    # regression: unrelated fused retained names unchanged
    ("C1=Cc2ccccc2OC1",   "2H-chromene"),
])
def test_phase263_coumarin_oxazolidine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
