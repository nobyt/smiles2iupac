"""Phase 741: chromane → 3,4-dihydro-2H-chromene and
coumarin → chromen-2-one (IUPAC 2013 PINs).

Retained names chromane/isochromane/coumarin are not PINs;
the systematic names are used instead.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # chromane → 3,4-dihydro-2H-chromene
    ("c1ccc2c(c1)CCCO2",              "3,4-dihydro-2H-chromene"),
    ("CC1CCc2ccccc2O1",               "2-methyl-3,4-dihydro-2H-chromene"),
    ("CC1COc2ccccc2C1",               "3-methyl-3,4-dihydro-2H-chromene"),
    ("CC1CCOc2ccccc21",               "4-methyl-3,4-dihydro-2H-chromene"),
    ("Cc1ccc2c(c1)OCCC2",             "7-methyl-3,4-dihydro-2H-chromene"),
    # isochromane → 3,4-dihydro-1H-isochromene
    ("c1ccc2c(c1)CCOC2",              "3,4-dihydro-1H-isochromene"),
    ("CC1OCCc2ccccc21",               "1-methyl-3,4-dihydro-1H-isochromene"),
    ("CC1COCc2ccccc21",               "4-methyl-3,4-dihydro-1H-isochromene"),
    # coumarin → chromen-2-one
    ("O=c1ccc2ccccc2o1",              "chromen-2-one"),
    ("Cc1cc2ccccc2oc1=O",             "3-methylchromen-2-one"),
    ("Cc1cc(=O)oc2ccccc12",           "4-methylchromen-2-one"),
    ("Cc1ccc2ccc(=O)oc2c1",           "7-methylchromen-2-one"),
    ("O=c1cc(O)c2ccccc2o1",           "4-hydroxychromen-2-one"),
    ("O=c1cc(N)c2ccccc2o1",           "4-aminochromen-2-one"),
    ("O=c1ccc2ccc(O)cc2o1",           "7-hydroxychromen-2-one"),
    # Regression: existing PINs unchanged
    ("O=c1ccoc2ccccc12",              "chromone"),            # chromone is a PIN
    ("O=c1occc2ccccc12",              "isochromen-1-one"),    # isocoumarin
    ("C1=Cc2ccccc2OC1",               "2H-chromene"),
    ("C1=COc2ccccc2C1",               "4H-chromene"),
])
def test_phase741_chromane_coumarin_pins(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
