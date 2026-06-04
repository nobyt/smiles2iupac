"""Phase 418: 1H-Indole-2,3-dione (isatin), benzofuran-2(3H)-one,
and benzo[b]thiophen-2(3H)-one.

IUPAC 2013 P-31.1.3: retained/systematic names for benzo-fused 5-membered
keto compounds with O or S heteroatoms.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-indole-2,3-dione (isatin) — two exo C=O, N1-H
    ("O=C1Nc2ccccc2C1=O",            "1H-indole-2,3-dione"),
    # benzofuran-2(3H)-one — O in ring, C2=O, C3 is CH2
    ("O=C1COc2ccccc21",              "benzofuran-2(3H)-one"),
    # benzo[b]thiophen-2(3H)-one — S in ring, C2=O, C3 is CH2
    ("O=C1CSc2ccccc21",              "benzo[b]thiophen-2(3H)-one"),
    # regression: 1H-indole unchanged
    ("c1ccc2[nH]ccc2c1",              "1H-indole"),
    # regression: benzofuran unchanged
    ("c1ccc2occc2c1",                  "benzofuran"),
    # regression: benzo[b]thiophene unchanged
    ("c1ccc2sccc2c1",                  "benzo[b]thiophene"),
    # regression: indan-2-one unchanged (Phase 411)
    ("O=C1Cc2ccccc2C1",               "indan-2-one"),
    # regression: benzene unchanged
    ("c1ccccc1",                        "benzene"),
])
def test_phase418_benzo_keto_heterocycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
