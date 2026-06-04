"""Phase 417: 1H-Indazol-3(2H)-one and phenanthridin-6(5H)-one.

IUPAC 2013 P-31.1.3: retained/systematic names for benzo-fused pyrazolone
(5-membered ring) and the tricyclic phenanthridinone.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-indazol-3(2H)-one (C3=O, N1-H, N2-H)
    ("O=c1[nH][nH]c2ccccc12",            "1H-indazol-3(2H)-one"),
    # phenanthridin-6(5H)-one (C6=O, N5-H, tricyclic)
    ("O=c1[nH]c2ccccc2c2ccccc12",        "phenanthridin-6(5H)-one"),
    # regression: 1H-indazole unchanged
    ("c1ccc2[nH]ncc2c1",                  "1H-indazole"),
    # regression: phenanthridine unchanged (canonical: c1ccc2c(c1)cnc1ccccc12)
    ("c1nc2ccccc2c2ccccc21",              "phenanthridine"),
    # regression: acridin-9(10H)-one unchanged (Phase 413)
    ("O=c1c2ccccc2[nH]c2ccccc12",        "acridin-9(10H)-one"),
    # regression: quinolin-4(1H)-one unchanged (Phase 416)
    ("O=c1cc[nH]c2ccccc12",              "quinolin-4(1H)-one"),
    # regression: benzene unchanged
    ("c1ccccc1",                           "benzene"),
])
def test_phase417_indazolone_phenanthridinone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
