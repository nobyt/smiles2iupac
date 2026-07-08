"""Phase 852: N-methyl thiolactams → -thione suffix (no indicated-H when N is substituted).

When an N-heterocycle thione has N substituted, the indicated-H in the suffix is
omitted (P-14.7.1). These compounds arrive via the sulfanyl path (RDKit converts
exo C=S to C-SH tautomer). Added a generic (nH) stripping step at the end of the
hydroxy/sulfanyl block in _apply_hetero_suffixes.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Quinoline: N-methyl thiones (no (1H))
    ("Cn1ccc(=S)c2ccccc21",          "1-methylquinolin-4-thione"),
    ("Cn1c(=S)ccc2ccccc21",          "1-methylquinolin-2-thione"),
    # Isoquinoline: N-methyl thiones (no (2H))
    ("Cn1ccc2ccccc2c1=S",            "2-methylisoquinolin-1-thione"),
    ("Cn1cc2ccccc2cc1=S",            "2-methylisoquinolin-3-thione"),
    # Phthalazine: N-2 methyl (no (2H))
    ("Cn1ncc2ccccc2c1=S",            "2-methylphthalazin-1-thione"),
    # Quinazoline: N-3 methyl (no (3H))
    ("Cn1cnc2ccccc2c1=S",            "3-methylquinazolin-4-thione"),
    # Cinnoline: N-2 methyl (no (2H))
    ("Cn1nc2ccccc2cc1=S",            "2-methylcinnolin-3-thione"),
    # Quinoxaline: N-1 methyl (no (1H))
    ("Cn1c(=S)cnc2ccccc21",          "1-methylquinoxalin-2-thione"),
    # Benzimidazole: N-3 methyl drops (3H); 1H- prefix retained
    ("Cn1c(=S)[nH]c2ccccc21",        "3-methyl-1H-benzimidazol-2-thione"),
    # Phenanthridine: N-5 methyl (no (5H))
    ("Cn1c(=S)c2ccccc2c2ccccc21",    "5-methylphenanthridin-6-thione"),
    # NH forms: indicated-H is retained (N has no substituent)
    ("Sc1ccc2ccccc2n1",              "quinolin-2(1H)-thione"),
    ("Sc1ccnc2ccccc12",              "quinolin-4(1H)-thione"),
    ("Sc1ccccn1",                    "pyridin-2(1H)-thione"),
    ("Sc1ccncc1",                    "pyridin-4(1H)-thione"),
])
def test_phase852(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
