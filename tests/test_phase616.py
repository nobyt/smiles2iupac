"""Phase 616: simple NH monocyclics — methyl-substituted pyrazole, triazoles, tetrazole."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-pyrazole
    ("c1cn[nH]c1",     "1H-pyrazole"),
    ("Cc1cn[nH]c1",    "4-methyl-1H-pyrazole"),
    ("Cc1cc[nH]n1",    "3-methyl-1H-pyrazole"),
    ("Cc1ccn[nH]1",    "5-methyl-1H-pyrazole"),
    # 1H-1,2,3-triazole
    ("c1c[nH]nn1",     "1H-1,2,3-triazole"),
    ("Cc1c[nH]nn1",    "4-methyl-1H-1,2,3-triazole"),
    ("Cc1cnn[nH]1",    "5-methyl-1H-1,2,3-triazole"),
    # 1H-1,2,4-triazole
    ("c1nc[nH]n1",     "1H-1,2,4-triazole"),
    ("Cc1nc[nH]n1",    "3-methyl-1H-1,2,4-triazole"),
    ("Cc1ncn[nH]1",    "5-methyl-1H-1,2,4-triazole"),
    # 1H-tetrazole
    ("c1nn[nH]n1",     "1H-tetrazole"),
    ("Cc1nn[nH]n1",    "5-methyl-1H-tetrazole"),
    ("c1nnn[nH]1",     "1H-tetrazole"),
    ("Cc1nnn[nH]1",    "5-methyl-1H-tetrazole"),
    # 2H-1,2,3-triazole
    ("c1cn[nH]n1",     "2H-1,2,3-triazole"),
    ("Cc1cn[nH]n1",    "4-methyl-2H-1,2,3-triazole"),
])
def test_phase616_nh_monocyclics(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
