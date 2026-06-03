"""Phase 231: silanol naming (IUPAC 2013 P-68.4.2).

R_n-Si(OH)_{4-n} → {alkyl(s)}silanol (for mono-OH case).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # trimethylsilanol: 3 methyl + 1 OH
    ("C[Si](C)(C)O",   "trimethylsilanol"),
    # dimethylsilanol: 2 methyl + 1 OH
    ("C[Si](O)(O)C",   "dimethylsilanol"),
    # regression: silane still works
    ("C[Si](C)(C)C",   "tetramethylsilane"),
    ("C[Si](C)(C)CC",  "ethyltrimethylsilane"),
])
def test_phase231_silanol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
