"""Phase 297: carbonohydrazide naming (IUPAC 2013 P-66.3.5.1.2).

H2N-NH-C(=O)-NH-NH2  →  carbonohydrazide
Retained trivial name; similar to urea/semicarbazide handling.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("NNC(=O)NN",          "carbonohydrazide"),
    # regressions: hydrazide unchanged
    ("CC(=O)NN",           "ethanohydrazide"),
    ("NNC(=O)N",           "semicarbazide"),
    ("NC(=O)N",            "urea"),
])
def test_phase297_carbonohydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
