"""Phase 316: E/Z stereodescriptor in ester names (IUPAC 2013 P-93.5).

Esters of α,β-unsaturated acids now carry the E/Z prefix.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("COC(=O)/C=C/C",     "(2E)-methyl but-2-enoate"),
    ("CCOC(=O)/C=C/C",    "(2E)-ethyl but-2-enoate"),
    ("COC(=O)/C=C\\C",    "(2Z)-methyl but-2-enoate"),
    ("COC(=O)/C=C/CC",    "(2E)-methyl pent-2-enoate"),
    # regressions: saturated esters unchanged
    ("COC(=O)CC",         "methyl propanoate"),
    ("CCOC(=O)C",         "ethyl acetate"),
    ("COC(=O)C",          "methyl acetate"),
])
def test_phase316_ez_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
