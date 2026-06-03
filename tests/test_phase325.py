"""Phase 325: E/Z stereodescriptor in hydrazide names (IUPAC 2013 P-93.5).

Hydrazides of alpha,beta-unsaturated acids now carry the E/Z prefix.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C/C=C/C(=O)NN",            "(2E)-but-2-enohydrazide"),
    (r"C/C=C\C(=O)NN",           "(2Z)-but-2-enohydrazide"),
    ("C/C=C/C(=O)NNC",           "(2E)-N'-methylbut-2-enohydrazide"),
    # regressions: saturated hydrazides unchanged
    ("CC(=O)NN",                  "ethanohydrazide"),
    ("CCC(=O)NN",                 "propanohydrazide"),
    ("CC(=O)NNC",                 "N'-methylethanohydrazide"),
])
def test_phase325_ez_hydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
