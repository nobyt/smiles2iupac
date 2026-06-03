"""Phase 317: E/Z stereodescriptor in acid halide and thioester names (IUPAC 2013 P-93.5).

Acid halides and thioesters of α,β-unsaturated acids now carry the E/Z prefix.
Also fixes thioester double-bond suffix (was 'anethioate', now 'enethioate').
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # acid halides
    ("C/C=C/C(=O)Cl",    "(2E)-but-2-enoyl chloride"),
    ("C/C=C/C(=O)Br",    "(2E)-but-2-enoyl bromide"),
    ("C/C=C\\C(=O)Cl",   "(2Z)-but-2-enoyl chloride"),
    # thioesters
    ("C/C=C/C(=O)SC",    "(2E)-S-methyl but-2-enethioate"),
    ("C/C=C/C(=O)SCC",   "(2E)-S-ethyl but-2-enethioate"),
    # regressions: saturated acid halides unchanged
    ("CC(=O)Cl",          "acetyl chloride"),
    ("CCC(=O)Cl",         "propanoyl chloride"),
    # regressions: saturated thioesters unchanged
    ("CC(=O)SC",          "S-methyl ethanethioate"),
    ("CCC(=O)SC",         "S-methyl propanethioate"),
])
def test_phase317_ez_acid_halide_thioester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
