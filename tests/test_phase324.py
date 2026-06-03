"""Phase 324: E/Z in imidic acid and thioamide names; N-sub imidic acid (IUPAC 2013 P-93.5).

Imidic acids and thioamides of alpha,beta-unsaturated chains carry E/Z prefix.
Imidic acid with N-substituent: N-{sub}{stem}imidic acid.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # E/Z imidic acid
    ("C/C=C/C(=N)O",            "(2E)-but-2-enimidic acid"),
    # N-substituted imidic acid
    ("CC(=NC)O",                "N-methylethanimidic acid"),
    ("CC(=NCC)O",               "N-ethylethanimidic acid"),
    # E/Z thioamide
    ("C/C=C/C(=S)N",            "(2E)-but-2-enethioamide"),
    (r"C/C=C\C(=S)N",           "(2Z)-but-2-enethioamide"),
    ("C/C=C/C(=S)NC",           "(2E)-N-methylbut-2-enethioamide"),
    # regressions: saturated imidic acid unchanged
    ("CC(=N)O",                 "ethanimidic acid"),
    ("CCC(=N)O",                "propanimidic acid"),
    # regressions: saturated thioamide unchanged
    ("CC(=S)N",                 "ethanethioamide"),
    ("CC(=S)NC",                "N-methylethanethioamide"),
])
def test_phase324_ez_imidic_thioamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
