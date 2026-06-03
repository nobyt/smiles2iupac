"""Phase 330: E/Z in telluramide; dicarboxylate locant optimization (IUPAC 2013).

Telluramide chains now carry E/Z prefix.
Unsaturated dicarboxylate dianions get the lowest possible locant for the double bond.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # telluramide E/Z
    ("C/C=C/C(=[Te])N",          "(2E)-but-2-eneteluramide"),
    (r"C/C=C\C(=[Te])N",         "(2Z)-but-2-eneteluramide"),
    # dicarboxylate locant optimization (no stereo marks)
    ("[O-]C(=O)CC=CC(=O)[O-]",   "pent-2-enedioate"),
    ("[O-]C(=O)C=CCC(=O)[O-]",   "pent-2-enedioate"),
    # regressions: saturated telluramide unchanged
    ("CC(=[Te])N",               "ethaneteluramide"),
    ("CCC(=[Te])N",              "propaneteluramide"),
    # regressions: saturated dicarboxylates use retained names
    ("[O-]C(=O)CCC(=O)[O-]",     "succinate"),
    ("[O-]C(=O)CC(=O)[O-]",      "malonate"),
    # regressions: stereo-marked dicarboxylate (phase 327)
    ("[O-]C(=O)/C=C/C(=O)[O-]",  "(2E)-but-2-enedioate"),
    (r"[O-]C(=O)/C=C\C(=O)[O-]", "(2Z)-but-2-enedioate"),
])
def test_phase330_telluramide_ez_dicarboxylate_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
