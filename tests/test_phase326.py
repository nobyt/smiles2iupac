"""Phase 326: E/Z in carboxylate, diacid halide names (IUPAC 2013 P-93.5).

Carboxylate anions and diacid halides of alpha,beta-unsaturated chains
now carry the E/Z prefix. Also adds unsaturated (enedioyl) suffix to
diacid halides.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # carboxylate E/Z
    ("C/C=C/C(=O)[O-]",             "(2E)-but-2-enoate"),
    (r"C/C=C\C(=O)[O-]",            "(2Z)-but-2-enoate"),
    # diacid halide E/Z + enedioyl suffix
    ("ClC(=O)/C=C/C(=O)Cl",         "(2E)-but-2-enedioyl dichloride"),
    (r"ClC(=O)/C=C\C(=O)Cl",        "(2Z)-but-2-enedioyl dichloride"),
    # regressions: saturated carboxylate unchanged
    ("CCC(=O)[O-]",                  "propanoate"),
    ("CC(=O)[O-]",                   "acetate"),
    # regressions: saturated diacid halide unchanged
    ("ClC(=O)CC(=O)Cl",             "propanedioyl dichloride"),
    ("ClC(=O)C(=O)Cl",              "ethanedioyl dichloride"),
    ("BrC(=O)CC(=O)Br",             "propanedioyl dibromide"),
])
def test_phase326_ez_carboxylate_diacid_halide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
