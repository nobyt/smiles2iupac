"""Phase 329: E/Z in amidine/selenoamide; diacid halide locant optimization (IUPAC 2013).

Amidimidamide (amidine) and selenoamide chains now carry E/Z prefix.
Unsaturated diacid halides get the lowest possible locant for the double bond.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # amidine E/Z
    ("C/C=C/C(=N)N",         "(2E)-but-2-enimidamide"),
    (r"C/C=C\C(=N)N",        "(2Z)-but-2-enimidamide"),
    ("C/C=C/CC(=N)N",        "(3E)-pent-3-enimidamide"),
    # selenoamide E/Z
    ("C/C=C/C(=[Se])N",      "(2E)-but-2-eneselenoamide"),
    (r"C/C=C\C(=[Se])N",     "(2Z)-but-2-eneselenoamide"),
    # diacid halide locant optimization: double bond gets lowest locant
    ("ClC(=O)CC=CC(=O)Cl",   "pent-2-enedioyl dichloride"),
    ("ClC(=O)C=CCC(=O)Cl",   "pent-2-enedioyl dichloride"),
    # regressions: saturated amidine/selenoamide unchanged
    ("CC(=N)N",              "ethanimidamide"),
    ("CCC(=N)N",             "propanimidamide"),
    ("CC(=[Se])N",           "ethaneselenoamide"),
    # regressions: saturated diacid halide unchanged
    ("ClC(=O)CC(=O)Cl",      "propanedioyl dichloride"),
    ("ClC(=O)CCC(=O)Cl",     "butanedioyl dichloride"),
])
def test_phase329_ez_amidine_selenoamide_diacid_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
