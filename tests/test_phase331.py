"""Phase 331: Unsaturated dithioamide/diselenoamide naming; diester locant optimization.

Dithioamide and diselenoamide chains now get proper unsaturated names with lowest locants.
Diester unsaturated chains get the lowest possible locant for the double bond.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # dithioamide unsaturated
    ("NC(=S)CC=CC(=S)N",           "pent-2-enedithioamide"),
    ("NC(=S)C=CCC(=S)N",           "pent-2-enedithioamide"),
    # diselenoamide unsaturated
    ("NC(=[Se])CC=CC(=[Se])N",     "pent-2-enediselenoamide"),
    ("NC(=[Se])C=CCC(=[Se])N",     "pent-2-enediselenoamide"),
    # diester locant optimization
    ("COC(=O)CC=CC(=O)OC",         "dimethyl pent-2-enedioate"),
    ("COC(=O)C=CCC(=O)OC",         "dimethyl pent-2-enedioate"),
    # regressions: saturated dithioamide/diselenoamide unchanged
    ("NC(=S)C(=S)N",               "ethanedithioamide"),
    ("NC(=S)CC(=S)N",              "propanedithioamide"),
    ("NC(=[Se])C(=[Se])N",         "ethanediselenoamide"),
    # regressions: diester with stereo / retained names
    ("COC(=O)/C=C/C(=O)OC",        "(2E)-dimethyl but-2-enedioate"),
    ("COC(=O)CCC(=O)OC",           "dimethyl succinate"),
    ("CCOC(=O)CC(=O)OCC",          "diethyl malonate"),
])
def test_phase331_unsaturated_dithio_diselenoamide_diester_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
