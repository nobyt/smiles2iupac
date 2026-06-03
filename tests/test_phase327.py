"""Phase 327: E/Z in dicarboxylate and imidate ester names (IUPAC 2013 P-93.5).

Dicarboxylate dianions and imidate esters of alpha,beta-unsaturated chains
now carry the E/Z prefix. Unsaturated dicarboxylates use systematic names
('but-2-enedioate') rather than retained names ('succinate').
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # dicarboxylate E/Z
    ("[O-]C(=O)/C=C/C(=O)[O-]",     "(2E)-but-2-enedioate"),
    (r"[O-]C(=O)/C=C\C(=O)[O-]",    "(2Z)-but-2-enedioate"),
    # imidate ester E/Z
    ("C/C=C/C(=N)OC",               "methyl (2E)-but-2-enimidate"),
    (r"C/C=C\C(=N)OC",              "methyl (2Z)-but-2-enimidate"),
    # regressions: saturated dicarboxylates use retained names
    ("[O-]C(=O)CCC(=O)[O-]",        "succinate"),
    ("[O-]C(=O)CC(=O)[O-]",         "malonate"),
    ("[O-]C(=O)C(=O)[O-]",          "oxalate"),
    # regressions: saturated imidate esters unchanged
    ("CC(=N)OC",                    "methyl ethanimidate"),
    ("CCC(=N)OCC",                  "ethyl propanimidate"),
])
def test_phase327_ez_dicarboxylate_imidate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
