"""Phase 362: Diketones on the same aromatic ring are not merged into dione.

When two ketone groups are both directly attached to the same benzene ring
(e.g. 1,4-diacetylbenzene), the functional group aggregator previously
incorrectly merged them into a "dione" group, producing garbage output like
"1-phenylethane-2,4-dione". The fix: do not merge two aryl ketones into a
dione when both C=O groups are adjacent to the same aromatic ring.

The correct IUPAC 2013 preferred name uses the ketone chain as parent and the
substituted phenyl ring as substituent (same logic as acetophenone →
1-phenylethan-1-one): "1-(4-acetylphenyl)ethan-1-one" for the para isomer.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # para diacetylbenzene
    ("CC(=O)c1ccc(cc1)C(=O)C",   "1-(4-acetylphenyl)ethan-1-one"),
    # ortho diacetylbenzene
    ("CC(=O)c1ccccc1C(=O)C",     "1-(2-acetylphenyl)ethan-1-one"),
    # meta diacetylbenzene
    ("CC(=O)c1cccc(C(=O)C)c1",   "1-(3-acetylphenyl)ethan-1-one"),
    # Regressions: chain diones unchanged
    ("CC(=O)CC(=O)C",            "pentane-2,4-dione"),
    ("CCC(=O)CC(=O)CC",          "heptane-3,5-dione"),
    ("CC(=O)CCCC(=O)C",          "heptane-2,6-dione"),
    # Single aryl ketone unchanged
    ("CC(=O)c1ccccc1",           "acetophenone"),
    ("CC(=O)c1ccc(C)cc1",        "1-(4-methylphenyl)ethan-1-one"),
    # Cyclic ketones unchanged
    ("O=C1CCCCC1",               "cyclohexanone"),
    ("O=C1CCCC1",                "cyclopentanone"),
])
def test_phase362_diaryl_dione_no_merge(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
