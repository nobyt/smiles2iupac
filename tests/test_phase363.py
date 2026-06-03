"""Phase 363: Chain substituents in secondary/tertiary amine naming.

When the parent chain of a secondary or tertiary amine has branch substituents
(e.g. a methyl group at position 2), those were previously missing from the
name. The fix: collect chain substituents and merge them with N-substituents,
sorted alphabetically.

E.g. CC(C)CNCC  → N-ethyl-2-methylpropan-1-amine  (methyl at C2, N-ethyl)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Chain methyl at position 2 (isobutyl-N backbone)
    ("CC(C)CNCC",          "N-ethyl-2-methylpropan-1-amine"),
    # Same-name substituents: numeric locant before N locant
    ("CC(C)CNC",           "2-methyl-N-methylpropan-1-amine"),
    # Same isobutyl backbone, longer N-chain: propyl wins as parent
    ("CC(C)CNCCC",         "2-methyl-N-propylpropan-1-amine"),
    # Longer chain (butyl=4C) wins over isobutyl (3C): 2-methylpropyl becomes N-sub
    ("CC(C)CNCCCC",        "N-(2-methylpropyl)butan-1-amine"),
    # Primary amine with branched chain (uses general path, not this one)
    ("CC(C)CN",            "2-methylpropan-1-amine"),
    # Regressions: linear chains unchanged
    ("CNC",                "N-methylmethanamine"),
    ("CNCC",               "N-methylethanamine"),
    ("CCNCC",              "N-ethylethanamine"),
    ("CCN(C)CC",           "N-ethyl-N-methylethanamine"),
    ("CNC(C)C",            "N-methylpropan-2-amine"),
    ("C/C=C/CNC",          "(2E)-N-methylbut-2-en-1-amine"),
])
def test_phase363_amine_chain_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
