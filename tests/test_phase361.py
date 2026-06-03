"""Phase 361: Correct chain selection for N-alkyl secondary amines.

When N is bonded to a secondary/tertiary carbon (e.g. propan-2-yl backbone),
the longest carbon chain *through* that carbon is used as the parent alkanamine,
and N is cited with its numeric locant (e.g. propan-2-amine).
Previously the chain was only collected in one direction, giving a shorter parent.

Also: N-substituents containing a locant digit are wrapped in parentheses per
IUPAC 2013 P-16.3.4 (e.g. N-(propan-2-yl)butan-2-amine, not N-propan-2-yl...).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N on interior carbon of 3-carbon chain
    ("CNC(C)C",        "N-methylpropan-2-amine"),
    ("CC(NC)C",        "N-methylpropan-2-amine"),
    # N on interior carbon of 4-carbon chain
    ("CC(NC)CC",       "N-methylbutan-2-amine"),
    ("CCC(NC)C",       "N-methylbutan-2-amine"),
    # Symmetric N-alkyl: both C neighbors give same interior chain
    ("CC(NC(C)C)C",    "N-(propan-2-yl)propan-2-amine"),
    # Different chain lengths: longest wins as parent
    ("CC(NC(C)C)CC",   "N-(propan-2-yl)butan-2-amine"),
    # Longest chain wins as parent (propan-2-amine > ethanamine)
    ("CCNC(C)C",       "N-ethylpropan-2-amine"),
    ("CNC(CC)C",       "N-methylbutan-2-amine"),
    # Regressions: terminal N unchanged
    ("CNC",            "N-methylmethanamine"),
    ("CNCC",           "N-methylethanamine"),
    ("CCNCC",          "N-ethylethanamine"),
    ("CCN(C)CC",       "N-ethyl-N-methylethanamine"),
    ("C/C=C/CNC",      "(2E)-N-methylbut-2-en-1-amine"),
    ("CNc1ccccc1",     "N-methylaniline"),
])
def test_phase361_secondary_amine_interior_chain(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
