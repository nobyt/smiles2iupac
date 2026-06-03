"""Phase 358: Gem-substituted cycloalkanol/cycloalkanone locant fix (IUPAC 2013).

When a substituent and the principal group (OH, =O) share the same ring locant,
the substituent locant must be shown explicitly.
Previously "1-methylcyclohexan-1-ol" was named "methylcyclohexan-1-ol".
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,1-gem: methyl and -ol at the same ring atom
    ("OC1(C)CCCCC1",             "1-methylcyclohexan-1-ol"),
    ("CC1(O)CCCCC1",             "1-methylcyclohexan-1-ol"),
    ("OC1(C)CCCC1",              "1-methylcyclopentan-1-ol"),
    ("OC1(C)CCC1",               "1-methylcyclobutan-1-ol"),
    # gem-dimethyl still works
    ("CC1(C)CCCCC1",             "1,1-dimethylcyclohexane"),
    # Regressions: mono-sub no suffix still omits locant
    ("CC1CCCCC1",                "methylcyclohexane"),
    ("OC1CCCCC1",                "cyclohexanol"),
    # Non-gem: locant stays
    ("OC1CCC(C)CC1",             "4-methylcyclohexan-1-ol"),
    ("OC1CCCCC1C",               "2-methylcyclohexan-1-ol"),
])
def test_phase358_gem_cycloalkanol_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
