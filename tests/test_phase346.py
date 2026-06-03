"""Phase 346: O-substituted hydroxylamine naming (IUPAC 2013).

Compounds of the form NH2-O-R are named as O-alkylhydroxylamine.
E/Z-bearing alkenyl O-substituents get proper outer parentheses.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # O-alkyl hydroxylamines
    ("NOC",                     "O-methylhydroxylamine"),
    ("NOCC",                    "O-ethylhydroxylamine"),
    ("NOCCC",                   "O-propylhydroxylamine"),
    ("NOc1ccccc1",              "O-phenylhydroxylamine"),
    # O-alkenyl hydroxylamine with E/Z
    ("C/C=C/CON",               "O-[(2E)-but-2-en-1-yl]hydroxylamine"),
    (r"C/C=C\CON",              "O-[(2Z)-but-2-en-1-yl]hydroxylamine"),
    # regressions: N-substituted hydroxylamines unchanged
    ("NO",                      "hydroxylamine"),
    ("CNO",                     "N-methylhydroxylamine"),
    ("CCNO",                    "N-ethylhydroxylamine"),
    ("CN(O)C",                  "N,N-dimethylhydroxylamine"),
])
def test_phase346_o_substituted_hydroxylamine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
