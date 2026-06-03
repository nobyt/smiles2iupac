"""Phase 295: propene double bond locant in substituted derivatives (IUPAC 2013 P-31.1.2.2).

Unsubstituted propene → locant omitted: "propene"
Substituted 3-carbon ene → locant cited: "3-chloroprop-1-ene"
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # substituted: locant must be cited
    ("ClCC=C",   "3-chloroprop-1-ene"),
    ("BrCC=C",   "3-bromoprop-1-ene"),
    ("FCC=C",    "3-fluoroprop-1-ene"),
    # regressions: unsubstituted propene keeps omitted locant
    ("C=CC",     "propene"),
    ("CC=C",     "propene"),
    # branch on C2 of propene: locant still omitted (2-methylprop-1-ene vs 2-methylpropene)
    # IUPAC 2013 P-31.1.2.2: locant omitted only for parent chain "propene" with no substituents
    ("CC(=C)C",  "2-methylprop-1-ene"),
    # existing correct cases: suffix forces locant numbering
    ("C=CCO",    "prop-2-en-1-ol"),
    ("C=CCC",    "but-1-ene"),
])
def test_phase295_propene_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
