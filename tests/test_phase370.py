"""Phase 370: Disulfide preferred IUPAC 2013 names (dialkyl format).

Disulfides (R-S-S-R') keep the functional-class 'dialkyl disulfide' form here
(no simple substitutive equivalent is commonly used for the -S-S- linkage as
a substituent chain the way -S- is for sulfides).

Sulfides (R-S-R', single S) were reverted to functional-class naming by this
phase at the time, but that was a mistake: IUPAC 2013 treats sulfanyl exactly
like alkoxy (P-63.1.5 parallels P-63.1.4) — the substitutive '(alkylsulfanyl)
parent' form is the preferred IUPAC name, and 'R R' sulfide' is only an
acceptable alternative, exactly as '(alkoxy)parent' (e.g. methoxymethane) is
preferred over 'dialkyl ether'. Phase 855 restored substitutive sulfide
naming; only the plain-sulfide rows below were corrected, the disulfide rows
are untouched.

Rules (disulfides only, same as for sulfone/sulfoxide, Phase 369):
  - Same groups: di<group> disulfide
  - Different groups: <group1> <group2> disulfide (alphabetical)
  - Complex group names (digits or parentheses) → wrap in parentheses/brackets
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Symmetric sulfides
    ("CSC",                   "(methylsulfanyl)methane"),
    ("CCSCC",                 "(ethylsulfanyl)ethane"),
    # Asymmetric sulfides (longer chain is parent)
    ("CSCC",                  "(methylsulfanyl)ethane"),
    ("CSCCC",                 "1-(methylsulfanyl)propane"),
    ("c1ccccc1Sc1ccccc1",     "(phenylsulfanyl)benzene"),
    ("c1ccccc1SC",            "(methylsulfanyl)benzene"),
    # Branched group substituent
    ("CC(C)SC",               "2-(methylsulfanyl)propane"),
    # Vinyl sulfides
    ("C=CSC",                 "(methylsulfanyl)ethene"),
    # Complex (E/Z) substituent
    ("C/C=C/CSC",             "(2E)-1-(methylsulfanyl)but-2-ene"),
    # Symmetric disulfides
    ("CSSC",                  "dimethyl disulfide"),
    ("CCSSCC",                "diethyl disulfide"),
    # Asymmetric disulfide
    ("CSSCC",                 "ethyl methyl disulfide"),
    # Complex disulfide
    ("C/C=C/CSSC",            "[(2E)-but-2-en-1-yl] methyl disulfide"),
    # Regressions: polysulfide (trisulfide) unchanged
    ("CSSSC",                 "dimethyl trisulfide"),
    # Ring thioether unchanged
    ("C1CCSC1",               "thiolane"),
])
def test_phase370_sulfide_disulfide_preferred(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
