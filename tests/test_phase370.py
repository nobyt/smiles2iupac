"""Phase 370: Sulfide/disulfide preferred IUPAC 2013 names (dialkyl format).

Previously _name_sulfide and _name_disulfide used substitutive forms like
'(methylsulfanyl)methane' and '(methyldisulfanyl)methane'. IUPAC 2013 P-63.6.1
/ P-63.7.1 require the functional-class form: group names + 'sulfide'/'disulfide'.

Rules (same as for sulfone/sulfoxide, Phase 369):
  - Same groups: di<group> sulfide/disulfide
  - Different groups: <group1> <group2> sulfide/disulfide (alphabetical)
  - Complex group names (digits or parentheses) → wrap in parentheses/brackets
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Symmetric sulfides
    ("CSC",                   "dimethyl sulfide"),
    ("CCSCC",                 "diethyl sulfide"),
    # Asymmetric sulfides (alphabetical order)
    ("CSCC",                  "ethyl methyl sulfide"),
    ("CSCCC",                 "methyl propyl sulfide"),
    ("c1ccccc1Sc1ccccc1",     "diphenyl sulfide"),
    ("c1ccccc1SC",            "methyl phenyl sulfide"),
    # Branched group → parens
    ("CC(C)SC",               "methyl (propan-2-yl) sulfide"),
    # Vinyl sulfides
    ("C=CSC",                 "ethenyl methyl sulfide"),
    # Complex (E/Z) groups → brackets
    ("C/C=C/CSC",             "[(2E)-but-2-en-1-yl] methyl sulfide"),
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
