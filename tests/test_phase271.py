"""Phase 271: chain direction tie-break — multiple bond locant before substituent locant
(IUPAC 2013 P-44.1 numbering priority: C→D→E→F order).

When the principal characteristic group locant is tied in both chain directions,
the double bond locant must be minimised before falling back to substituent locants.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # prop-1-en-2-ol: OH at C2 tied, db at C1-C2 (position 1) wins over C2-C3 (position 2)
    ("CC(O)=C",      "prop-1-en-2-ol"),
    ("C=C(C)O",      "prop-1-en-2-ol"),   # same molecule, different SMILES atom order
    ("OC(=C)C",      "prop-1-en-2-ol"),
    # but-1-en-2-ol (= methyl vinyl carbinol)
    ("C=C(O)CC",     "but-1-en-2-ol"),
    # ethenol: no tie (trivial, regression)
    ("C=CO",         "ethenol"),
    # regression: prop-2-en-1-ol unchanged
    ("C=CCO",        "prop-2-en-1-ol"),
    ("OCC=C",        "prop-2-en-1-ol"),
    # regression: alkynes unchanged
    ("C#CCO",        "prop-2-yn-1-ol"),
    # regression: ketones unchanged
    ("CC(=O)CO",     "1-hydroxypropan-2-one"),
    ("CC(=O)CBr",    "1-bromopropan-2-one"),
])
def test_phase271_chain_direction_tiebreak(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
