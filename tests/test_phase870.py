"""Phase 870: N-hydroxy thio/seleno/telluramides.

The N-hydroxy substituent on a thioamide/selenoamide/telluramide was dropped:
CC(=S)NO was named "ethanethioamide" instead of "N-hydroxyethanethioamide".
The oxygen analogue CC(=O)NO already gives "N-hydroxyacetamide", so this brings
the chalcogen amides in line by collecting the N-OH as an "N-hydroxy"
substituent alongside any N-alkyl groups (alphabetically ordered, so
N-hydroxy-N-methyl...).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N-hydroxy chalcogen amides
    ("CC(=S)NO",       "N-hydroxyethanethioamide"),
    ("CC(=[Se])NO",    "N-hydroxyethaneselenoamide"),
    ("CC(=[Te])NO",    "N-hydroxyethanetelluramide"),
    ("CCC(=S)NO",      "N-hydroxypropanethioamide"),
    # N-hydroxy + N-alkyl (alphabetical: hydroxy before methyl)
    ("CC(=S)N(C)O",    "N-hydroxy-N-methylethanethioamide"),
    ("CC(=[Se])N(C)O", "N-hydroxy-N-methylethaneselenoamide"),
    # regression: oxygen amide unchanged
    ("CC(=O)NO",       "N-hydroxyacetamide"),
    ("CC(=O)N(C)O",    "N-hydroxy-N-methylacetamide"),
    # regression: plain and N-alkyl chalcogen amides unchanged
    ("CC(=S)N",        "ethanethioamide"),
    ("CC(=S)NC",       "N-methylethanethioamide"),
    ("CC(=S)NCC",      "N-ethylethanethioamide"),
    ("CC(=[Se])N",     "ethaneselenoamide"),
    ("CC(=[Te])N",     "ethanetelluramide"),
    ("NC(=S)N",        "thiourea"),
])
def test_phase870_n_hydroxy_chalcogen_amides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
