"""Phase 364: Heteroatom spiro compounds (azaspiro, oxaspiro, thiaspiro).

Previously _try_spiro rejected any ring containing N, O, or S atoms, falling
back to an incorrect substituted-ring name. Now it mirrors the _try_bicyclo
logic: allows C/N/O/S ring atoms, identifies the spiro carbon (must be C),
picks the numbering that gives heteroatoms the lowest possible locants, and
prepends the aza/oxa/thia prefix.

Also fixed _het_a_prefix to join multiple heteroatom-type parts with '-'
instead of '' (e.g. '1-oxa-2-aza' not '1-oxa2-aza').
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N in smaller ring, adjacent to spiro → locant 1
    ("C1CCC2(CC1)CCCN2",    "1-azaspiro[4.5]decane"),
    # N in larger ring, equidistant from spiro → locant 8
    ("N1CCC2(CC1)CCCC2",    "8-azaspiro[4.5]decane"),
    # O in larger ring
    ("O1CCC2(CC1)CCCC2",    "8-oxaspiro[4.5]decane"),
    # O in 7-membered ring, adjacent to spiro → locant 6
    ("C1CCOC2(CC1)CCCC2",   "6-oxaspiro[4.6]undecane"),
    # Both O and N in the 5-membered ring (O adjacent to spiro)
    ("C1CCC2(CC1)CCNO2",    "1-oxa-2-azaspiro[4.5]decane"),
    # Equal-sized rings: N in ring → numbered as small to get lower locant
    ("C1CC2(CCCC2)CN1",     "2-azaspiro[4.4]nonane"),
    # N in 5-ring, other ring is 4-membered → 6-azaspiro[3.4]octane
    ("C1CC2(CCC2)CN1",      "6-azaspiro[3.4]octane"),
    # Substituted hetero-spiro
    ("CC1CCC2(CC1)CCCN2",   "8-methyl-1-azaspiro[4.5]decane"),
    # Regressions: pure-carbon spiro unchanged
    ("C1CCC2(CC1)CCCC2",    "spiro[4.5]decane"),
    ("C1CCC2(CC1)CCC2",     "spiro[3.5]nonane"),
    ("O=C1CCC2(CC1)CCCC2",  "spiro[4.5]decan-8-one"),
])
def test_phase364_heteroatom_spiro(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
