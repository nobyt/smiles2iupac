"""Phase 879: fix Phase 877 regression -- biphenyl-substituent detection was
treating the PARENT ring as part of the substituent's own ring system.

Oc1ccc(-c2ccccc2)cc1 (4-phenylphenol) has phenol's own ring as the parent
(principal group -OH -> phenol) and a single plain phenyl ring as the
substituent at position 4. _find_benzene_biphenyl(graph) searches the WHOLE
molecule for any bonded pair of unsubstituted 6-membered aromatic rings --
which matches (phenol-ring, phenyl-ring) globally, even though phenol-ring
is the parent, not part of the substituent. The Phase 877 code then treated
the parent ring as though it belonged to the substituent, emitting the
nonsensical "[1,1'-biphenyl]-1-yl" (attachment at the ring-fusion carbon,
which already has 3 bonds -- chemically impossible for an aromatic carbon).

Fix: before taking the biphenyl-substituent path, require that the "other"
ring (the one NOT containing the substituent root) is disjoint from
`excluded` (the parent/already-claimed atoms). If the other ring overlaps
excluded, it's the parent itself, not part of this substituent -- fall
through to plain phenyl/substituted-phenyl naming.

NOTE: at the time this fix landed, that fallback produced the non-PIN
"4-phenylphenol"/"4-phenylaniline" form. Phase 880 (same session) then added
the biphenyl ol/amine PIN suffix, so these two cases now correctly resolve
to "[1,1'-biphenyl]-4-ol"/"[1,1'-biphenyl]-4-amine" via a *different* code
path (_name_biphenyl_carbo_suffix in ring_handler.py, dispatched before the
substituent-naming code this test file is about). The assertions below were
updated to match; the regression this file guards against -- the parent
ring being absorbed into a substituent name -- is still covered by the
"no chemically impossible attachment" check.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # the regression: OH/NH2-ring is parent, phenyl is a plain substituent.
    # (now resolves to the Phase 880 PIN form, not the "4-phenylphenol" this
    # fix originally restored -- see NOTE above)
    ("Oc1ccc(-c2ccccc2)cc1", "[1,1'-biphenyl]-4-ol"),
    ("Nc1ccc(-c2ccccc2)cc1", "[1,1'-biphenyl]-4-amine"),
    # regression: genuine biphenyl-as-substituent (Phase 877) still works
    ("CC(=O)c1ccc(-c2ccccc2)cc1", "1-([1,1'-biphenyl]-4-yl)ethan-1-one"),
    # regression: biphenyl as parent (Phase 874-876) unchanged
    ("c1ccccc1-c1ccccc1",         "1,1'-biphenyl"),
    ("Cc1ccc(-c2ccc(C)cc2)cc1",   "4,4'-dimethyl-1,1'-biphenyl"),
    ("OC(=O)c1ccc(-c2ccccc2)cc1", "[1,1'-biphenyl]-4-carboxylic acid"),
])
def test_phase879_biphenyl_parent_not_absorbed_into_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase879_no_chemically_impossible_biphenylyl_1yl():
    # attachment at the ring-fusion carbon of a biphenyl is impossible
    # (that carbon already has 3 substituents: 2 ring bonds + 1 inter-ring bond)
    result = smiles_to_iupac("Oc1ccc(-c2ccccc2)cc1")
    assert "biphenyl]-1-yl" not in result
    assert "(" not in result  # no substituent-style parenthetical at all
