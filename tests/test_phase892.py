"""Phase 892: large-ring heterocycles (11+ membered) -- crown ethers and
other macrocyclic polyethers/polyamines.

Found via the same fresh probe sweep as Phase 889-891.
C1COCCOCCOCCOCCO1 (a 15-membered ring with 5 oxygens, the pentaoxa
analogue of the crown-ether family) was named "cyclopentadecane" -- all 5
ring oxygens silently dropped, reporting the compound as a plain
hydrocarbon. The same bug affects 18-crown-6 itself
(C1COCCOCCOCCOCCOCCO1), one of the most famous macrocyclic ligands in
chemistry.

Root cause: _match_multi_het_ring (the general "a-nomenclature" heteroatom-
replacement namer for polyheteroatom saturated rings) only had suffix
entries for ring sizes 7-10 (epane/ocane/onane/ecane, the classical
Hantzsch-Widman-derived stems). Any larger ring returned None, which
propagated all the way up through name_heterocycle also returning None,
dropping the molecule into the generic carbocycle path that has no
representation for ring heteroatoms at all.

Fixed by extending the same function to ring sizes 11+, using the real
IUPAC convention for that range: "cyclo{alkane-stem}ane" (e.g.
cyclopentadecane for a 15-ring) instead of a dedicated HW stem, with the
same heteroatom-replacement locant/prefix logic already used for 7-10
rings. Caught (and fixed) an elision bug while implementing: the existing
7-10 code always dropped the trailing "a" of "oxa"/"aza" before the ring
stem (correct for epane/ocane/onane/ecane, which all start with a vowel)
but the new "cyclo..." stem starts with a consonant and must NOT be
elided ("oxacyclopentadecane", not "oxcyclopentadecane"). Verified via
OPSIN, including 18-crown-6 itself.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C1COCCOCCOCCOCCO1",     "1,4,7,10,13-pentaoxacyclopentadecane"),
    ("C1COCCOCCOCCOCCOCCO1",  "1,4,7,10,13,16-hexaoxacyclooctadecane"),  # 18-crown-6
    ("C1CNCCOCCOCC1",         "1,4-dioxa-7-azacycloundecane"),
    # regression: plain large carbocycle (no heteroatom) unchanged
    ("C1CCCCCCCCCC1", "cycloundecane"),
    # regression: existing 7-10 membered multi-heteroatom rings unchanged
    ("C1COCCOCO1",    "1,3,6-trioxocane"),
    ("C1CCOCCOCC1",   "1,4-dioxonane"),
])
def test_phase892_large_ring_heterocycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase892_18_crown_6_not_misnamed_as_hydrocarbon():
    result = smiles_to_iupac("C1COCCOCCOCCOCCOCCO1")
    assert "oxa" in result
    assert result != "cyclooctadecane"
