"""Phase 880: biphenyl ol/amine PIN suffix ([1,1'-biphenyl]-N-ol / -N-amine).

Extends Phase 876's carbo-family biphenyl suffix (COOH/CHO/CN) to -OH and
-NH2 attached directly to a biphenyl ring, per IUPAC 2013 P-28.2.1 (biphenyl
as parent hydride). Previously these fell to the non-PIN "4-phenylphenol" /
"4-phenylaniline" functional-class-substituent form.

Implementation: generalized _name_biphenyl_carbo_suffix's suffix table and
fixed its anchor-detection loop, which assumed the principal-group atoms
were always exocyclic (true for COOH/CHO/CN, where pgrp.atom_indices holds
only the acid/aldehyde/nitrile atoms one hop from the ring). For alcohol/
amine, pgrp.atom_indices ALSO includes the ring carbon itself (it's directly
bonded to O/N), so the old loop searched that ring carbon's neighbors for a
"ring member" and could pick an adjacent ring atom instead of itself,
silently shifting the computed locant off by one (para wrongly came out as
"3" instead of "4"). Fixed by checking `if a in ring_all: anchor = a` before
falling back to the neighbor search.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # para
    ("Oc1ccc(-c2ccccc2)cc1", "[1,1'-biphenyl]-4-ol"),
    ("Nc1ccc(-c2ccccc2)cc1", "[1,1'-biphenyl]-4-amine"),
    # meta
    ("Oc1cccc(-c2ccccc2)c1", "[1,1'-biphenyl]-3-ol"),
    ("Nc1cccc(-c2ccccc2)c1", "[1,1'-biphenyl]-3-amine"),
    # ortho
    ("Oc1ccccc1-c1ccccc1",   "[1,1'-biphenyl]-2-ol"),
    ("Nc1ccccc1-c1ccccc1",   "[1,1'-biphenyl]-2-amine"),
    # regression: plain phenol/aniline (single ring) unchanged
    ("Oc1ccccc1",     "phenol"),
    ("Nc1ccccc1",     "aniline"),
    ("Oc1ccc(C)cc1",  "4-methylphenol"),
    ("Nc1ccc(C)cc1",  "4-methylaniline"),
    # regression: biphenyl carbo-suffix (Phase 876) unaffected -- same anchor
    # code path, must still compute the correct para/meta/ortho locant
    ("OC(=O)c1ccc(-c2ccccc2)cc1", "[1,1'-biphenyl]-4-carboxylic acid"),
    ("O=Cc1ccc(-c2ccccc2)cc1",    "[1,1'-biphenyl]-4-carbaldehyde"),
    ("N#Cc1ccc(-c2ccccc2)cc1",    "[1,1'-biphenyl]-4-carbonitrile"),
    # regression: biphenyl as substituent (Phase 877) unaffected
    ("CC(=O)c1ccc(-c2ccccc2)cc1", "1-([1,1'-biphenyl]-4-yl)ethan-1-one"),
    # regression: biphenyl as plain parent (Phase 874/875) unaffected
    ("c1ccccc1-c1ccccc1",         "1,1'-biphenyl"),
    ("Cc1ccc(-c2ccc(C)cc2)cc1",   "4,4'-dimethyl-1,1'-biphenyl"),
])
def test_phase880_biphenyl_ol_amine_suffix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase880_ol_amine_scoped_to_unsubstituted_ring():
    # extra ring substituents are out of scope for this PIN path; must not
    # crash and must not silently claim a wrong biphenyl-ol/amine name.
    result = smiles_to_iupac("Oc1ccc(-c2ccccc2)cc1C")
    assert result != "[1,1'-biphenyl]-4-ol"
