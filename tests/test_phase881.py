"""Phase 881: unsubstituted linear terphenyl PIN (o-/m-/p-terphenyl).

Extends the biphenyl ring-assembly work (Phase 874-880) to 3-ring linear
assemblies. Previously a plain para-terphenyl (three benzene rings, single
bonds central-A and central-B) was named "1,4-diphenylbenzene" (a valid but
non-PIN functional-substituent form); the IUPAC 2013 P-28.2.1 PIN uses the
terphenyl parent hydride: "1,1':4',1''-terphenyl" (para),
"1,1':3',1''-terphenyl" (meta), "1,1':2',1''-terphenyl" (ortho).

Implementation: _find_terphenyl (ring_handler.py) detects 3 non-fused
all-carbon aromatic 6-rings where one ring (central) bridges to each of the
other two (outer) via exactly one single bond apiece, and the two outer
rings are not bonded to each other (i.e. a genuine linear chain, not a
triangle/branch). _name_terphenyl_assembly then numbers the central ring
starting at its first bridge carbon (locant 1'), walks the ring, and reports
where the second bridge carbon falls (2'/3'/4' for ortho/meta/para) -- the
same "count ring positions apart, take the smaller of the two directions"
formula already used for biphenyl locants. Scoped to the fully unsubstituted
case only (matching the narrow-scoping precedent from the biphenyl phases);
any additional ring substituent falls through unchanged (no crash, no wrong
claim).
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # para (linear, 1,4 on central ring)
    ("c1ccc(-c2ccc(-c3ccccc3)cc2)cc1", "1,1':4',1''-terphenyl"),
    # meta (1,3 on central ring)
    ("c1ccc(-c2cccc(-c3ccccc3)c2)cc1", "1,1':3',1''-terphenyl"),
    # ortho (1,2 on central ring)
    ("c1ccc(-c2ccccc2-c2ccccc2)cc1",   "1,1':2',1''-terphenyl"),
    # regression: plain biphenyl (2 rings) unchanged
    ("c1ccccc1-c1ccccc1",              "1,1'-biphenyl"),
    ("Cc1ccc(-c2ccc(C)cc2)cc1",        "4,4'-dimethyl-1,1'-biphenyl"),
    # regression: naphthalene/anthracene (fused, not ring-assemblies) unchanged
    ("c1ccc2ccccc2c1",                 "naphthalene"),
    ("c1ccc2cc3ccccc3cc2c1",           "anthracene"),
    # regression: triphenylmethane (3 rings, but branched off sp3 C, not a
    # direct ring-to-ring chain) must not be misdetected as terphenyl
    ("C(c1ccccc1)(c1ccccc1)c1ccccc1",  "(diphenylmethyl)benzene"),
])
def test_phase881_terphenyl_assembly(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase881_substituted_terphenyl_falls_through():
    # substituted-ring terphenyl is out of scope for this phase; must not
    # crash and must not silently claim a wrong terphenyl name.
    result = smiles_to_iupac("Cc1ccc(-c2ccc(-c3ccccc3)cc2)cc1")
    assert "terphenyl" not in result
