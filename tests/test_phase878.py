"""Phase 878: naphthalene as a substituent -> naphthalen-N-yl (stop ring drop).

Same info-loss bug class as Phase 877, but for a fused ring instead of a
linked one. When a naphthalene was a substituent on a chain/principal group,
the second (fused) ring was dropped entirely:
CC(=O)c1ccc2ccccc2c1 (2-acetonaphthone) was named "1-phenylethanone" --
identical to acetophenone (CC(=O)c1ccccc1). Two different compounds collapsed
to one name.

Root cause: the aryl-substituent BFS in _name_carbon_substituent traversed
the ring-fusion bond, collected 10 aromatic carbons, and (len != 6) fell back
to bare "phenyl". Fix: detect the naphthalene-as-substituent (two fused
all-carbon aromatic 6-rings, isolated -- not part of a larger fused system
like anthracene) and emit the naphthalen-N-yl PIN via the existing
_assign_naphthalene_locants numbering (attachment carbon gets the lowest
locant). Only pure (unsubstituted-ring) naphthalene substituents use this
path; substituted-ring / larger fused systems fall through unchanged.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 2-acetonaphthone (beta position)
    ("CC(=O)c1ccc2ccccc2c1",   "1-(naphthalen-2-yl)ethan-1-one"),
    # 1-acetonaphthone (alpha position)
    ("CC(=O)c1cccc2ccccc12",   "1-(naphthalen-1-yl)ethan-1-one"),
    # naphthalen-2-ylmethanol
    ("OCc1ccc2ccccc2c1",       "(naphthalen-2-yl)methanol"),
    # regression: acetophenone (single ring) unchanged
    ("CC(=O)c1ccccc1",         "1-phenylethanone"),
    # regression: naphthalene as PARENT (principal-group anchored) unchanged
    ("c1ccc2ccccc2c1",         "naphthalene"),
    ("Cc1ccc2ccccc2c1",        "2-methylnaphthalene"),
    ("OC(=O)c1ccc2ccccc2c1",   "naphthalene-2-carboxylic acid"),
    # regression: biphenyl (Phase 877) unchanged
    ("CC(=O)c1ccc(-c2ccccc2)cc1", "1-([1,1'-biphenyl]-4-yl)ethan-1-one"),
])
def test_phase878_naphthalenyl_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase878_not_confused_with_acetophenone():
    # 2-acetonaphthone and acetophenone must get distinct names.
    assert smiles_to_iupac("CC(=O)c1ccc2ccccc2c1") != smiles_to_iupac("CC(=O)c1ccccc1")


def test_phase878_substituted_ring_falls_through():
    # Substituted-ring naphthalenyl is out of scope for this phase; it must
    # not crash and must not silently claim a wrong naphthalenyl name.
    result = smiles_to_iupac("CC(=O)c1ccc2ccccc2c1C")
    assert result != "1-(naphthalen-2-yl)ethan-1-one"


def test_phase878_anthracene_substituent_not_misnamed_naphthalenyl():
    # A 3-ring fused substituent must not be (mis)claimed as naphthalenyl.
    result = smiles_to_iupac("CC(=O)c1ccc2cc3ccccc3cc2c1")
    assert "naphthalen" not in result
