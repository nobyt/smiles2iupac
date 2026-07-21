"""Phase 872: hetero-ring assemblies (bipyridine etc.) + heteroaryl substituent
lowest-locant fix.

Two linked defects:

1. name_heterocycle bailed out (returned None) whenever there was more than one
   hetero ring, so a non-fused assembly of two pyridine rings fell through to
   the carbocycle path, which named the parent pyridine as "benzene" -> the
   name "(pyridin-4-yl)benzene" referenced a benzene ring absent from the
   molecule. Now a non-fused assembly of identical hetero rings is named with
   one ring as parent and the rest as heteroaryl substituents.

2. name_substituent numbered a heteroaryl substituent's attachment point from a
   single fixed rotation, so meta/ortho attachments came out as pyridin-5-yl /
   pyridin-6-yl instead of the lowest-locant pyridin-3-yl / pyridin-2-yl. The
   attachment locant now also considers the mirror numbering that preserves the
   heteroatom locants.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # hetero-ring assemblies: no phantom benzene, correct locants
    ("c1ccncc1-c1ccncc1", "3-(pyridin-4-yl)pyridine"),
    ("c1cccnc1-c1cccnc1", "2-(pyridin-3-yl)pyridine"),
    # heteroaryl substituent lowest-locant (parent = benzoic acid)
    ("OC(=O)c1ccc(-c2cccnc2)cc1", "4-(pyridin-3-yl)benzoic acid"),
    ("OC(=O)c1ccc(-c2ccncc2)cc1", "4-(pyridin-4-yl)benzoic acid"),
    ("OC(=O)c1ccc(-c2ccccn2)cc1", "4-(pyridin-2-yl)benzoic acid"),
    ("OC(=O)c1ccc(-c2ccco2)cc1",  "4-(furan-2-yl)benzoic acid"),
    ("OC(=O)c1ccc(-c2cccs2)cc1",  "4-(thiophen-2-yl)benzoic acid"),
    # regression: no phantom benzene anywhere in these
    ("c1ccncc1-c1ccccc1", "3-phenylpyridine"),
    ("c1ccc(-c2ccncc2)cc1", "4-phenylpyridine"),
    # regression: single heterocycles and benzene biphenyl unchanged
    ("c1ccncc1", "pyridine"),
    ("Cc1ccncc1", "4-methylpyridine"),
    ("Cc1cccnc1", "3-methylpyridine"),
    ("c1ccc2ccccc2c1", "naphthalene"),
    ("c1ccc2ncccc2c1", "quinoline"),
    ("c1ccccc1-c1ccccc1", "1,1'-biphenyl"),
    ("c1ccoc1", "furan"),
    ("c1ccncn1", "pyrimidine"),
])
def test_phase872_hetero_ring_assemblies(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase872_no_phantom_benzene():
    # A molecule with no benzene ring must never be named with "benzene".
    for smiles in ("c1ccncc1-c1ccncc1", "c1cccnc1-c1cccnc1"):
        name = smiles_to_iupac(smiles)
        assert "benzene" not in name, f"{smiles} -> {name} (phantom benzene)"
