"""Phase 874: substituted biphenyl PIN (X-...-1,1'-biphenyl).

The unsubstituted parent already gave the PIN "1,1'-biphenyl", but substituted
biphenyls fell back to substituted-benzene naming ("1-methyl-4-phenylbenzene").
Per IUPAC 2013 P-28.2.1 the ring assembly of two benzene rings uses biphenyl as
the parent hydride, so the PIN is e.g. 4-methyl-1,1'-biphenyl.

Fix: in _build_ring_name, when the benzene parent carries exactly one
(unsubstituted) phenyl substituent plus other substituents, renumber so the
phenyl-bearing carbon is position 1 (choosing the direction giving the other
substituents the lowest locants) and emit "{prefix}-1,1'-biphenyl". If the
other ring is itself substituted, its substituent is named "(4-methylphenyl)"
(not bare "phenyl") and this branch is skipped -- both-ring-substituted
assemblies keep their existing naming.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # single substituent, various positions
    ("Cc1ccc(-c2ccccc2)cc1",   "4-methyl-1,1'-biphenyl"),
    ("Cc1ccccc1-c1ccccc1",     "2-methyl-1,1'-biphenyl"),
    ("Cc1cccc(-c2ccccc2)c1",   "3-methyl-1,1'-biphenyl"),
    ("Clc1ccc(-c2ccccc2)cc1",  "4-chloro-1,1'-biphenyl"),
    ("Brc1ccc(-c2ccccc2)cc1",  "4-bromo-1,1'-biphenyl"),
    ("Fc1ccc(-c2ccccc2)cc1",   "4-fluoro-1,1'-biphenyl"),
    # substituent on the ring written second (parent selection)
    ("c1ccc(-c2ccc(C)cc2)cc1", "4-methyl-1,1'-biphenyl"),
    # multiple substituents on one ring
    ("Cc1cc(C)cc(-c2ccccc2)c1", "3,5-dimethyl-1,1'-biphenyl"),
    ("Clc1ccc(-c2ccccc2)cc1C", "4-chloro-3-methyl-1,1'-biphenyl"),
    # regression: unsubstituted parent unchanged
    ("c1ccc(-c2ccccc2)cc1",    "1,1'-biphenyl"),
    # regression: not biphenyls
    ("Cc1ccccc1C",             "1,2-dimethylbenzene"),
    ("Cc1ccc(C)cc1",           "1,4-dimethylbenzene"),
    ("Clc1ccccc1",             "chlorobenzene"),
    # regression: hetero ring assembly & suffix biphenyl keep their naming
    ("c1ccc(-c2ccncc2)cc1",    "4-phenylpyridine"),
    ("OC(=O)c1ccc(-c2ccccc2)cc1", "[1,1'-biphenyl]-4-carboxylic acid"),  # Phase 876
    # regression: linear terphenyl now gets its PIN (Phase 881, was the
    # non-PIN "1,4-diphenylbenzene" until then)
    ("c1ccc(-c2ccc(-c3ccccc3)cc2)cc1", "1,1':4',1''-terphenyl"),
])
def test_phase874_substituted_biphenyl_pin(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
