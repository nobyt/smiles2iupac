"""Phase 875: both-ring-substituted biphenyls (X,Y'-...-1,1'-biphenyl).

Phase 874 handled biphenyls substituted on one ring only. When both rings bear
substituents the substituent on the second ring was named "(4-methylphenyl)"
and the molecule fell back to "1-(4-methylphenyl)-4-methylbenzene".

Phase 875 replaces the string-based biphenyl branch with an atom-level namer
(_name_biphenyl_assembly) that handles single- and both-ring cases uniformly:
each ring is numbered with its inter-ring carbon at position 1, substituents on
the two rings get unprimed / primed locants, and the primed/unprimed assignment
plus numbering directions are chosen by the lowest-locant rule (unprimed ranks
below primed, with an alphabetical tiebreak):

  Cc1ccc(-c2ccc(C)cc2)cc1  -> 4,4'-dimethyl-1,1'-biphenyl
  Cc1ccc(-c2cccc(C)c2)cc1  -> 3,4'-dimethyl-1,1'-biphenyl   (meta ring unprimed)
  Clc1ccc(-c2ccc(C)cc2)cc1 -> 4-chloro-4'-methyl-1,1'-biphenyl

Only non-fused all-carbon two-benzene-ring assemblies are handled; hetero
assemblies, fused systems, 3+ ring assemblies and suffix biphenyls fall
through unchanged.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # both rings substituted
    ("Cc1ccc(-c2ccc(C)cc2)cc1",  "4,4'-dimethyl-1,1'-biphenyl"),
    ("Cc1ccc(-c2cccc(C)c2)cc1",  "3,4'-dimethyl-1,1'-biphenyl"),
    ("Cc1cccc(-c2cccc(C)c2)c1",  "3,3'-dimethyl-1,1'-biphenyl"),
    ("Cc1ccccc1-c1ccccc1C",      "2,2'-dimethyl-1,1'-biphenyl"),
    # mixed substituents, both SMILES orderings converge to the canonical name
    ("Clc1ccc(-c2ccc(C)cc2)cc1", "4-chloro-4'-methyl-1,1'-biphenyl"),
    ("Cc1ccc(-c2ccc(Cl)cc2)cc1", "4-chloro-4'-methyl-1,1'-biphenyl"),
    # regression: Phase 874 single-ring cases (now via the same atom-level path)
    ("Cc1ccc(-c2ccccc2)cc1",     "4-methyl-1,1'-biphenyl"),
    ("Cc1ccccc1-c1ccccc1",       "2-methyl-1,1'-biphenyl"),
    ("Clc1ccc(-c2ccccc2)cc1C",   "4-chloro-3-methyl-1,1'-biphenyl"),
    ("c1ccccc1-c1ccccc1",        "1,1'-biphenyl"),
    # regression: non-biphenyl and out-of-scope assemblies unchanged
    ("Cc1ccccc1C",               "1,2-dimethylbenzene"),
    ("c1ccc(-c2ccncc2)cc1",      "4-phenylpyridine"),
    ("OC(=O)c1ccc(-c2ccccc2)cc1", "[1,1'-biphenyl]-4-carboxylic acid"),  # Phase 876
    # Phase 881: linear terphenyl now gets its PIN, was non-PIN before then
    ("c1ccc(-c2ccc(-c3ccccc3)cc2)cc1", "1,1':4',1''-terphenyl"),
    ("c1ccc2ccccc2c1",           "naphthalene"),
])
def test_phase875_both_ring_biphenyl(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
