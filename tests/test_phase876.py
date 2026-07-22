"""Phase 876: biphenyl carbo-suffix PIN ([1,1'-biphenyl]-X-carboxylic acid etc.).

A carboxylic acid / aldehyde / nitrile on a biphenyl was named as a substituted
benzene ("4-phenylbenzoic acid"). Per IUPAC 2013 P-28.2.1 the ring assembly of
two benzene rings is the biphenyl parent, so the PIN uses the bracketed
biphenyl name with an exocyclic suffix:

  OC(=O)c1ccc(-c2ccccc2)cc1 -> [1,1'-biphenyl]-4-carboxylic acid
  O=Cc1ccc(-c2ccccc2)cc1    -> [1,1'-biphenyl]-4-carbaldehyde
  N#Cc1ccc(-c2ccccc2)cc1    -> [1,1'-biphenyl]-4-carbonitrile

The suffix ring is numbered with its inter-ring carbon at position 1 and the
suffix takes the lowest locant. Scoped to the carbo-family (COOH/CHO/CN) with
no other substituents; biphenyls that also bear substituents, and other group
types (ol/amine/ketone), fall through to their existing naming.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # carboxylic acid at para / meta / ortho
    ("OC(=O)c1ccc(-c2ccccc2)cc1", "[1,1'-biphenyl]-4-carboxylic acid"),
    ("OC(=O)c1cccc(-c2ccccc2)c1", "[1,1'-biphenyl]-3-carboxylic acid"),
    ("OC(=O)c1ccccc1-c1ccccc1",   "[1,1'-biphenyl]-2-carboxylic acid"),
    # aldehyde / nitrile
    ("O=Cc1ccc(-c2ccccc2)cc1",    "[1,1'-biphenyl]-4-carbaldehyde"),
    ("N#Cc1ccc(-c2ccccc2)cc1",    "[1,1'-biphenyl]-4-carbonitrile"),
    ("O=Cc1ccccc1-c1ccccc1",      "[1,1'-biphenyl]-2-carbaldehyde"),
    # regression: plain benzene carbo-suffixes (retained names) unchanged
    ("OC(=O)c1ccccc1", "benzoic acid"),
    ("O=Cc1ccccc1",    "benzaldehyde"),
    ("N#Cc1ccccc1",    "benzonitrile"),
    ("OC(=O)c1ccc(C)cc1", "4-methylbenzoic acid"),
    ("OC(=O)CCC(=O)O", "butanedioic acid"),
    # regression: alkane biphenyls unchanged (Phase 874/875)
    ("c1ccccc1-c1ccccc1",       "1,1'-biphenyl"),
    ("Cc1ccc(-c2ccc(C)cc2)cc1", "4,4'-dimethyl-1,1'-biphenyl"),
])
def test_phase876_biphenyl_carbo_suffix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
