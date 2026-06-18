"""Phase 621: R/S stereodescriptor locant-1 omission (IUPAC 2013 P-93.5.2.3).

When there is exactly one R/S stereocenter and its locant is 1, the locant
is omitted from the descriptor.  Multi-center and locant>1 cases keep the
locant.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # single stereocenter at C-1: locant omitted
    ("[C@@H](F)(Cl)Br",           "(R)-bromochlorofluoromethane"),
    ("[C@H](F)(Cl)Br",            "(S)-bromochlorofluoromethane"),
    ("[C@@H](O)(N)C",             "(S)-1-aminoethanol"),
    ("[C@H](O)(N)C",              "(R)-1-aminoethanol"),
    # single stereocenter at C-1, chain length > 1: locant still omitted
    ("[C@@H](F)(Cl)CCl",          "(R)-1,2-dichloro-1-fluoroethane"),
    # single stereocenter at C-2: locant retained (only C-1 is omitted)
    ("C[C@@H](F)CC",              "(2R)-2-fluorobutane"),
    ("C[C@H](F)CC",               "(2S)-2-fluorobutane"),
    # two stereocenters: both locants retained even if one is 1
    ("[C@@H](F)([C@@H](F)Cl)Cl", "(1S,2S)-1,2-dichloro-1,2-difluoroethane"),
    ("[C@@H](F)([C@H](F)Cl)Cl",  "(1S,2R)-1,2-dichloro-1,2-difluoroethane"),
    # E/Z always retains locant (P-93.5.2.2)
    ("C(/C=C/C)CC",               "(2E)-hex-2-ene"),
    ("C(/C=C\\C)CC",              "(2Z)-hex-2-ene"),
])
def test_phase621_rs_locant_omission(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
