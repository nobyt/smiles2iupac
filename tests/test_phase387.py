"""Phase 387: Nitroso group on cyclic compounds (IUPAC 2013 P-61.7).

O=N-C(ring) should be named as nitrosocyclohexane, not as a chain compound.

Bug: O=NC1CCCCC1 was giving '6-nitrosoundecane' instead of 'nitrosocyclohexane'.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("O=NC1CCCCC1",   "nitrosocyclohexane"),
    ("O=NC1CCCC1",    "nitrosocyclopentane"),
    ("O=NC1CCCCCC1",  "nitrosocycloheptane"),
    ("O=NC1CCC1",     "nitrosocyclobutane"),
    # nitroso on benzene
    ("O=Nc1ccccc1",   "nitrosobenzene"),
    # regression: nitroso on chain still works
    ("O=NCC",         "nitrosoethane"),
    ("O=NCCC",        "1-nitrosopropane"),
    ("CNO",           "N-methylhydroxylamine"),
])
def test_phase387_nitroso_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
