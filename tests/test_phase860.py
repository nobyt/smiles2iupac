"""Phase 860: guard the loosened phosphine/arsane/bismuthane/stibane oxide &
sulfide detectors so they only claim compounds whose central atom carries
*only* C/H substituents besides the one =O / =S.

Phase 859 loosened these detectors from len(c_neighbors) >= 3 to >= 1 so that
primary/secondary oxides/sulfides (RPH2=O, R2PH=O, ...) would be named. But
the only heteroatom guard was `not o_single`, which does not exclude N/Cl/F
substituents. As a result compounds like the phosphonic diamide CP(=O)(N)N or
the phosphonic dichloride CP(=O)(Cl)Cl matched phosphine_oxide and were
mis-named "methylphosphane oxide", silently dropping the amide/halide groups.

Phase 860 adds an "all non-chalcogen neighbours are C or H" guard to each
loosened branch. A real phosphine oxide is R_nP(H)_{3-n}=O with only C/H on P,
so the guard is definitionally correct: legitimate oxides/sulfides still match,
while N/Cl/F-substituted derivatives fall through (they are separate functional
classes -- phosphonamides, phosphonic dihalides, etc. -- not yet modelled, so
falling through is honest rather than a confident wrong name).
"""

import pytest
from smiles2iupac import smiles_to_iupac


# The Phase 859 behaviour that must be preserved (only C/H + one =O/=S).
@pytest.mark.parametrize("smiles,expected", [
    ("CP(C)=O",          "dimethylphosphane oxide"),
    ("CP(C)=S",          "dimethylphosphane sulfide"),
    ("CP=O",             "methylphosphane oxide"),
    ("CP=S",             "methylphosphane sulfide"),
    ("C[AsH](C)=O",      "dimethylarsane oxide"),
    ("C[AsH](C)=S",      "dimethylarsane sulfide"),
    ("C[Bi](C)=O",       "dimethylbismuthane oxide"),
    ("C[Sb](C)=S",       "dimethylstibane sulfide"),
    ("C[P](C)(C)=O",     "trimethylphosphane oxide"),
    ("C[As](C)(C)=S",    "trimethylarsane sulfide"),
    ("C[Bi](C)(C)=O",    "trimethylbismuthane oxide"),
])
def test_phase860_real_oxides_still_named(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


# Hetero-substituted P/As/Bi central atoms must NOT be mis-claimed as a plain
# oxide/sulfide (which would drop the N/Cl/F substituents). We assert the
# negative rather than an exact fall-through string, so the test stays valid if
# these classes later get real names.
@pytest.mark.parametrize("smiles,forbidden", [
    ("C[P](=O)(N)N",     "methylphosphane oxide"),    # phosphonic diamide
    ("C[P](=O)(NC)NC",   "methylphosphane oxide"),
    ("C[P](=O)(Cl)Cl",   "methylphosphane oxide"),    # phosphonic dichloride
    ("C[P](=O)(F)F",     "methylphosphane oxide"),
    ("C[P](=S)(N)N",     "methylphosphane sulfide"),  # phosphonothioic diamide
    ("C[As](=O)(N)N",    "methylarsane oxide"),
    ("C[Bi](C)(=O)N",    "dimethylbismuthane oxide"),
])
def test_phase860_hetero_substituted_not_misclaimed(smiles, forbidden):
    assert smiles_to_iupac(smiles) != forbidden
