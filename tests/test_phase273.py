"""Phase 273: substituted bicyclo[l.m.n]alkane and spiro[m.n]alkane (IUPAC 2013 P-31.1.2).

Von Baeyer bicyclic and spiro compounds with substituents.
Numbering rules (P-31.1.2.2 / P-31.1.3.1): choose the orientation
that gives the lowest locant set to substituents.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── unsubstituted (regression) ─────────────────────────────────────────
    ("C1CC2CCC1CC2",         "bicyclo[2.2.2]octane"),
    ("C1CC2CCC1C2",          "bicyclo[2.2.1]heptane"),
    ("C1CC2CCCC1CC2",        "bicyclo[3.2.2]nonane"),
    ("C1CCC2(CC1)CCCC2",     "spiro[4.5]decane"),
    ("C1CCC2(CC1)CCC2",      "spiro[3.5]nonane"),

    # ── substituted bicyclo[2.2.2]octane ───────────────────────────────────
    # methyl on bridge → position 2
    ("CC1CC2CCC1CC2",        "2-methylbicyclo[2.2.2]octane"),
    # methyl at bridgehead → position 1
    ("CC12CCC(CC1)CC2",      "1-methylbicyclo[2.2.2]octane"),
    # bromo on bridge
    ("BrC1CC2CCC1CC2",       "2-bromobicyclo[2.2.2]octane"),
    # chloromethyl substituent
    ("ClCC1CC2CCC1CC2",      "2-(chloromethyl)bicyclo[2.2.2]octane"),

    # ── substituted bicyclo[2.2.1]heptane (norbornane) ─────────────────────
    ("CC1CC2CCC1C2",         "2-methylbicyclo[2.2.1]heptane"),
    ("CC12CCC(C1)CC2",       "1-methylbicyclo[2.2.1]heptane"),
    ("BrC1CC2CCC1C2",        "2-bromobicyclo[2.2.1]heptane"),

    # ── substituted spiro[4.5]decane ────────────────────────────────────────
    # methyl in the larger ring → position 8
    ("CC1CCC2(CC1)CCCC2",    "8-methylspiro[4.5]decane"),
    # gem-dimethyl on the large ring at spiro atom
    ("CC1(C)CCC2(CC1)CCCC2", "8,8-dimethylspiro[4.5]decane"),
])
def test_phase273_substituted_polycyclic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
