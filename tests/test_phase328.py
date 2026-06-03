"""Phase 328: E/Z in diester names; dioic acid locant optimization (IUPAC 2013 P-44/P-93.5).

Diesters of alpha,beta-unsaturated diacids carry E/Z prefix.
Dioic acids with double bonds now get the lowest possible locant for the double bond.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # diester E/Z
    ("COC(=O)/C=C/C(=O)OC",         "(2E)-dimethyl but-2-enedioate"),
    (r"COC(=O)/C=C\C(=O)OC",        "(2Z)-dimethyl but-2-enedioate"),
    ("CCOC(=O)/C=C/C(=O)OCC",       "(2E)-diethyl but-2-enedioate"),
    # dioic acid locant optimization: double bond gets lowest locant
    ("OC(=O)CC=CC(=O)O",            "pent-2-enedioic acid"),
    ("OC(=O)C=CCC(=O)O",            "pent-2-enedioic acid"),
    # regressions: saturated diesters use retained names
    ("COC(=O)CCC(=O)OC",            "dimethyl succinate"),
    ("CCOC(=O)CC(=O)OCC",           "diethyl malonate"),
    # regressions: benzene diester unchanged
    ("COC(=O)c1ccc(C(=O)OC)cc1",   "dimethyl benzene-1,4-dicarboxylate"),
])
def test_phase328_ez_diester_dioic_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
