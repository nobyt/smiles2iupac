"""Phase 432: Pyrido[b]indole retained names — β- and α-carboline isomers (IUPAC 2013 P-31.1.3).

9H-pyrido[3,4-b]indole (β-carboline) and 9H-pyrido[2,3-b]indole (α-carboline)
are C11H8N2 tricyclic pyridoindoles that currently output 'benzene' (wrong).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 9H-pyrido[3,4-b]indole (beta-carboline) — N at beta position of tricyclic
    ("c1ccc2c(c1)[nH]c1ccncc12",       "9H-pyrido[3,4-b]indole"),
    # 9H-pyrido[2,3-b]indole (alpha-carboline) — N at alpha position
    ("c1ccc2c(c1)[nH]c1cccnc12",       "9H-pyrido[2,3-b]indole"),
    # 9H-pyrido[4,3-b]indole (delta-carboline) — N adjacent to indole junction
    ("c1ccc2c(c1)[nH]c1ncccc12",       "9H-pyrido[4,3-b]indole"),
    # 9H-pyrido[3,4-c]indole (gamma-carboline) — N at gamma position
    ("c1ccc2c(c1)[nH]c1cnccc12",       "9H-pyrido[3,4-c]indole"),
    # regression: 1H-indole unchanged
    ("c1ccc2[nH]ccc2c1",               "1H-indole"),
    # regression: quinoline unchanged
    ("c1ccc2ncccc2c1",                 "quinoline"),
    # regression: 9H-carbazole unchanged (Phase 134)
    ("c1ccc2[nH]c3ccccc3c2c1",         "9H-carbazole"),
])
def test_phase432_pyridoindole(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
