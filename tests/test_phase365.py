"""Phase 365: Multi-heteroatom monocyclic ring naming (7-10 membered).

Previously _match_hantzsch_widman returned None for rings with 2+ heteroatoms,
causing 7-10 membered heterocycles to fall back to their carbocyclic names
(e.g. 'cycloheptane' instead of '1,4-dioxepane').

Added _match_multi_het_ring which applies a-nomenclature: locants + multiplied
prefix ('ox'/'az'/'thi') + ring suffix ('epane'/'ocane'/'onane'/'ecane').
Multiple heteroatom types are separated with '-' (e.g. '1-oxa-4-azepane').
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 7-membered same-type heteroatoms
    ("C1COCCOC1",   "1,4-dioxepane"),
    ("C1CNCCNC1",   "1,4-diazepane"),
    ("C1CSCCSC1",   "1,4-dithiepane"),
    # 7-membered mixed heteroatoms (O cited first per priority)
    ("C1CNCCOC1",   "1-oxa-4-azepane"),
    ("O1CCNCCC1",   "1-oxa-4-azepane"),
    # 8-membered
    ("C1CCOCCOC1",  "1,4-dioxocane"),
    # 9-membered
    ("C1CCNCCNCC1", "1,4-diazonane"),
    # Substituted 7-membered diaza ring
    ("CC1CNCCNC1",  "6-methyl-1,4-diazepane"),
    # Regressions: smaller rings with retained names unchanged
    ("C1CNCCO1",    "morpholine"),
    ("C1COCCO1",    "1,4-dioxane"),
    ("C1CCNCC1",    "piperidine"),
    ("C1CCOC1",     "oxolane"),
    ("C1CCNC1",     "pyrrolidine"),
])
def test_phase365_multi_het_ring(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
