"""Phase 247: heteroatom ring substituent naming (oxiran-2-yl, thiiran-2-yl, etc.).

When a small heterocyclic ring (oxirane, thiirane, azetidine, etc.) is a substituent,
it is named as {ring_name}-{locant}-yl using the Hantzsch-Widman ring name.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # oxiranyl substituents
    ("OCC1CO1",          "(oxiran-2-yl)methanol"),
    ("OCC1CS1",          "(thiiran-2-yl)methanol"),
    ("OCC1CCO1",         "(oxetan-2-yl)methanol"),
    ("OCC1CCCO1",        "(oxolan-2-yl)methanol"),
    ("NCC1CO1",          "(oxiran-2-yl)methanamine"),
    # regression: 2-substituted oxirane still works
    ("CC1CO1",           "2-methyloxirane"),
    ("ClCC1CO1",         "2-(chloromethyl)oxirane"),
    # regression: oxirane/thiirane as parent unchanged
    ("C1CO1",            "oxirane"),
    ("C1CS1",            "thiirane"),
    # regression: cycloalkyl (all-carbon ring) unchanged
    ("OCC1CC1",          "cyclopropylmethanol"),
    ("OCC1CCCC1",        "cyclopentylmethanol"),
])
def test_phase247_heteroatom_ring_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
