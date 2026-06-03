"""Phase 368: Imine locant omission for 2-carbon chains (IUPAC 2013 P-62.3).

Previously 'ethan-1-imine' was given instead of 'ethanimine' for CC=N.
The locant 1 is unambiguous for a 2-carbon chain (the imine must be at C1)
so it is omitted, mirroring the amine rule (ethan-1-amine → ethanamine).

This fix applies to both chain imines (via name_assembler) and nitrone
imine-name construction (via _name_nitrone in group_namers).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 2C imine: locant omitted
    ("CC=N",           "ethanimine"),
    ("CC=NC",          "N-methylethanimine"),
    ("CC=NCC",         "N-ethylethanimine"),
    # 1C imine unchanged
    ("C=N",            "methanimine"),
    # 3C+ imines still carry locant
    ("CCC=N",          "propan-1-imine"),
    ("CC(=N)C",        "propan-2-imine"),
    ("CCC=NC",         "N-methylpropan-1-imine"),
    # Nitrone (imine N-oxide): 2C imine also omits locant
    ("CC=[N+]([O-])C", "N-methylethanimine N-oxide"),
    # Regressions: non-imine names unaffected
    ("CC=O",           "acetaldehyde"),
    ("CC(=O)C",        "acetone"),
])
def test_phase368_imine_locant_omission(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
