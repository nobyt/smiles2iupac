"""Phase 254: acetyl/formyl retained names for all halides (IUPAC 2013 P-65.1.1.4).

Previously only 'acetyl chloride' and 'formyl chloride' were retained;
now 'acetyl fluoride', 'acetyl bromide', 'acetyl iodide' etc. also use
the retained acyl name prefix.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # acetyl halides (2-carbon acid chain)
    ("CC(=O)F",      "acetyl fluoride"),
    ("CC(=O)Cl",     "acetyl chloride"),
    ("CC(=O)Br",     "acetyl bromide"),
    ("CC(=O)I",      "acetyl iodide"),
    # formyl halides (1-carbon acid chain)
    ("O=CF",         "formyl fluoride"),
    ("O=CCl",        "formyl chloride"),
    ("O=CBr",        "formyl bromide"),
    # regression: longer chains still use systematic names
    ("CCC(=O)Cl",    "propanoyl chloride"),
    ("CCC(=O)F",     "propanoyl fluoride"),
    ("CCCC(=O)Cl",   "butanoyl chloride"),
])
def test_phase254_acyl_halide_retained(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
