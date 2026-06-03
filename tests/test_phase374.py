"""Phase 374: Thiohydrazide and selenohydrazide naming (IUPAC 2013 P-65.1.4).

RC(=S)-NHNH2 → {stem}anethiohydrazide
RC(=[Se])-NHNH2 → {stem}aneselenohydrazide

Previously C(=S)-NH-NH2 was swallowed by the thioamide detector (which
checked only for C(=S)-N without excluding N-N patterns), so the compound
was named as a thioamide instead.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Thiohydrazides
    ("CC(=S)NN",          "ethanethiohydrazide"),
    ("CCC(=S)NN",         "propanethiohydrazide"),
    ("CCCC(=S)NN",        "butanethiohydrazide"),
    ("C(=S)NN",           "methanethiohydrazide"),
    # Aryl thiohydrazide
    ("c1ccccc1C(=S)NN",   "benzothiohydrazide"),
    # N'-substituted thiohydrazide
    ("CC(=S)NNC",         "N'-methylethanethiohydrazide"),
    # Selenohydrazides
    ("CC(=[Se])NN",       "ethaneselenohydrazide"),
    ("CCC(=[Se])NN",      "propaneselenohydrazide"),
    # Regressions: thioamide unchanged
    ("CC(=S)N",           "ethanethioamide"),
    ("CC(=S)NC",          "N-methylethanethioamide"),
    # Regressions: hydrazide unchanged
    ("CC(=O)NN",          "ethanohydrazide"),
    ("CCC(=O)NN",         "propanohydrazide"),
])
def test_phase374_thiohydrazide_selenohydrazide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
