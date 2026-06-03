"""Phase 371: Selenide/telluride/diselenide preferred IUPAC 2013 names.

Previously _name_selenide_telluride and _name_diselenide_ditelluride used
substitutive forms like '(methylselanyl)methane'. IUPAC 2013 P-66.x requires
the functional-class form analogous to sulfide (Phase 370): group names +
'selenide'/'telluride'/'diselenide'/'ditelluride'.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Selenides
    ("C[Se]C",               "dimethyl selenide"),
    ("CC[Se]CC",             "diethyl selenide"),
    ("C[Se]CC",              "ethyl methyl selenide"),
    # Tellurides
    ("C[Te]C",               "dimethyl telluride"),
    # Diselenides
    ("C[Se][Se]C",           "dimethyl diselenide"),
    ("CC[Se][Se]CC",         "diethyl diselenide"),
    # Ditellurides
    ("C[Te][Te]C",           "dimethyl ditelluride"),
    # E/Z groups → brackets
    ("C/C=C/C[Se]C",         "[(2E)-but-2-en-1-yl] methyl selenide"),
    ("C/C=C/C[Te]C",         "[(2E)-but-2-en-1-yl] methyl telluride"),
    # Regressions: selenol/tellurole unchanged
    ("CC[SeH]",              "ethaneselenol"),
])
def test_phase371_selenide_telluride_preferred(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
