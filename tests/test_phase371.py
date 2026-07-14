"""Phase 371: Diselenide/ditelluride preferred IUPAC 2013 names.

Diselenides/ditellurides keep the functional-class 'dialkyl diselenide'
form (no simple substitutive equivalent commonly used for -Se-Se-/-Te-Te-
as a substituent chain, same reasoning as disulfide in Phase 370).

Selenides/tellurides (single Se/Te) were reverted to functional-class
naming by this phase at the time, mirroring Phase 370's sulfide mistake --
same fix (Phase 855/856): substitutive '(alkylselanyl)parent' /
'(alkyltellanyl)parent' is the IUPAC 2013 preferred name.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Selenides
    ("C[Se]C",               "methylselanylmethane"),
    ("CC[Se]CC",             "ethylselanylethane"),
    ("C[Se]CC",              "methylselanylethane"),
    # Tellurides
    ("C[Te]C",               "methyltellanylmethane"),
    # Diselenides
    ("C[Se][Se]C",           "dimethyl diselenide"),
    ("CC[Se][Se]CC",         "diethyl diselenide"),
    # Ditellurides
    ("C[Te][Te]C",           "dimethyl ditelluride"),
    # E/Z substituent
    ("C/C=C/C[Se]C",         "(2E)-1-methylselanylbut-2-ene"),
    ("C/C=C/C[Te]C",         "(2E)-1-methyltellanylbut-2-ene"),
    # Regressions: selenol/tellurole unchanged
    ("CC[SeH]",              "ethaneselenol"),
])
def test_phase371_selenide_telluride_preferred(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
