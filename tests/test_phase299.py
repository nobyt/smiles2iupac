"""Phase 299: geminal di-substituents on methane — locants omitted (IUPAC 2013 P-14.5).

For methane (1-carbon chain), all substituents are at C1. Locants are redundant
and must not be cited: "methanediol" not "methane-1,1-diol".
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("OC(O)",              "methanediol"),
    ("OCO",                "methanediol"),
    ("NCN",                "methanediamine"),
    ("NC(N)",              "methanediamine"),
    ("SC(S)",              "methanedithiol"),
    # regressions: 2+ carbon chains keep locants
    ("OCC(O)",             "ethane-1,2-diol"),
    ("NCC(N)",             "ethane-1,2-diamine"),
    ("OCCO",               "ethane-1,2-diol"),
])
def test_phase299_methane_geminal(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
