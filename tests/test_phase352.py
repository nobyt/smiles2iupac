"""Phase 352: E/Z stereo descriptors on N=N bonds (azo compounds, IUPAC 2013).

When an azo compound R-N=N-R carries defined E/Z geometry in the SMILES,
the descriptor is prepended as "(E)-" or "(Z)-".
Unspecified geometry produces no prefix.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Symmetric azo compounds with E/Z
    ("C/N=N/C",                "(E)-azomethane"),
    (r"C/N=N\C",               "(Z)-azomethane"),
    ("CC/N=N/CC",              "(E)-azoethane"),
    (r"CC/N=N\CC",             "(Z)-azoethane"),
    ("CCC/N=N/CCC",            "(E)-azopropane"),
    (r"CCC/N=N\CCC",           "(Z)-azopropane"),
    ("CCCC/N=N/CCCC",          "(E)-azobutane"),
    (r"CCCC/N=N\CCCC",         "(Z)-azobutane"),
    # Unspecified geometry: no prefix
    ("CN=NC",                  "azomethane"),
    ("CCN=NCC",                "azoethane"),
    # Regression: azobenzene unchanged
    ("c1ccc(N=Nc2ccccc2)cc1",  "azobenzene"),
])
def test_phase352_ez_azo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
