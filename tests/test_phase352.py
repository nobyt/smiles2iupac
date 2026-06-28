"""Phase 352: E/Z stereo descriptors on N=N bonds (azo compounds, IUPAC 2013).

When an azo compound R-N=N-R carries defined E/Z geometry in the SMILES,
the descriptor is prepended as "(E)-" or "(Z)-".
Unspecified geometry produces no prefix.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Symmetric diazenes with E/Z
    ("C/N=N/C",                "(E)-dimethyldiazene"),
    (r"C/N=N\C",               "(Z)-dimethyldiazene"),
    ("CC/N=N/CC",              "(E)-diethyldiazene"),
    (r"CC/N=N\CC",             "(Z)-diethyldiazene"),
    ("CCC/N=N/CCC",            "(E)-dipropyldiazene"),
    (r"CCC/N=N\CCC",           "(Z)-dipropyldiazene"),
    ("CCCC/N=N/CCCC",          "(E)-dibutyldiazene"),
    (r"CCCC/N=N\CCCC",         "(Z)-dibutyldiazene"),
    # Unspecified geometry: no prefix
    ("CN=NC",                  "dimethyldiazene"),
    ("CCN=NCC",                "diethyldiazene"),
    # Regression: diphenyldiazene
    ("c1ccc(N=Nc2ccccc2)cc1",  "diphenyldiazene"),
])
def test_phase352_ez_azo(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
