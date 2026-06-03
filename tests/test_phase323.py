"""Phase 323: N-substituted imidate ester naming (IUPAC 2013 P-65.1.2.3).

Imidate esters with N-substituents are named '{alkyl} N-{sub}{stem}imidate'.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CCOC(=NC)CC",    "ethyl N-methylpropanimidate"),
    ("COC(=NC)C",      "methyl N-methylethanimidate"),
    ("COC(=NCC)C",     "methyl N-ethylethanimidate"),
    ("CCOC(=NCC)CC",   "ethyl N-ethylpropanimidate"),
    # regressions: unsubstituted imidate esters unchanged
    ("CC(=N)OCC",      "ethyl ethanimidate"),
    ("CCC(=N)OC",      "methyl propanimidate"),
    ("C(=N)OC",        "methyl methanimidate"),
])
def test_phase323_n_sub_imidate_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
