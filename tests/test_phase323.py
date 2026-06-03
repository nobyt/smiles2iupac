"""Phase 323: N-substituted imidate ester naming (IUPAC 2013 P-65.1.2.3).

Imidate esters with N-substituents are named '{alkyl} N-{sub}{stem}imidoate'.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CCOC(=NC)CC",    "ethyl N-methylpropanimidoate"),
    ("COC(=NC)C",      "methyl N-methylethanimidoate"),
    ("COC(=NCC)C",     "methyl N-ethylethanimidoate"),
    ("CCOC(=NCC)CC",   "ethyl N-ethylpropanimidoate"),
    # regressions: unsubstituted imidoate esters unchanged
    ("CC(=N)OCC",      "ethyl ethanimidoate"),
    ("CCC(=N)OC",      "methyl propanimidoate"),
    ("C(=N)OC",        "methyl methanimidoate"),
])
def test_phase323_n_sub_imidate_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
