"""Phase 518: phosphanium / sulfonium / arsonium onium cation naming
(IUPAC 2013 P-73.2, P-73.3, P-73.4: cationic P, S, As with organic substituents).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phosphanium [P+] with 4 C substituents
    ("C[P+](C)(C)C",    "tetramethylphosphanium"),
    ("CC[P+](CC)(CC)CC","tetraethylphosphanium"),
    # sulfonium [S+] with 3 C substituents
    ("C[S+](C)C",       "trimethylsulfonium"),
    ("CC[S+](CC)CC",    "triethylsulfonium"),
    # arsonium [As+] with 4 C substituents
    ("[As+](C)(C)(C)C", "tetramethylarsonium"),
    # regression: ammonium [N+] still works
    ("C[N+](C)(C)C",    "tetramethylazanium"),
    # regression: neutral phosphane still works
    ("CP(C)C",          "trimethylphosphane"),
    # regression: thioether still works
    ("CSC",             "dimethyl sulfide"),
])
def test_phase518_onium_cations(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
