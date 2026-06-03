"""Phase 204: nitrone (imine N-oxide) naming (IUPAC 2013 P-62.5.1)

  C=[N+]([O-])C  → N-methylmethanimine N-oxide
  CC=[N+]([O-])C → N-methylethan-1-imine N-oxide
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Nitrone: CH2=N+(CH3)-O-  → N-methylmethanimine N-oxide
    ("C=[N+]([O-])C",    "N-methylmethanimine N-oxide"),
    # Nitrone: CH3CH=N+(CH3)-O- → N-methylethan-1-imine N-oxide
    ("CC=[N+]([O-])C",   "N-methylethan-1-imine N-oxide"),
    # regression: amine N-oxide still works
    ("C[N+](C)(C)[O-]",  "N,N-dimethylmethanamine N-oxide"),
])
def test_phase204_nitrone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
