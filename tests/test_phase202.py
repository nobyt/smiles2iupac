"""Phase 202: N-substituted hydroxylamine naming (IUPAC 2013 P-62.3)

  CNO   → N-methylhydroxylamine
  CCNO  → N-ethylhydroxylamine
  CN(O)C → N,N-dimethylhydroxylamine
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # N-monosubstituted hydroxylamines
    ("CNO",    "N-methylhydroxylamine"),
    ("CCNO",   "N-ethylhydroxylamine"),
    ("CCCNO",  "N-propylhydroxylamine"),
    # N,N-disubstituted hydroxylamine
    ("CN(O)C", "N,N-dimethylhydroxylamine"),
    # regression: hydroxylamine itself retained
    ("NO",     "hydroxylamine"),
    # regression: normal amines unaffected
    ("NC",     "methanamine"),
    ("NCC",    "ethanamine"),
])
def test_phase202_n_substituted_hydroxylamine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
