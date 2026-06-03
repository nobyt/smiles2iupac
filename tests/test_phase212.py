"""Phase 212: amidoxime naming (IUPAC 2013 P-66.4.1.2)

  CC(=N)NO    → N-hydroxyethanimidamide
  C(=N)NO     → N-hydroxymethanimidamide
  CCC(=N)NO   → N-hydroxypropanimidamide

Amidoximes R-C(=NH)-NHOH have the -OH on the amine N (single bond).
In amidine nomenclature, this amine N is 'N' (unprimed).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Amidoximes
    ("C(=N)NO",    "N-hydroxymethanimidamide"),
    ("CC(=N)NO",   "N-hydroxyethanimidamide"),
    ("CCC(=N)NO",  "N-hydroxypropanimidamide"),
    # regression: N-substituted hydroxylamine still works
    ("CNO",        "N-methylhydroxylamine"),
    ("CCNO",       "N-ethylhydroxylamine"),
    # regression: hydroxamic acid still works
    ("CC(=O)NO",   "N-hydroxyacetamide"),
])
def test_phase212_amidoxime(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
