"""Phase 558: 2,6- and 2,7-naphthyridine naming.
Both isomers have C2 symmetry (C1≡C5 in 2,6; C1≡C8 in 2,7), so each has
only 3 unique substitution positions; lower locant (1, 3, 4) is preferred.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 2,6-naphthyridine (C2 symmetric: C1≡C5, C3≡C7, C4≡C8)
    ("c1cc2cnccc2cn1",    "2,6-naphthyridine"),
    ("Cc1nccc2cnccc12",   "1-methyl-2,6-naphthyridine"),
    ("Cc1cc2cnccc2cn1",   "3-methyl-2,6-naphthyridine"),
    ("Cc1cncc2ccncc12",   "4-methyl-2,6-naphthyridine"),
    # 2,7-naphthyridine (C2 symmetric: C1≡C8, C3≡C6, C4≡C5)
    ("c1cc2ccncc2cn1",    "2,7-naphthyridine"),
    ("Cc1nccc2ccncc12",   "1-methyl-2,7-naphthyridine"),
    ("Cc1cc2ccncc2cn1",   "3-methyl-2,7-naphthyridine"),
    ("Cc1cncc2cnccc12",   "4-methyl-2,7-naphthyridine"),
])
def test_phase558_26_27_naphthyridines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
