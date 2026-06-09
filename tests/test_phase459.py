"""Phase 459: benzo[c][1,6]naphthyridine, benzo[c][1,7]naphthyridine,
benzo[g][1,8]naphthyridine (IUPAC 2013 P-31.1.3 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzene fused at bond c (C3-C4) of [1,6]naphthyridine ring A
    ("c1ccc2c(c1)cnc1ccncc12",    "benzo[c][1,6]naphthyridine"),
    # benzene fused at bond c (C3-C4) of [1,7]naphthyridine ring A
    ("c1ccc2c(c1)cnc1cnccc12",    "benzo[c][1,7]naphthyridine"),
    # benzene fused at bond g (C6-C7) of [1,8]naphthyridine ring B
    ("c1ccc2c(c1)ncc1ncccc12",    "benzo[g][1,8]naphthyridine"),
])
def test_phase459_benzo_naphthyridines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
