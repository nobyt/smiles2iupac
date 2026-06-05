"""Phase 448: benzo[b][1,5/1,6/1,8]naphthyridine, 1H-naphtho[2,3-d]imidazole
(IUPAC 2013 P-31.1.2 fusion nomenclature).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[b][1,8]naphthyridine (both N adjacent to same ring junction)
    ("c1ccc2nc3ncccc3cc2c1",    "benzo[b][1,8]naphthyridine"),
    # benzo[b][1,5]naphthyridine (N in each ring, each adjacent to different junction)
    ("c1ccc2nc3cccnc3cc2c1",    "benzo[b][1,5]naphthyridine"),
    # benzo[b][1,6]naphthyridine (one N adjacent to junction, other not)
    ("c1ccc2nc3ccncc3cc2c1",    "benzo[b][1,6]naphthyridine"),
    # 1H-naphtho[2,3-d]imidazole (naphthalene C2-C3 fused to imidazole bond d)
    ("c1ccc2cc3[nH]cnc3cc2c1",  "1H-naphtho[2,3-d]imidazole"),
    # regressions
    ("c1ccc2c(c1)ccc1nnccc12",  "benzo[h]cinnoline"),
    ("c1ccc2nc3ccoc3cc2c1",     "furo[3,2-b]quinoline"),
    ("c1ccc2nc3cocc3cc2c1",     "furo[3,4-b]quinoline"),
    ("c1ccc2nc3cnccc3cc2c1",    "benzo[b][1,7]naphthyridine"),
])
def test_phase448_benzo_naphthyridine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
