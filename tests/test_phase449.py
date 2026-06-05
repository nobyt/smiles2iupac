"""Phase 449: fix benzo[c]cinnoline (was mislabeled), benzo[b][1,7]naphthyridine,
1H-naphtho[2,3-d]pyrazole, naphtho[2,3-d]oxazole, naphtho[2,3-d]thiazole.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # benzo[c]cinnoline (angular, N-N bond in middle ring)
    ("c1ccc2c(c1)nnc1ccccc12",  "benzo[c]cinnoline"),
    # 1H-naphtho[2,3-d]pyrazole (N-H at pyrazole position 1)
    ("c1ccc2cc3[nH]ncc3cc2c1",  "1H-naphtho[2,3-d]pyrazole"),
    # naphtho[2,3-d]oxazole
    ("c1ccc2cc3ocnc3cc2c1",     "naphtho[2,3-d]oxazole"),
    # naphtho[2,3-d]thiazole
    ("c1ccc2cc3scnc3cc2c1",     "naphtho[2,3-d]thiazole"),
    # regressions: benzo[b][1,7]naphthyridine (previously mislabeled)
    ("c1ccc2nc3cnccc3cc2c1",    "benzo[b][1,7]naphthyridine"),
    ("c1ccc2nc3ncccc3cc2c1",    "benzo[b][1,8]naphthyridine"),
    ("c1ccc2nc3cccnc3cc2c1",    "benzo[b][1,5]naphthyridine"),
    ("c1ccc2nc3ccncc3cc2c1",    "benzo[b][1,6]naphthyridine"),
    ("c1ccc2cc3[nH]cnc3cc2c1",  "1H-naphtho[2,3-d]imidazole"),
    ("c1ccc2cc3nnccc3cc2c1",    "benzo[g]cinnoline"),
    ("c1ccc2c(c1)ccc1ccnnc12",  "benzo[f]cinnoline"),
])
def test_phase449_naphtho_fused(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
