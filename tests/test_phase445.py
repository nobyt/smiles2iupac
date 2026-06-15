"""Phase 445: phenanthridine, benzo[h]isoquinoline, phenanthroline isomers,
benzo[c/f/g]cinnoline retained names (IUPAC 2013 P-31.1.3).
Note: 4,7-phenanthroline entry corrected in Phase 456 to pyrido[3,4-g]quinoline.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phenanthridine (= benzo[f]isoquinoline)
    ("c1ccc2c(c1)ccc1cccnc12",   "phenanthridine"),
    # benzo[h]isoquinoline
    ("c1ccc2c(c1)ccc1ccncc12",   "benzo[h]isoquinoline"),
    # 1,10-phenanthroline
    ("c1cnc2c(c1)ccc1cccnc12",   "1,10-phenanthroline"),
    # 1,8-phenanthroline
    ("c1cnc2c(c1)ccc1cnccc12",   "1,8-phenanthroline"),
    # 1,7-phenanthroline
    ("c1cnc2c(c1)ccc1ncccc12",   "1,7-phenanthroline"),
    # pyrido[3,4-g]quinoline (linear; was wrongly called 4,7-phenanthroline)
    ("c1cnc2cc3ccncc3cc2c1",     "pyrido[3,4-g]quinoline"),
    # benzo[b][1,7]naphthyridine (was mislabeled benzo[c]cinnoline in earlier session)
    ("c1ccc2nc3cnccc3cc2c1",     "benzo[b][1,7]naphthyridine"),
    # benzo[g]cinnoline
    ("c1ccc2cc3nnccc3cc2c1",     "benzo[g]cinnoline"),
    # benzo[f]cinnoline
    ("c1ccc2c(c1)ccc1ccnnc12",   "benzo[f]cinnoline"),
    # regressions
    ("c1ccc2nc3ccccc3cc2c1",     "acridine"),
    ("c1ccc2c(c1)ccc1ncccc12",   "benzo[h]quinoline"),
    ("c1ccc2cc3ncccc3cc2c1",     "benzo[f]quinoline"),
    ("c1ccc2cc3cnccc3cc2c1",     "benzo[f]isoquinoline"),
    ("c1ccc2c(c1)ccc1cnccc12",   "benzo[g]isoquinoline"),
])
def test_phase445_benzo_tricyclics(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
