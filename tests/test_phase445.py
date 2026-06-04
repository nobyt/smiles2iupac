"""Phase 445: phenanthridine, benzo[h]isoquinoline, phenanthroline isomers,
benzo[c]cinnoline, and benzo[g]phthalazine retained names (IUPAC 2013 P-31.1.3).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phenanthridine (= benzo[f]isoquinoline)
    ("c1ccc2c(c1)ccc1cccnc12",   "phenanthridine"),
    # benzo[h]isoquinoline
    ("c1ccc2c(c1)ccc1ccncc12",   "benzo[h]isoquinoline"),
    # 1,10-phenanthroline
    ("c1cnc2ccc3ncccc3c2c1",     "1,10-phenanthroline"),
    # 1,8-phenanthroline
    ("c1cnc2c(c1)ccc1ncccc12",   "1,8-phenanthroline"),
    # 1,7-phenanthroline
    ("c1cnc2c(c1)cnc1ccccc12",   "1,7-phenanthroline"),
    # 4,7-phenanthroline
    ("c1cnc2cc3ccncc3cc2c1",     "4,7-phenanthroline"),
    # benzo[c]cinnoline
    ("c1ccc2nc3cnccc3cc2c1",     "benzo[c]cinnoline"),
    # benzo[g]phthalazine
    ("c1ccc2cc3nnccc3cc2c1",     "benzo[g]phthalazine"),
    # regressions
    ("c1ccc2nc3ccccc3cc2c1",     "acridine"),
    ("c1ccc2c(c1)ccc1ncccc12",   "benzo[h]quinoline"),
    ("c1ccc2cc3ncccc3cc2c1",     "benzo[f]quinoline"),
    ("c1ccc2cc3cnccc3cc2c1",     "benzo[f]isoquinoline"),
    ("c1ccc2c(c1)ccc1cnccc12",   "benzo[g]isoquinoline"),
])
def test_phase445_benzo_tricyclics(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
