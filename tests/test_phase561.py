"""Phase 561: Substituted 4,7-, 1,7-, and 1,8-phenanthroline naming.
4,7-phenanthroline has C2 symmetry → 4 unique positions (locants 1,2,3,5).
1,7- and 1,8-phenanthroline are asymmetric → 8 unique C positions each.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 4,7-phenanthroline (C2 symmetric: C1≡C10, C2≡C9, C3≡C8, C5≡C6)
    ("c1cnc2ccc3ncccc3c2c1",    "4,7-phenanthroline"),
    ("Cc1cnc2ccc3ncccc3c2c1",   "2-methyl-4,7-phenanthroline"),
    ("Cc1ccc2c(ccc3ncccc32)n1", "3-methyl-4,7-phenanthroline"),
    ("Cc1cc2ncccc2c2cccnc12",   "5-methyl-4,7-phenanthroline"),
    ("Cc1ccnc2ccc3ncccc3c12",   "1-methyl-4,7-phenanthroline"),
    # 1,7-phenanthroline (no symmetry; 8 unique C positions)
    ("c1cnc2c(c1)ccc1ncccc12",  "1,7-phenanthroline"),
    ("Cc1ccc2ccc3ncccc3c2n1",   "2-methyl-1,7-phenanthroline"),
    ("Cc1cnc2c(ccc3ncccc32)c1", "3-methyl-1,7-phenanthroline"),
    ("Cc1ccnc2c1ccc1ncccc12",   "4-methyl-1,7-phenanthroline"),
    ("Cc1cc2ncccc2c2ncccc12",   "5-methyl-1,7-phenanthroline"),
    ("Cc1cc2cccnc2c2cccnc12",   "6-methyl-1,7-phenanthroline"),
    ("Cc1ccc2c(ccc3cccnc32)n1", "8-methyl-1,7-phenanthroline"),
    ("Cc1cnc2ccc3cccnc3c2c1",   "9-methyl-1,7-phenanthroline"),
    ("Cc1ccnc2ccc3cccnc3c12",   "10-methyl-1,7-phenanthroline"),
    # 1,8-phenanthroline (no symmetry; 8 unique C positions)
    ("c1cnc2c(c1)ccc1cnccc12",  "1,8-phenanthroline"),
    ("Cc1ccc2ccc3cnccc3c2n1",   "2-methyl-1,8-phenanthroline"),
    ("Cc1cnc2c(ccc3cnccc32)c1", "3-methyl-1,8-phenanthroline"),
    ("Cc1ccnc2c1ccc1cnccc12",   "4-methyl-1,8-phenanthroline"),
    ("Cc1cc2cnccc2c2ncccc12",   "5-methyl-1,8-phenanthroline"),
    ("Cc1cc2cccnc2c2ccncc12",   "6-methyl-1,8-phenanthroline"),
    ("Cc1nccc2c1ccc1cccnc12",   "7-methyl-1,8-phenanthroline"),
    ("Cc1cc2c(ccc3cccnc32)cn1", "9-methyl-1,8-phenanthroline"),
    ("Cc1cncc2ccc3cccnc3c12",   "10-methyl-1,8-phenanthroline"),
])
def test_phase561_phenanthrolines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
