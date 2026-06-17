"""Phase 562: Substituted 2,7-, 2,6-, 1,6-, 3,6-, 3,5-, and 4,5-phenanthroline naming."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 2,7-phenanthroline (N at 2 and 7; 8 unique C positions: 1,3,4,5,6,8,9,10)
    ("c1cnc2ccc3ccncc3c2c1",    "2,7-phenanthroline"),
    ("Cc1cnc2ccc3ccncc3c2c1",   "10-methyl-2,7-phenanthroline"),
    ("Cc1ccc2c(ccc3ccncc32)n1", "1-methyl-2,7-phenanthroline"),
    ("Cc1cc2ccncc2c2cccnc12",   "3-methyl-2,7-phenanthroline"),
    ("Cc1cc2ncccc2c2cnccc12",   "4-methyl-2,7-phenanthroline"),
    ("Cc1cncc2c1ccc1ncccc12",   "5-methyl-2,7-phenanthroline"),
    ("Cc1cc2ccc3ncccc3c2cn1",   "6-methyl-2,7-phenanthroline"),
    ("Cc1nccc2ccc3ncccc3c12",   "8-methyl-2,7-phenanthroline"),
    ("Cc1ccnc2ccc3ccncc3c12",   "9-methyl-2,7-phenanthroline"),
    # 2,6-phenanthroline (N at 2 and 6; 8 unique C positions: 1,3,4,5,7,8,9,10)
    ("c1ccc2c(c1)ncc1ccncc12",  "2,6-phenanthroline"),
    ("Cc1ccc2c(c1)ncc1ccncc12", "8-methyl-2,6-phenanthroline"),
    ("Cc1ccc2ncc3ccncc3c2c1",   "9-methyl-2,6-phenanthroline"),
    ("Cc1cccc2ncc3ccncc3c12",   "10-methyl-2,6-phenanthroline"),
    ("Cc1cccc2c1ncc1ccncc12",   "7-methyl-2,6-phenanthroline"),
    ("Cc1nc2ccccc2c2cnccc12",   "5-methyl-2,6-phenanthroline"),
    ("Cc1cncc2c1cnc1ccccc12",   "4-methyl-2,6-phenanthroline"),
    ("Cc1cc2cnc3ccccc3c2cn1",   "3-methyl-2,6-phenanthroline"),
    ("Cc1nccc2cnc3ccccc3c12",   "1-methyl-2,6-phenanthroline"),
    # 1,6-phenanthroline (N at 1 and 6; 8 unique C positions: 2,3,4,5,7,8,9,10)
    ("c1cnc2c(c1)cnc1ccccc12",  "1,6-phenanthroline"),
    ("Cc1cnc2c(cnc3ccccc32)c1", "3-methyl-1,6-phenanthroline"),
    ("Cc1ccc2cnc3ccccc3c2n1",   "2-methyl-1,6-phenanthroline"),
    ("Cc1ccnc2c1cnc1ccccc12",   "4-methyl-1,6-phenanthroline"),
    ("Cc1nc2ccccc2c2ncccc12",   "5-methyl-1,6-phenanthroline"),
    ("Cc1cccc2c1ncc1cccnc12",   "7-methyl-1,6-phenanthroline"),
    ("Cc1ccc2c(c1)ncc1cccnc12", "8-methyl-1,6-phenanthroline"),
    ("Cc1ccc2ncc3cccnc3c2c1",   "9-methyl-1,6-phenanthroline"),
    ("Cc1cccc2ncc3cccnc3c12",   "10-methyl-1,6-phenanthroline"),
    # 3,6-phenanthroline (N at 3 and 6; 8 unique C positions: 1,2,4,5,7,8,9,10)
    ("c1ccc2c(c1)ncc1cnccc12",  "3,6-phenanthroline"),
    ("Cc1ccc2c(c1)ncc1cnccc12", "8-methyl-3,6-phenanthroline"),
    ("Cc1ccc2ncc3cnccc3c2c1",   "9-methyl-3,6-phenanthroline"),
    ("Cc1cccc2ncc3cnccc3c12",   "10-methyl-3,6-phenanthroline"),
    ("Cc1cccc2c1ncc1cnccc12",   "7-methyl-3,6-phenanthroline"),
    ("Cc1nc2ccccc2c2ccncc12",   "5-methyl-3,6-phenanthroline"),
    ("Cc1nccc2c1cnc1ccccc12",   "4-methyl-3,6-phenanthroline"),
    ("Cc1cc2c(cn1)cnc1ccccc12", "2-methyl-3,6-phenanthroline"),
    ("Cc1cncc2cnc3ccccc3c12",   "1-methyl-3,6-phenanthroline"),
    # 3,5-phenanthroline (N at 3 and 5; 8 unique C positions: 1,2,4,6,7,8,9,10)
    ("c1ccc2c(c1)cnc1cnccc12",  "3,5-phenanthroline"),
    ("Cc1ccc2c(cnc3cnccc32)c1", "10-methyl-3,5-phenanthroline"),
    ("Cc1ccc2cnc3cnccc3c2c1",   "9-methyl-3,5-phenanthroline"),
    ("Cc1cccc2cnc3cnccc3c12",   "8-methyl-3,5-phenanthroline"),
    ("Cc1cccc2c1cnc1cnccc12",   "1-methyl-3,5-phenanthroline"),
    ("Cc1nc2cnccc2c2ccccc12",   "2-methyl-3,5-phenanthroline"),
    ("Cc1nccc2c1ncc1ccccc12",   "4-methyl-3,5-phenanthroline"),
    ("Cc1cc2c(cn1)ncc1ccccc12", "6-methyl-3,5-phenanthroline"),
    ("Cc1cncc2ncc3ccccc3c12",   "7-methyl-3,5-phenanthroline"),
    # 4,5-phenanthroline (N at 4 and 5; 8 unique C positions: 1,2,3,6,7,8,9,10)
    ("c1ccc2c(c1)cnc1ncccc12",  "4,5-phenanthroline"),
    ("Cc1ccc2c(cnc3ncccc32)c1", "1-methyl-4,5-phenanthroline"),
    ("Cc1ccc2cnc3ncccc3c2c1",   "10-methyl-4,5-phenanthroline"),
    ("Cc1cccc2cnc3ncccc3c12",   "9-methyl-4,5-phenanthroline"),
    ("Cc1cccc2c1cnc1ncccc12",   "2-methyl-4,5-phenanthroline"),
    ("Cc1nc2ncccc2c2ccccc12",   "3-methyl-4,5-phenanthroline"),
    ("Cc1ccc2c(ncc3ccccc32)n1", "6-methyl-4,5-phenanthroline"),
    ("Cc1cnc2ncc3ccccc3c2c1",   "7-methyl-4,5-phenanthroline"),
    ("Cc1ccnc2ncc3ccccc3c12",   "8-methyl-4,5-phenanthroline"),
])
def test_phase562_phenanthrolines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
