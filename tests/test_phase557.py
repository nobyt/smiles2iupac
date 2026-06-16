"""Phase 557: 1,x-naphthyridine naming (bicyclic diazines, 10 atoms).
1,5- and 1,8-naphthyridine have C2 symmetry so C2≡C7 etc.; only lower locants used.
1,6- and 1,7-naphthyridine are asymmetric; all 6 C positions are distinct.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1,5-naphthyridine (C2 symmetric: C2≡C8, C3≡C7, C4≡C6)
    ("c1cnc2cccnc2c1",      "1,5-naphthyridine"),
    ("Cc1ccc2ncccc2n1",     "2-methyl-1,5-naphthyridine"),
    ("Cc1cnc2cccnc2c1",     "3-methyl-1,5-naphthyridine"),
    ("Cc1ccnc2cccnc12",     "4-methyl-1,5-naphthyridine"),
    # 1,6-naphthyridine (asymmetric; all 6 C positions distinct)
    ("c1cnc2ccncc2c1",      "1,6-naphthyridine"),
    ("Cc1ccc2cnccc2n1",     "2-methyl-1,6-naphthyridine"),
    ("Cc1cnc2ccncc2c1",     "3-methyl-1,6-naphthyridine"),
    ("Cc1ccnc2ccncc12",     "4-methyl-1,6-naphthyridine"),
    ("Cc1nccc2ncccc12",     "5-methyl-1,6-naphthyridine"),
    ("Cc1cc2ncccc2cn1",     "7-methyl-1,6-naphthyridine"),
    ("Cc1cncc2cccnc12",     "8-methyl-1,6-naphthyridine"),
    # 1,7-naphthyridine (asymmetric; all 6 C positions distinct)
    ("c1cnc2cnccc2c1",      "1,7-naphthyridine"),
    ("Cc1ccc2ccncc2n1",     "2-methyl-1,7-naphthyridine"),
    ("Cc1cnc2cnccc2c1",     "3-methyl-1,7-naphthyridine"),
    ("Cc1ccnc2cnccc12",     "4-methyl-1,7-naphthyridine"),
    ("Cc1cncc2ncccc12",     "5-methyl-1,7-naphthyridine"),
    ("Cc1cc2cccnc2cn1",     "6-methyl-1,7-naphthyridine"),
    ("Cc1nccc2cccnc12",     "8-methyl-1,7-naphthyridine"),
    # 1,8-naphthyridine (C2 symmetric: C2≡C7, C3≡C6, C4≡C5)
    ("c1cnc2ncccc2c1",      "1,8-naphthyridine"),
    ("Cc1ccc2cccnc2n1",     "2-methyl-1,8-naphthyridine"),
    ("Cc1cnc2ncccc2c1",     "3-methyl-1,8-naphthyridine"),
    ("Cc1ccnc2ncccc12",     "4-methyl-1,8-naphthyridine"),
])
def test_phase557_naphthyridines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
