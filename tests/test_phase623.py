"""Phase 623: methyl-substituted diazanaphthalenes and acridine/pteridine (IUPAC 2013).

Covers all carbon positions on quinoxaline, quinazoline, phthalazine,
1,5-/1,7-/1,8-/2,6-naphthyridine, acridine, and pteridine.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # quinoxaline (N at 1,4; C2≡C3 symmetric, C5≡C8, C6≡C7)
    ("Cc1ccc2nccnc2c1", "6-methylquinoxaline"),
    ("Cc1cccc2nccnc12", "5-methylquinoxaline"),
    ("Cc1cnc2ccccc2n1", "2-methylquinoxaline"),
    # quinazoline (N at 1,3)
    ("Cc1ccc2ncncc2c1", "6-methylquinazoline"),
    ("Cc1ccc2cncnc2c1", "7-methylquinazoline"),
    ("Cc1cccc2cncnc12", "8-methylquinazoline"),
    ("Cc1ncc2ccccc2n1", "2-methylquinazoline"),
    ("Cc1ncnc2ccccc12", "4-methylquinazoline"),
    ("Cc1cccc2ncncc12", "5-methylquinazoline"),
    # phthalazine (N at 2,3)
    ("Cc1ccc2cnncc2c1", "6-methylphthalazine"),
    ("Cc1cccc2cnncc12", "5-methylphthalazine"),
    ("Cc1nncc2ccccc12", "1-methylphthalazine"),
    # 1,5-naphthyridine
    ("Cc1ccc2ncccc2n1", "2-methyl-1,5-naphthyridine"),
    ("Cc1cnc2cccnc2c1", "3-methyl-1,5-naphthyridine"),
    ("Cc1ccnc2cccnc12", "4-methyl-1,5-naphthyridine"),
    # 1,7-naphthyridine
    ("Cc1cnc2cnccc2c1", "3-methyl-1,7-naphthyridine"),
    ("Cc1ccnc2cnccc12", "4-methyl-1,7-naphthyridine"),
    ("Cc1cncc2ncccc12", "5-methyl-1,7-naphthyridine"),
    ("Cc1cc2cccnc2cn1", "6-methyl-1,7-naphthyridine"),
    ("Cc1nccc2cccnc12", "8-methyl-1,7-naphthyridine"),
    ("Cc1ccc2ccncc2n1", "2-methyl-1,7-naphthyridine"),
    # 1,8-naphthyridine
    ("Cc1cnc2ncccc2c1", "3-methyl-1,8-naphthyridine"),
    ("Cc1ccc2cccnc2n1", "2-methyl-1,8-naphthyridine"),
    ("Cc1ccnc2ncccc12", "4-methyl-1,8-naphthyridine"),
    # 2,6-naphthyridine
    ("Cc1cncc2ccncc12", "4-methyl-2,6-naphthyridine"),
    ("Cc1cc2cnccc2cn1", "3-methyl-2,6-naphthyridine"),
    ("Cc1nccc2cnccc12", "1-methyl-2,6-naphthyridine"),
    # acridine (N at 10)
    ("Cc1ccc2nc3ccccc3cc2c1", "3-methylacridine"),
    ("Cc1ccc2cc3ccccc3nc2c1", "2-methylacridine"),
    ("Cc1cccc2cc3ccccc3nc12", "1-methylacridine"),
    ("Cc1cccc2nc3ccccc3cc12", "4-methylacridine"),
    ("Cc1c2ccccc2nc2ccccc12",  "9-methylacridine"),
    # pteridine (N at 1,3,5,8)
    ("Cc1cnc2ncncc2n1", "6-methylpteridine"),
    ("Cc1cnc2cncnc2n1", "7-methylpteridine"),
    ("Cc1ncc2nccnc2n1", "2-methylpteridine"),
    ("Cc1ncnc2nccnc12", "4-methylpteridine"),
])
def test_phase623_methyl_diazanaphthalenes_acridine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
