"""Phase 630: benzo-fused 7-membered partially saturated rings —
6,7,8,9-tetrahydro-5H-benzo[7]annulene, benzazepines, benzoxepines,
benzothiepine, benzodiazepines, and their lactam/imine forms,
with methyl substituents at unique positions (IUPAC 2013).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("c1ccc2c(c1)CCCCC2", "6,7,8,9-tetrahydro-5H-benzo[7]annulene"),
    ("c1ccc2c(c1)CCCCN2", "2,3,4,5-tetrahydro-1-benzazepine"),
    ("c1ccc2c(c1)CCCNC2", "2,3,4,5-tetrahydro-1H-2-benzazepine"),
    ("c1ccc2c(c1)CCNCC2", "2,3,4,5-tetrahydro-1H-3-benzazepine"),
    ("c1ccc2c(c1)CCCCO2", "2,3,4,5-tetrahydro-1-benzoxepine"),
    ("c1ccc2c(c1)CCCOC2", "2,3,4,5-tetrahydro-1H-2-benzoxepine"),
    ("c1ccc2c(c1)CCOCC2", "2,3,4,5-tetrahydro-1H-3-benzoxepine"),
    ("c1ccc2c(c1)CCCCS2", "2,3,4,5-tetrahydro-1-benzothiepine"),
    ("c1ccc2c(c1)NCCCN2", "2,3,4,5-tetrahydro-1H-1,5-benzodiazepine"),
    ("C1=Nc2ccccc2NCC1",  "2,3-dihydro-1H-1,5-benzodiazepine"),
    ("O=C1CCCc2ccccc2N1", "3,4-dihydro-1H-1-benzazepin-2(5H)-one"),
    ("O=C1NCCNc2ccccc21", "2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one"),
    # 6,7,8,9-tetrahydro-5H-benzo[7]annulene: unique positions 1, 2, 5, 6, 7
    ("Cc1cccc2c1CCCCC2",  "1-methyl-6,7,8,9-tetrahydro-5H-benzo[7]annulene"),
    ("Cc1ccc2c(c1)CCCCC2", "2-methyl-6,7,8,9-tetrahydro-5H-benzo[7]annulene"),
    ("CC1CCCCc2ccccc21",  "5-methyl-6,7,8,9-tetrahydro-5H-benzo[7]annulene"),
    ("CC1CCCc2ccccc2C1",  "6-methyl-6,7,8,9-tetrahydro-5H-benzo[7]annulene"),
    ("CC1CCc2ccccc2CC1",  "7-methyl-6,7,8,9-tetrahydro-5H-benzo[7]annulene"),
    # 2,3,4,5-tetrahydro-1-benzazepine: N at 1; unique 1–9
    ("CN1CCCCc2ccccc21",  "1-methyl-2,3,4,5-tetrahydro-1-benzazepine"),
    ("CC1CCCc2ccccc2N1",  "2-methyl-2,3,4,5-tetrahydro-1-benzazepine"),
    ("CC1CCc2ccccc2NC1",  "3-methyl-2,3,4,5-tetrahydro-1-benzazepine"),
    ("CC1CCNc2ccccc2C1",  "4-methyl-2,3,4,5-tetrahydro-1-benzazepine"),
    ("CC1CCCNc2ccccc21",  "5-methyl-2,3,4,5-tetrahydro-1-benzazepine"),
    ("Cc1cccc2c1CCCCN2",  "6-methyl-2,3,4,5-tetrahydro-1-benzazepine"),
    ("Cc1ccc2c(c1)CCCCN2", "7-methyl-2,3,4,5-tetrahydro-1-benzazepine"),
    ("Cc1ccc2c(c1)NCCCC2", "8-methyl-2,3,4,5-tetrahydro-1-benzazepine"),
    ("Cc1cccc2c1NCCCC2",  "9-methyl-2,3,4,5-tetrahydro-1-benzazepine"),
    # 2,3,4,5-tetrahydro-1H-2-benzazepine: N at 2; unique 1, 3–9
    ("CC1NCCCc2ccccc21",  "1-methyl-2,3,4,5-tetrahydro-1H-2-benzazepine"),
    ("CC1CCc2ccccc2CN1",  "3-methyl-2,3,4,5-tetrahydro-1H-2-benzazepine"),
    ("CC1CNCc2ccccc2C1",  "4-methyl-2,3,4,5-tetrahydro-1H-2-benzazepine"),
    ("CC1CCNCc2ccccc21",  "5-methyl-2,3,4,5-tetrahydro-1H-2-benzazepine"),
    ("Cc1cccc2c1CCCNC2",  "6-methyl-2,3,4,5-tetrahydro-1H-2-benzazepine"),
    ("Cc1ccc2c(c1)CCCNC2", "7-methyl-2,3,4,5-tetrahydro-1H-2-benzazepine"),
    ("Cc1ccc2c(c1)CNCCC2", "8-methyl-2,3,4,5-tetrahydro-1H-2-benzazepine"),
    ("Cc1cccc2c1CNCCC2",  "9-methyl-2,3,4,5-tetrahydro-1H-2-benzazepine"),
    # 2,3,4,5-tetrahydro-1H-3-benzazepine: N at 3; C2-symmetric (1≡5, 2≡4, 6≡9, 7≡8)
    ("CC1CNCCc2ccccc21",  "1-methyl-2,3,4,5-tetrahydro-1H-3-benzazepine"),
    ("CC1Cc2ccccc2CCN1",  "2-methyl-2,3,4,5-tetrahydro-1H-3-benzazepine"),
    ("Cc1cccc2c1CCNCC2",  "6-methyl-2,3,4,5-tetrahydro-1H-3-benzazepine"),
    ("Cc1ccc2c(c1)CCNCC2", "7-methyl-2,3,4,5-tetrahydro-1H-3-benzazepine"),
    # 2,3,4,5-tetrahydro-1-benzoxepine: O at 1; unique 2–9
    ("CC1CCCc2ccccc2O1",  "2-methyl-2,3,4,5-tetrahydro-1-benzoxepine"),
    ("CC1CCc2ccccc2OC1",  "3-methyl-2,3,4,5-tetrahydro-1-benzoxepine"),
    ("CC1CCOc2ccccc2C1",  "4-methyl-2,3,4,5-tetrahydro-1-benzoxepine"),
    ("CC1CCCOc2ccccc21",  "5-methyl-2,3,4,5-tetrahydro-1-benzoxepine"),
    ("Cc1cccc2c1CCCCO2",  "6-methyl-2,3,4,5-tetrahydro-1-benzoxepine"),
    ("Cc1ccc2c(c1)CCCCO2", "7-methyl-2,3,4,5-tetrahydro-1-benzoxepine"),
    ("Cc1ccc2c(c1)OCCCC2", "8-methyl-2,3,4,5-tetrahydro-1-benzoxepine"),
    ("Cc1cccc2c1OCCCC2",  "9-methyl-2,3,4,5-tetrahydro-1-benzoxepine"),
    # 2,3,4,5-tetrahydro-1H-2-benzoxepine: O at 2; unique 1, 3–9
    ("CC1OCCCc2ccccc21",  "1-methyl-2,3,4,5-tetrahydro-1H-2-benzoxepine"),
    ("CC1CCc2ccccc2CO1",  "3-methyl-2,3,4,5-tetrahydro-1H-2-benzoxepine"),
    ("CC1COCc2ccccc2C1",  "4-methyl-2,3,4,5-tetrahydro-1H-2-benzoxepine"),
    ("CC1CCOCc2ccccc21",  "5-methyl-2,3,4,5-tetrahydro-1H-2-benzoxepine"),
    ("Cc1cccc2c1CCCOC2",  "6-methyl-2,3,4,5-tetrahydro-1H-2-benzoxepine"),
    ("Cc1ccc2c(c1)CCCOC2", "7-methyl-2,3,4,5-tetrahydro-1H-2-benzoxepine"),
    ("Cc1ccc2c(c1)COCCC2", "8-methyl-2,3,4,5-tetrahydro-1H-2-benzoxepine"),
    ("Cc1cccc2c1COCCC2",  "9-methyl-2,3,4,5-tetrahydro-1H-2-benzoxepine"),
    # 2,3,4,5-tetrahydro-1H-3-benzoxepine: O at 3; C2-symmetric (1≡5, 2≡4, 6≡9, 7≡8)
    ("CC1COCCc2ccccc21",  "1-methyl-2,3,4,5-tetrahydro-1H-3-benzoxepine"),
    ("CC1Cc2ccccc2CCO1",  "2-methyl-2,3,4,5-tetrahydro-1H-3-benzoxepine"),
    ("Cc1cccc2c1CCOCC2",  "6-methyl-2,3,4,5-tetrahydro-1H-3-benzoxepine"),
    ("Cc1ccc2c(c1)CCOCC2", "7-methyl-2,3,4,5-tetrahydro-1H-3-benzoxepine"),
    # 2,3,4,5-tetrahydro-1-benzothiepine: S at 1; unique 2–9
    ("CC1CCCc2ccccc2S1",  "2-methyl-2,3,4,5-tetrahydro-1-benzothiepine"),
    ("CC1CCc2ccccc2SC1",  "3-methyl-2,3,4,5-tetrahydro-1-benzothiepine"),
    ("CC1CCSc2ccccc2C1",  "4-methyl-2,3,4,5-tetrahydro-1-benzothiepine"),
    ("CC1CCCSc2ccccc21",  "5-methyl-2,3,4,5-tetrahydro-1-benzothiepine"),
    ("Cc1cccc2c1CCCCS2",  "6-methyl-2,3,4,5-tetrahydro-1-benzothiepine"),
    ("Cc1ccc2c(c1)CCCCS2", "7-methyl-2,3,4,5-tetrahydro-1-benzothiepine"),
    ("Cc1ccc2c(c1)SCCCC2", "8-methyl-2,3,4,5-tetrahydro-1-benzothiepine"),
    ("Cc1cccc2c1SCCCC2",  "9-methyl-2,3,4,5-tetrahydro-1-benzothiepine"),
    # 2,3,4,5-tetrahydro-1H-1,5-benzodiazepine: C2-symmetric (N1≡N5, C2≡C4, C6≡C9, C7≡C8)
    ("CC1CCNc2ccccc2N1",  "2-methyl-2,3,4,5-tetrahydro-1H-1,5-benzodiazepine"),
    ("CC1CNc2ccccc2NC1",  "3-methyl-2,3,4,5-tetrahydro-1H-1,5-benzodiazepine"),
    ("Cc1cccc2c1NCCCN2",  "6-methyl-2,3,4,5-tetrahydro-1H-1,5-benzodiazepine"),
    ("Cc1ccc2c(c1)NCCCN2", "7-methyl-2,3,4,5-tetrahydro-1H-1,5-benzodiazepine"),
    # 2,3-dihydro-1H-1,5-benzodiazepine: N1-C2-C3-C4=N5; unique 1–4, 6–9
    ("CN1CCC=Nc2ccccc21",  "1-methyl-2,3-dihydro-1H-1,5-benzodiazepine"),
    ("CC1CC=Nc2ccccc2N1",  "2-methyl-2,3-dihydro-1H-1,5-benzodiazepine"),
    ("CC1C=Nc2ccccc2NC1",  "3-methyl-2,3-dihydro-1H-1,5-benzodiazepine"),
    ("CC1=Nc2ccccc2NCC1",  "4-methyl-2,3-dihydro-1H-1,5-benzodiazepine"),
    ("Cc1cccc2c1N=CCCN2",  "6-methyl-2,3-dihydro-1H-1,5-benzodiazepine"),
    ("Cc1ccc2c(c1)N=CCCN2", "7-methyl-2,3-dihydro-1H-1,5-benzodiazepine"),
    ("Cc1ccc2c(c1)NCCC=N2", "8-methyl-2,3-dihydro-1H-1,5-benzodiazepine"),
    ("Cc1cccc2c1NCCC=N2",  "9-methyl-2,3-dihydro-1H-1,5-benzodiazepine"),
    # 3,4-dihydro-1H-1-benzazepin-2(5H)-one: N1; C2=O; unique 1, 3–9
    ("CN1C(=O)CCCc2ccccc21", "1-methyl-3,4-dihydro-1H-1-benzazepin-2(5H)-one"),
    ("CC1CCc2ccccc2NC1=O",   "3-methyl-3,4-dihydro-1H-1-benzazepin-2(5H)-one"),
    ("CC1CC(=O)Nc2ccccc2C1", "4-methyl-3,4-dihydro-1H-1-benzazepin-2(5H)-one"),
    ("CC1CCC(=O)Nc2ccccc21", "5-methyl-3,4-dihydro-1H-1-benzazepin-2(5H)-one"),
    ("Cc1cccc2c1CCCC(=O)N2", "6-methyl-3,4-dihydro-1H-1-benzazepin-2(5H)-one"),
    ("Cc1ccc2c(c1)CCCC(=O)N2", "7-methyl-3,4-dihydro-1H-1-benzazepin-2(5H)-one"),
    ("Cc1ccc2c(c1)NC(=O)CCC2", "8-methyl-3,4-dihydro-1H-1-benzazepin-2(5H)-one"),
    ("Cc1cccc2c1NC(=O)CCC2",   "9-methyl-3,4-dihydro-1H-1-benzazepin-2(5H)-one"),
    # 2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one: N1-C2-C3-N4-C5=O; unique 1–4, 6–9
    ("CN1CCNC(=O)c2ccccc21", "1-methyl-2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one"),
    ("CC1CNC(=O)c2ccccc2N1", "2-methyl-2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one"),
    ("CC1CNc2ccccc2C(=O)N1", "3-methyl-2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one"),
    ("CN1CCNc2ccccc2C1=O",   "4-methyl-2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one"),
    ("Cc1cccc2c1C(=O)NCCN2", "6-methyl-2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one"),
    ("Cc1ccc2c(c1)C(=O)NCCN2", "7-methyl-2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one"),
    ("Cc1ccc2c(c1)NCCNC2=O", "8-methyl-2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one"),
    ("Cc1cccc2c1NCCNC2=O",   "9-methyl-2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one"),
])
def test_phase630_benzo_fused_7membered(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
