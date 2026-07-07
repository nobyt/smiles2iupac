"""Phase 846: drop mid-string -nH- indicated-H when N at that locant is substituted.

IUPAC 2013 P-14.7.1: when the atom bearing the indicated hydrogen is substituted,
the indicated hydrogen is omitted.  Extends Phases 844/845 to names where the nH
marker appears mid-string as -nH- (e.g. tetrahydro-1H-indole, dihydro-1H-benzimidazole).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 4,5,6,7-tetrahydro-1H-indole: N1-H → keep -1H-
    ("c1cc2c([nH]1)CCCC2",   "4,5,6,7-tetrahydro-1H-indole"),
    # N1-methyl → drop -1H-
    ("CN1C=CC=2CCCCC12",     "1-methyl-4,5,6,7-tetrahydroindole"),
    # C2-methyl: N1 still has H → keep -1H-
    ("CC=1NC=2CCCCC2C1",     "2-methyl-4,5,6,7-tetrahydro-1H-indole"),
    # 4,5,6,7-tetrahydro-1H-benzimidazole: N1-methyl → drop -1H-
    ("CN1C=NC2=C1CCCC2",     "1-methyl-4,5,6,7-tetrahydrobenzimidazole"),
    # C2-methyl → keep -1H-
    ("CC1=NC2=C(N1)CCCC2",   "2-methyl-4,5,6,7-tetrahydro-1H-benzimidazole"),
    # 4,5,6,7-tetrahydro-1H-indazole: N1-methyl → drop -1H-
    ("CN1N=CC=2CCCCC12",     "1-methyl-4,5,6,7-tetrahydroindazole"),
    # C3-methyl → keep -1H-
    ("CC1=NNC=2CCCCC12",     "3-methyl-4,5,6,7-tetrahydro-1H-indazole"),
    # 4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole: N1-methyl → drop -1H-
    ("CN1N=NC2=C1CCCC2",     "1-methyl-4,5,6,7-tetrahydrobenzo[d][1,2,3]triazole"),
    # 2,3-dihydro-1H-benzimidazole: N1-methyl → drop -1H-
    ("CN1CNC2=C1C=CC=C2",    "1-methyl-2,3-dihydrobenzimidazole"),
    # C2-methyl → keep -1H-
    ("CC1NC2=C(N1)C=CC=C2",  "2-methyl-2,3-dihydro-1H-benzimidazole"),
    # 2,3-dihydro-1H-indazole: N1-methyl → drop -1H-
    ("CN1NCC2=CC=CC=C12",    "1-methyl-2,3-dihydroindazole"),
    # C1-methyl in 1H-isochromene (C position!) → keep -1H-
    ("CC1OCCc2ccccc21",      "1-methyl-3,4-dihydro-1H-isochromene"),
    # Separator retained when -nH- precedes a locant digit
    # 2,3-dihydro-1H-1,5-benzodiazepine: N1-methyl → drop -1H-, keep separator before 1,5
    ("CN1CCC=Nc2ccccc21",    "1-methyl-2,3-dihydro-1,5-benzodiazepine"),
    # unsubstituted → keep -1H-
    ("C1=Nc2ccccc2NCC1",     "2,3-dihydro-1H-1,5-benzodiazepine"),
    # 3,4-dihydro-1H-1-benzazepin-2(5H)-one: N1-methyl → drop -1H-
    ("CN1C(=O)CCCc2ccccc21", "1-methyl-3,4-dihydro-1-benzazepin-2(5H)-one"),
    # 2,3-dihydro-1H-1,4-benzodiazepin-5(4H)-one: N1-methyl → drop -1H-
    ("CN1CCNC(=O)c2ccccc21", "1-methyl-2,3-dihydro-1,4-benzodiazepin-5(4H)-one"),
])
def test_phase846(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
