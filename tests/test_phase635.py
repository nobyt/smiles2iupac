"""Phase 635: 2,3-dihydrobenzo[d]-fused 5-membered heterocycles —
2,3-dihydrobenzo[d]oxazole, thiazole, isoxazole, isothiazole,
2,3-dihydro-1H-benzimidazole, and 2,3-dihydro-1H-indazole (IUPAC 2013).
The benzene ring is aromatic; the 5-membered ring has 2 sp3 atoms.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("O1CNC2=C1C=CC=C2",  "2,3-dihydrobenzo[d]oxazole"),
    ("S1CNC2=C1C=CC=C2",  "2,3-dihydrobenzo[d]thiazole"),
    ("O1NCC2=C1C=CC=C2",  "2,3-dihydrobenzo[d]isoxazole"),
    ("S1NCC2=C1C=CC=C2",  "2,3-dihydrobenzo[d]isothiazole"),
    ("N1CNC2=C1C=CC=C2",  "2,3-dihydro-1H-benzimidazole"),
    ("N1NCC2=CC=CC=C12",  "2,3-dihydro-1H-indazole"),
    # 2,3-dihydrobenzo[d]oxazole: O(1,None)-C(2,sp3)-N(3,NH)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    ("CC1OC2=C(N1)C=CC=C2",  "2-methyl-2,3-dihydrobenzo[d]oxazole"),
    ("CN1COC2=C1C=CC=C2",    "3-methyl-2,3-dihydrobenzo[d]oxazole"),
    ("CC1=CC=CC2=C1NCO2",    "4-methyl-2,3-dihydrobenzo[d]oxazole"),
    ("CC=1C=CC2=C(NCO2)C1",  "5-methyl-2,3-dihydrobenzo[d]oxazole"),
    ("CC1=CC2=C(NCO2)C=C1",  "6-methyl-2,3-dihydrobenzo[d]oxazole"),
    ("CC1=CC=CC=2NCOC21",    "7-methyl-2,3-dihydrobenzo[d]oxazole"),
    # 2,3-dihydrobenzo[d]thiazole: S(1,None)-C(2,sp3)-N(3,NH)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    ("CC1SC2=C(N1)C=CC=C2",  "2-methyl-2,3-dihydrobenzo[d]thiazole"),
    ("CN1CSC2=C1C=CC=C2",    "3-methyl-2,3-dihydrobenzo[d]thiazole"),
    ("CC1=CC=CC2=C1NCS2",    "4-methyl-2,3-dihydrobenzo[d]thiazole"),
    ("CC=1C=CC2=C(NCS2)C1",  "5-methyl-2,3-dihydrobenzo[d]thiazole"),
    ("CC1=CC2=C(NCS2)C=C1",  "6-methyl-2,3-dihydrobenzo[d]thiazole"),
    ("CC1=CC=CC=2NCSC21",    "7-methyl-2,3-dihydrobenzo[d]thiazole"),
    # 2,3-dihydrobenzo[d]isoxazole: O(1,None)-N(2,NH)-C(3,sp3)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    ("CN1OC2=C(C1)C=CC=C2",  "2-methyl-2,3-dihydrobenzo[d]isoxazole"),
    ("CC1NOC2=C1C=CC=C2",    "3-methyl-2,3-dihydrobenzo[d]isoxazole"),
    ("CC1=CC=CC2=C1CNO2",    "4-methyl-2,3-dihydrobenzo[d]isoxazole"),
    ("CC=1C=CC2=C(CNO2)C1",  "5-methyl-2,3-dihydrobenzo[d]isoxazole"),
    ("CC1=CC2=C(CNO2)C=C1",  "6-methyl-2,3-dihydrobenzo[d]isoxazole"),
    ("CC1=CC=CC=2CNOC21",    "7-methyl-2,3-dihydrobenzo[d]isoxazole"),
    # 2,3-dihydrobenzo[d]isothiazole: S(1,None)-N(2,NH)-C(3,sp3)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    ("CN1SC2=C(C1)C=CC=C2",  "2-methyl-2,3-dihydrobenzo[d]isothiazole"),
    ("CC1NSC2=C1C=CC=C2",    "3-methyl-2,3-dihydrobenzo[d]isothiazole"),
    ("CC1=CC=CC2=C1CNS2",    "4-methyl-2,3-dihydrobenzo[d]isothiazole"),
    ("CC=1C=CC2=C(CNS2)C1",  "5-methyl-2,3-dihydrobenzo[d]isothiazole"),
    ("CC1=CC2=C(CNS2)C=C1",  "6-methyl-2,3-dihydrobenzo[d]isothiazole"),
    ("CC1=CC=CC=2CNSC21",    "7-methyl-2,3-dihydrobenzo[d]isothiazole"),
    # 2,3-dihydro-1H-benzimidazole: C2-symmetric (N1≡N3→1, C4≡C7→4, C5≡C6→5)
    ("CN1CNC2=C1C=CC=C2",    "1-methyl-2,3-dihydro-1H-benzimidazole"),
    ("CC1NC2=C(N1)C=CC=C2",  "2-methyl-2,3-dihydro-1H-benzimidazole"),
    ("CC1=CC=CC=2NCNC21",    "4-methyl-2,3-dihydro-1H-benzimidazole"),
    ("CC1=CC2=C(NCN2)C=C1",  "5-methyl-2,3-dihydro-1H-benzimidazole"),
    ("CN1CNC2=C1C=CC=C2",    "1-methyl-2,3-dihydro-1H-benzimidazole"),  # 3-methyl→1-methyl (C2-sym)
    # 2,3-dihydro-1H-indazole: N(1)-N(2)-C(3,sp3)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    ("CN1NCC2=CC=CC=C12",    "1-methyl-2,3-dihydro-1H-indazole"),
    ("CN1NC2=CC=CC=C2C1",    "2-methyl-2,3-dihydro-1H-indazole"),
    ("CC1NNC2=CC=CC=C12",    "3-methyl-2,3-dihydro-1H-indazole"),
    ("CC1=C2CNNC2=CC=C1",    "4-methyl-2,3-dihydro-1H-indazole"),
    ("CC=1C=C2CNNC2=CC1",    "5-methyl-2,3-dihydro-1H-indazole"),
    ("CC1=CC=C2CNNC2=C1",    "6-methyl-2,3-dihydro-1H-indazole"),
    ("CC=1C=CC=C2CNNC12",    "7-methyl-2,3-dihydro-1H-indazole"),
])
def test_phase635_dihydrobenzod_5membered(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
