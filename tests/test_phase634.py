"""Phase 634: 4,5,6,7-tetrahydrobenzo[d]-fused 5-membered heterocycles —
benzo[d]oxazole, benzo[d]thiazole, benzo[d]isoxazole, benzo[d]isothiazole,
1H-benzo[d][1,2,3]triazole, and 2H-benzo[d][1,2,3]triazole (IUPAC 2013).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("O1C=NC2=C1CCCC2",   "4,5,6,7-tetrahydrobenzo[d]oxazole"),
    ("S1C=NC2=C1CCCC2",   "4,5,6,7-tetrahydrobenzo[d]thiazole"),
    ("O1N=CC2=C1CCCC2",   "4,5,6,7-tetrahydrobenzo[d]isoxazole"),
    ("S1N=CC2=C1CCCC2",   "4,5,6,7-tetrahydrobenzo[d]isothiazole"),
    ("N1N=NC2=C1CCCC2",   "4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole"),
    ("N=1NN=C2C1CCCC2",   "4,5,6,7-tetrahydro-2H-benzo[d][1,2,3]triazole"),
    # 4,5,6,7-tetrahydrobenzo[d]oxazole: O(1)-C(2)-N(3,None)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    ("CC=1OC2=C(N1)CCCC2",  "2-methyl-4,5,6,7-tetrahydrobenzo[d]oxazole"),
    ("CC1CCCC2=C1N=CO2",    "4-methyl-4,5,6,7-tetrahydrobenzo[d]oxazole"),
    ("CC1CCC2=C(N=CO2)C1",  "5-methyl-4,5,6,7-tetrahydrobenzo[d]oxazole"),
    ("CC1CC2=C(N=CO2)CC1",  "6-methyl-4,5,6,7-tetrahydrobenzo[d]oxazole"),
    ("CC1CCCC=2N=COC21",    "7-methyl-4,5,6,7-tetrahydrobenzo[d]oxazole"),
    # 4,5,6,7-tetrahydrobenzo[d]thiazole: S(1)-C(2)-N(3,None)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    ("CC=1SC2=C(N1)CCCC2",  "2-methyl-4,5,6,7-tetrahydrobenzo[d]thiazole"),
    ("CC1CCCC2=C1N=CS2",    "4-methyl-4,5,6,7-tetrahydrobenzo[d]thiazole"),
    ("CC1CCC2=C(N=CS2)C1",  "5-methyl-4,5,6,7-tetrahydrobenzo[d]thiazole"),
    ("CC1CC2=C(N=CS2)CC1",  "6-methyl-4,5,6,7-tetrahydrobenzo[d]thiazole"),
    ("CC1CCCC=2N=CSC21",    "7-methyl-4,5,6,7-tetrahydrobenzo[d]thiazole"),
    # 4,5,6,7-tetrahydrobenzo[d]isoxazole: O(1)-N(2,None)-C(3)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    ("CC1=NOC2=C1CCCC2",    "3-methyl-4,5,6,7-tetrahydrobenzo[d]isoxazole"),
    ("CC1CCCC2=C1C=NO2",    "4-methyl-4,5,6,7-tetrahydrobenzo[d]isoxazole"),
    ("CC1CCC2=C(C=NO2)C1",  "5-methyl-4,5,6,7-tetrahydrobenzo[d]isoxazole"),
    ("CC1CC2=C(C=NO2)CC1",  "6-methyl-4,5,6,7-tetrahydrobenzo[d]isoxazole"),
    ("CC1CCCC=2C=NOC21",    "7-methyl-4,5,6,7-tetrahydrobenzo[d]isoxazole"),
    # 4,5,6,7-tetrahydrobenzo[d]isothiazole: S(1)-N(2,None)-C(3)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    ("CC1=NSC2=C1CCCC2",    "3-methyl-4,5,6,7-tetrahydrobenzo[d]isothiazole"),
    ("CC1CCCC2=C1C=NS2",    "4-methyl-4,5,6,7-tetrahydrobenzo[d]isothiazole"),
    ("CC1CCC2=C(C=NS2)C1",  "5-methyl-4,5,6,7-tetrahydrobenzo[d]isothiazole"),
    ("CC1CC2=C(C=NS2)CC1",  "6-methyl-4,5,6,7-tetrahydrobenzo[d]isothiazole"),
    ("CC1CCCC=2C=NSC21",    "7-methyl-4,5,6,7-tetrahydrobenzo[d]isothiazole"),
    # 1H-benzo[d][1,2,3]triazole: N(1,H)-N(2)-N(3)-C(3a,junc)-C(7a,junc)
    # C7a=adj to N1; C7(adj C7a), C6, C5, C4(adj C3a)
    ("CN1N=NC2=C1CCCC2",    "1-methyl-4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole"),
    ("CC1CCCC=2NN=NC21",    "4-methyl-4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole"),
    ("CC1CC2=C(NN=N2)CC1",  "5-methyl-4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole"),
    ("CC1CCC2=C(NN=N2)C1",  "6-methyl-4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole"),
    ("CC1CCCC2=C1NN=N2",    "7-methyl-4,5,6,7-tetrahydro-1H-benzo[d][1,2,3]triazole"),
    # 2H-benzo[d][1,2,3]triazole: C2-symmetric (4≡7→4, 5≡6→5), N(2,H) at atom5
    ("CC1CCCc2n[nH]nc21",   "4-methyl-4,5,6,7-tetrahydro-2H-benzo[d][1,2,3]triazole"),
    ("CC1CCCc2n[nH]nc21",   "4-methyl-4,5,6,7-tetrahydro-2H-benzo[d][1,2,3]triazole"),  # 7-methyl→4-methyl
    ("CC1CCc2n[nH]nc2C1",   "5-methyl-4,5,6,7-tetrahydro-2H-benzo[d][1,2,3]triazole"),
    ("CC1CCc2n[nH]nc2C1",   "5-methyl-4,5,6,7-tetrahydro-2H-benzo[d][1,2,3]triazole"),  # 6-methyl→5-methyl
])
def test_phase634_tetrahydrobenzod_5membered(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
