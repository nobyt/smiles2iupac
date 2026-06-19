"""Phase 633: 4,5,6,7-tetrahydrobenzo-fused 5-membered aromatic heterocycles —
4,5,6,7-tetrahydrobenzofuran, 4,5,6,7-tetrahydrobenzothiophene,
4,5,6,7-tetrahydro-1H-indole, 4,5,6,7-tetrahydro-1H-benzimidazole,
and 4,5,6,7-tetrahydro-1H-indazole, with methyl substituents at all positions
including N-methyl variants (IUPAC 2013).
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # parent compounds
    ("O1C=CC2=C1CCCC2",    "4,5,6,7-tetrahydrobenzofuran"),
    ("S1C=CC2=C1CCCC2",    "4,5,6,7-tetrahydrobenzothiophene"),
    ("N1C=CC=2CCCCC12",    "4,5,6,7-tetrahydro-1H-indole"),
    ("N1C=NC2=C1CCCC2",    "4,5,6,7-tetrahydro-1H-benzimidazole"),
    ("N1N=CC=2CCCCC12",    "4,5,6,7-tetrahydro-1H-indazole"),
    # 4,5,6,7-tetrahydrobenzofuran: O(1)-C(2)-C(3)-C(3a,junc)-C(4)-C(5)-C(6)-C(7)-C(7a,junc)
    ("CC=1OC2=C(C1)CCCC2",  "2-methyl-4,5,6,7-tetrahydrobenzofuran"),
    ("CC1=COC2=C1CCCC2",    "3-methyl-4,5,6,7-tetrahydrobenzofuran"),
    ("CC1CCCC2=C1C=CO2",    "4-methyl-4,5,6,7-tetrahydrobenzofuran"),
    ("CC1CCC2=C(C=CO2)C1",  "5-methyl-4,5,6,7-tetrahydrobenzofuran"),
    ("CC1CC2=C(C=CO2)CC1",  "6-methyl-4,5,6,7-tetrahydrobenzofuran"),
    ("CC1CCCC=2C=COC21",    "7-methyl-4,5,6,7-tetrahydrobenzofuran"),
    # 4,5,6,7-tetrahydrobenzothiophene: S(1)-C(2)-C(3)-C(3a)-C(4)-C(5)-C(6)-C(7)-C(7a)
    ("CC=1SC2=C(C1)CCCC2",  "2-methyl-4,5,6,7-tetrahydrobenzothiophene"),
    ("CC1=CSC2=C1CCCC2",    "3-methyl-4,5,6,7-tetrahydrobenzothiophene"),
    ("CC1CCCC2=C1C=CS2",    "4-methyl-4,5,6,7-tetrahydrobenzothiophene"),
    ("CC1CCC2=C(C=CS2)C1",  "5-methyl-4,5,6,7-tetrahydrobenzothiophene"),
    ("CC1CC2=C(C=CS2)CC1",  "6-methyl-4,5,6,7-tetrahydrobenzothiophene"),
    ("CC1CCCC=2C=CSC21",    "7-methyl-4,5,6,7-tetrahydrobenzothiophene"),
    # 4,5,6,7-tetrahydro-1H-indole: N(1H)-C(2)-C(3)-C(3a)-C(4)-C(5)-C(6)-C(7)-C(7a)
    ("CN1C=CC=2CCCCC12",    "1-methyl-4,5,6,7-tetrahydro-1H-indole"),
    ("CC=1NC=2CCCCC2C1",    "2-methyl-4,5,6,7-tetrahydro-1H-indole"),
    ("CC1=CNC=2CCCCC12",    "3-methyl-4,5,6,7-tetrahydro-1H-indole"),
    ("CC1C=2C=CNC2CCC1",    "4-methyl-4,5,6,7-tetrahydro-1H-indole"),
    ("CC1CC=2C=CNC2CC1",    "5-methyl-4,5,6,7-tetrahydro-1H-indole"),
    ("CC1CCC=2C=CNC2C1",    "6-methyl-4,5,6,7-tetrahydro-1H-indole"),
    ("CC1CCCC=2C=CNC12",    "7-methyl-4,5,6,7-tetrahydro-1H-indole"),
    # 4,5,6,7-tetrahydro-1H-benzimidazole: N(1H)-C(2)=N(3)-C(3a)-C(4)-C(5)-C(6)-C(7)-C(7a)
    # C2-symmetric (4≡7, 5≡6), positions 2, 4(≡7), 5(≡6) unique + N1, N3
    ("CN1C=NC2=C1CCCC2",    "1-methyl-4,5,6,7-tetrahydro-1H-benzimidazole"),
    ("CC1=NC2=C(N1)CCCC2",  "2-methyl-4,5,6,7-tetrahydro-1H-benzimidazole"),
    ("CC1CCCC=2NC=NC21",    "4-methyl-4,5,6,7-tetrahydro-1H-benzimidazole"),
    ("CC1CC2=C(NC=N2)CC1",  "5-methyl-4,5,6,7-tetrahydro-1H-benzimidazole"),
    # 4,5,6,7-tetrahydro-1H-indazole: N(1H)-N(2)=C(3)-C(3a)-C(4)-C(5)-C(6)-C(7)-C(7a)
    ("CN1N=CC=2CCCCC12",    "1-methyl-4,5,6,7-tetrahydro-1H-indazole"),
    ("CC1=NNC=2CCCCC12",    "3-methyl-4,5,6,7-tetrahydro-1H-indazole"),
    ("CC1C=2C=NNC2CCC1",    "4-methyl-4,5,6,7-tetrahydro-1H-indazole"),
    ("CC1CC=2C=NNC2CC1",    "5-methyl-4,5,6,7-tetrahydro-1H-indazole"),
    ("CC1CCC=2C=NNC2C1",    "6-methyl-4,5,6,7-tetrahydro-1H-indazole"),
    ("CC1CCCC=2C=NNC12",    "7-methyl-4,5,6,7-tetrahydro-1H-indazole"),
])
def test_phase633_tetrahydrobenzo_5membered(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
