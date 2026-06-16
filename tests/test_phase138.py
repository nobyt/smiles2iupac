"""Phase 138: 多環芳香族炭化水素 (PAH) 保留名 (IUPAC 2013 P-31.1.2)

azulene, 1H-indene, acenaphthylene, acenaphthene, fluoranthene,
chrysene, triphenylene, pyrene, perylene, benzo[a]pyrene, coronene
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # azulene (5+7 fused carbocyclic, C10H8)
    ("C1=CC2=CC=CC=CC2=C1", "azulene"),
    # 1H-indene (5+6 partially unsaturated, C9H8)
    ("C1=CC2=CC=CC=C2C1", "1H-indene"),
    # acenaphthylene (3-ring with 5-membered, C12H8)
    ("C1=CC2=CC=CC3=CC=CC1=C23", "acenaphthylene"),
    # acenaphthene (partially saturated acenaphthylene, C12H10)
    ("C1CC2=CC=CC3=CC=CC1=C23", "acenaphthene"),
    # fluoranthene (C16H10, fused with 5-membered ring)
    ("c1ccc2c(c1)-c1cccc3cccc-2c13", "fluoranthene"),
    # chrysene (4 fused 6-membered rings, C18H12)
    ("c1ccc2c(c1)ccc1c3ccccc3ccc21", "chrysene"),
    # triphenylene (4 fused 6-membered rings, symmetric, C18H12)
    ("c1ccc2c(c1)c1ccccc1c1ccccc21", "triphenylene"),
    # pyrene (4 fused rings, C16H10)
    ("c1cc2ccc3cccc4ccc(c1)c2c34", "pyrene"),
    # perylene (5 fused rings, C20H12)
    ("c1cc2cccc3c4cccc5cccc(c(c1)c23)c54", "perylene"),
    # benzo[a]pyrene (C20H12)
    ("c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34", "benzo[a]pyrene"),
    # coronene (C24H12)
    ("c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61", "coronene"),
    # 回帰: previously named PAH unchanged
    ("c1ccc2ccccc2c1",         "naphthalene"),
    ("c1ccc2cc3ccccc3cc2c1",   "anthracene"),
    ("c1ccc2c(c1)ccc1ccccc12", "phenanthrene"),
    # 回帰: fused heteroaromatics unchanged
    ("c1ccc2ncccc2c1",   "quinoline"),
    ("c1ccc2[nH]ccc2c1", "1H-indole"),
    # 回帰: Phase 133/134 partially saturated unchanged
    ("C1CCc2ccccc21", "indane"),
    ("C1COc2ccccc2C1", "chromane"),
])
def test_phase138_pah_retained(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
