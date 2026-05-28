"""Phase 133: 部分飽和縮合環保留名 (IUPAC 2013 P-31.1.2, P-31.1.6)

indane, 1,2,3,4-tetrahydronaphthalene, chromane, indoline,
indan-1-one, indolin-2-one, and related partially saturated fused systems
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # indane (2,3-dihydro-1H-indene)
    ("C1CCc2ccccc21", "indane"),
    # 1,2,3,4-tetrahydronaphthalene (tetralin)
    ("C1CCCc2ccccc21", "1,2,3,4-tetrahydronaphthalene"),
    # chromane (3,4-dihydro-2H-chromene)
    ("C1COc2ccccc2C1", "chromane"),
    # indoline (2,3-dihydroindole)
    ("C1Cc2ccccc2N1", "indoline"),
    # indan-1-one (1-indanone)
    ("O=C1CCc2ccccc21", "indan-1-one"),
    # indolin-2-one (oxindole)
    ("O=C1Cc2ccccc2N1", "indolin-2-one"),
    # 1,2,3,4-tetrahydroquinoline
    ("C1CCNc2ccccc21", "1,2,3,4-tetrahydroquinoline"),
    # 回帰: aromatic fused systems still work
    ("c1ccc2ncccc2c1", "quinoline"),
    ("c1ccc2occc2c1", "benzofuran"),
    ("c1ccc2[nH]ccc2c1", "1H-indole"),
    # 回帰: simple rings still work
    ("C1CCCCC1", "cyclohexane"),
    ("c1ccccc1", "benzene"),
])
def test_phase133_partially_saturated_fused(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
