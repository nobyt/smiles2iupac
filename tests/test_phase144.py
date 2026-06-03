"""Phase 144: 無機イオン・ヒドロキシルアミン 保留名 (IUPAC 2013 P-73)

hydroxylamine, ammonium, hydroxide, halide anions, common metal cations
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ヒドロキシルアミン
    ("NO",      "hydroxylamine"),
    # アンモニウム・水酸化物
    ("[NH4+]",  "ammonium"),
    ("[OH-]",   "hydroxide"),
    # ハロゲンアニオン
    ("[F-]",    "fluoride"),
    ("[Cl-]",   "chloride"),
    ("[Br-]",   "bromide"),
    ("[I-]",    "iodide"),
    # アルカリ金属カチオン
    ("[Na+]",   "sodium"),
    ("[K+]",    "potassium"),
    ("[Li+]",   "lithium"),
    # アルカリ土類・遷移金属
    ("[Ca+2]",  "calcium"),
    ("[Mg+2]",  "magnesium"),
    ("[Zn+2]",  "zinc"),
    ("[Fe+2]",  "iron(2+)"),
    ("[Fe+3]",  "iron(3+)"),
    ("[Al+3]",  "aluminium"),
    # 回帰: 通常化合物
    ("NN",      "hydrazine"),
    ("CC",      "ethane"),
    ("c1ccncc1", "pyridine"),
])
def test_phase144_ions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
