"""Phase 146: カルボキシレートアニオン・アンモニウムイオン・無機オキソ酸 (IUPAC 2013)

carboxylate anions (formate, acetate, propanoate...),
ammonium ions (methylazanium, tetramethylazanium...),
inorganic oxoacids (sulfuric acid, nitric acid, nitrous acid)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # カルボキシレートアニオン
    ("O=C[O-]",           "formate"),
    ("CC(=O)[O-]",         "acetate"),
    ("CCC(=O)[O-]",        "propanoate"),
    ("CCCC(=O)[O-]",       "butanoate"),
    ("CCCCC(=O)[O-]",      "pentanoate"),
    # アンモニウムイオン
    ("C[NH3+]",            "methylazanium"),
    ("C[NH2+]C",           "dimethylazanium"),
    ("C[NH+](C)C",         "trimethylazanium"),
    ("C[N+](C)(C)C",       "tetramethylazanium"),
    ("CC[N+](CC)(CC)CC",   "tetraethylazanium"),
    # 無機オキソ酸
    ("O=S(=O)(O)O",        "sulfuric acid"),
    ("O=S(O)O",            "sulfurous acid"),
    ("O=[N+]([O-])O",      "nitric acid"),
    ("O=NO",               "nitrous acid"),
    ("O=C(O)O",            "carbonic acid"),
    # 回帰: 中性アミン・カルボン酸
    ("CCN",                "ethanamine"),
    ("CC(=O)O",            "acetic acid"),
    ("C[N+](C)(C)[O-]",    "N,N-dimethylmethanamine N-oxide"),
    ("CCCC[N+](=O)[O-]",   "1-nitrobutane"),
])
def test_phase146_ions_and_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
