"""Phase 155: 芳香族エステル命名 — ヘテロ芳香環・置換ベンゼン (IUPAC 2013)

ベンゼン環エステル: 置換基を正しくロカント付きで命名。
ヘテロ芳香環エステル: pyridine/furan/thiophene/pyrrole に -carboxylate suffix。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ヘテロ芳香環エステル
    ("COC(=O)c1ccncc1",    "methyl pyridine-4-carboxylate"),
    ("COC(=O)c1cccnc1",    "methyl pyridine-3-carboxylate"),
    ("COC(=O)c1ccccn1",    "methyl pyridine-2-carboxylate"),
    ("COC(=O)c1ccco1",     "methyl furan-2-carboxylate"),
    ("CCOC(=O)c1cccs1",    "ethyl thiophene-2-carboxylate"),
    ("COC(=O)c1ccc[nH]1",  "methyl 1H-pyrrole-2-carboxylate"),
    # 置換ベンゼン環エステル
    ("COC(=O)c1ccc(cc1)N", "methyl 4-aminobenzoate"),
    ("COC(=O)c1ccc(Cl)cc1","methyl 4-chlorobenzoate"),
    # 回帰: 純ベンゼンエステル (Phase 47)
    ("COC(=O)c1ccccc1",    "methyl benzoate"),
    ("CCOC(=O)c1ccccc1",   "ethyl benzoate"),
    # 回帰: 非芳香族エステル
    ("CCOC(=O)CC",         "ethyl propanoate"),
    ("CC(=O)O",            "acetic acid"),
])
def test_phase155_aromatic_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
