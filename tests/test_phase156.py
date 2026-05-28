"""Phase 156: 環状化合物のロカント同順位アルファベット規則 (IUPAC 2013 P-14.5.2)

置換基ロカント集合が同一の場合、アルファベット順で早い置換基が
より低いロカントを得る。
例: 1-chloro-4-methylbenzene (chloro < methyl)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ハロゲン + メチル (ロカント集合が等しい場合のアルファベット順)
    ("Cc1ccc(Cl)cc1",    "1-chloro-4-methylbenzene"),
    ("Cc1ccc(F)cc1",     "1-fluoro-4-methylbenzene"),
    ("Cc1ccc(Br)cc1",    "1-bromo-4-methylbenzene"),
    # 回帰: 対称置換 (どちらの方向でも同一)
    ("Clc1ccc(Cl)cc1",   "1,4-dichlorobenzene"),
    ("Cc1ccc(C)cc1",     "1,4-dimethylbenzene"),
    ("Cc1cccc(C)c1",     "1,3-dimethylbenzene"),
    # アミノ基 (アルファベット順規則は不要な場合)
    ("Nc1ccc(Cl)cc1",    "4-chloroaniline"),
    ("Nc1ccc(C)cc1",     "4-methylaniline"),
    # 回帰
    ("CC(=O)O",          "acetic acid"),
    ("Nc1ccccc1",        "aniline"),
])
def test_phase156_ring_alpha_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
