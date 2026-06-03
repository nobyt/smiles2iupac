"""Phase 153: テトラゾール・1,3,5-トリアジナン・hexahydropyrimidine (IUPAC 2013)

1H-tetrazole: 5員 aromatic 4-N 環の保留名。
1,3,5-triazinane-2,4,6-trione: 環内 3 つの exo C=O を持つ6員三窒素環 (-trione は子音始まり→ e 保持)。
hexahydropyrimidine-2,4-dione: 6員 N1,N3 飽和環のジオン (-dione は子音始まり→ e 保持)。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # テトラゾール (5員 aromatic, 4 N atoms)
    ("c1nnn[nH]1",           "1H-tetrazole"),
    # 1,3,5-triazinane-2,4,6-trione (isocyanuric acid 非芳香型; 'd' consonant → no elision)
    ("O=C1NC(=O)NC(=O)N1",   "1,3,5-triazinane-2,4,6-trione"),
    # hexahydropyrimidine-2,4-dione (dihydrouracil 型)
    ("O=C1NC(=O)NCC1",        "hexahydropyrimidine-2,4-dione"),
    # 既存回帰: 類似 dione/trione
    ("O=C1CCC(=O)N1",         "pyrrolidine-2,5-dione"),
    ("O=C1NC(=O)CN1",         "imidazolidine-2,4-dione"),
    ("O=C1CCC(=O)O1",         "oxolane-2,5-dione"),
    # 回帰
    ("CC(=O)O",               "acetic acid"),
    ("c1ncc2[nH]cnc2n1",      "9H-purine"),
])
def test_phase153_tetrazole_triazinane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
