"""Phase 199: N-ハロアミン命名 (IUPAC 2013 P-62.3.1)

  ClNC  → N-chloromethanamine
  ClNCC → N-chloroethanamine
  BrNC  → N-bromomethanamine

N に C + ハロゲンが付く場合は N-{halide} 接頭辞で命名 (置換命名法)。
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 一級 N-ハロアミン (N に C 1 個 + ハロゲン 1 個)
    ("ClNC",    "N-chloromethanamine"),
    ("BrNC",    "N-bromomethanamine"),
    ("FNC",     "N-fluoromethanamine"),
    ("ClNCC",   "N-chloroethanamine"),
    # 二級 N-ハロアミン (N に C 2 個 + ハロゲン 1 個)
    ("ClN(C)C", "N-chloro-N-methylmethanamine"),
    # 回帰: 通常アミンは変わらない
    ("NC",      "methanamine"),
    ("NCC",     "ethanamine"),
    ("CNCC",    "N-methylethanamine"),
])
def test_phase199_n_haloamine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
