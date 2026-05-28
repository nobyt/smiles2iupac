"""Phase 121: チオアルデヒド・チオケトン命名 (IUPAC 2013 P-65.3)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # thioaldehydes (thial suffix)
    ("C=S", "methanethial"),
    ("CC=S", "ethanethial"),
    ("CCC=S", "propanethial"),
    ("CCCC=S", "butanethial"),
    # thioketones (thione suffix)
    ("CC(=S)C", "propane-2-thione"),
    ("CCC(=S)C", "butane-2-thione"),
    ("CCC(=S)CC", "pentane-3-thione"),
    # 回帰: isothiocyanate (N=C=S) は影響なし
    ("CCN=C=S", "ethyl isothiocyanate"),
    # 回帰: thioamide (C(=S)-N) は影響なし
    ("CC(=S)N", "ethanethioamide"),
    # 回帰: thioether (C-S-C)
    ("CSC", "dimethyl sulfide"),
])
def test_phase121_thioketone_thioaldehyde(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
