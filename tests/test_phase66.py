"""Phase 66: 環状ジオール / ジオンの立体化学記述子"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 環状ジオール + 立体化学
    ("O[C@@H]1CCCC[C@H]1O", "(1R,2R)-cyclohexane-1,2-diol"),
    ("O[C@H]1CCCC[C@@H]1O", "(1S,2S)-cyclohexane-1,2-diol"),
    ("O[C@@H]1CCCC[C@@H]1O", "(1R,2S)-cyclohexane-1,2-diol"),
    # 立体化学なし
    ("OC1CCCCC1O", "cyclohexane-1,2-diol"),
    # 環状ヒドロキシケトン + 立体化学
    ("O=C1CCCC[C@@H]1O", "(2S)-2-hydroxycyclohexan-1-one"),
    ("O=C1CCC[C@@H]1O", "(2S)-2-hydroxycyclopentan-1-one"),
    # シクロペンタンジオール
    ("O[C@@H]1CCC[C@H]1O", "(1R,2R)-cyclopentane-1,2-diol"),
    ("OC1CCCC1O", "cyclopentane-1,2-diol"),
])
def test_phase66_ring_stereo_diol(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
