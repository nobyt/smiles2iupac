"""Phase 126: ピペラジン・追加ヘテロ環保留名 (IUPAC 2013 P-31.1.3)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # piperazine (1,4-diazacyclohexane)
    ("C1CNCCN1", "piperazine"),
    ("N1CCNCC1", "piperazine"),
    ("C1NCCNC1", "piperazine"),
    # N-substituted piperazine
    ("CN1CCNCC1", "1-methylpiperazine"),
    ("CCN1CCNCC1", "1-ethylpiperazine"),
    # 回帰: piperidine unchanged
    ("C1CCCCN1", "piperidine"),
    # 回帰: morpholine unchanged
    ("C1CNCCO1", "morpholine"),
    # 回帰: pyrazine (aromatic) unchanged
    ("c1cnccn1", "pyrazine"),
])
def test_phase126_piperazine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
