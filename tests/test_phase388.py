"""Phase 388: Mixed numeric/N-locant sort crash fix in heterocycle_handler (IUPAC 2013).

When a saturated N-heterocycle has substituents on both carbon ring atoms (int
locants) and on the ring nitrogen (string 'N' locant), comparing the two types
caused a TypeError crash in _subs_key.  Numeric locants sort before letter
locants per IUPAC 2013 P-14.5.2.

Bug: CC1CCN(C)CC1 raised TypeError: '<' not supported between instances of
     'str' and 'int'.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # C-substituent + N-substituent in same ring (caused crash before fix)
    ("CC1CCN(C)CC1",        "4,N-dimethylpiperidine"),
    ("CC1CCCCN1C",          "2,N-dimethylpiperidine"),
    ("CC1CCCN1C",           "2,N-dimethylpyrrolidine"),
    ("CCC1CCCN1C",          "2-ethyl-N-methylpyrrolidine"),
    # N-only substitution (regression: no crash, N-locant preserved)
    ("CN1CCCCC1",           "N-methylpiperidine"),
    ("CN1CCCC1",            "N-methylpyrrolidine"),
    ("CN1CCCCCCC1",         "N-methylazocane"),
    # C-only substitution (regression: numeric locants)
    ("CC1CCNCC1",           "4-methylpiperidine"),
    ("CC1CCCN1",            "2-methylpyrrolidine"),
    # Piperazine with mixed C/N substitution (lower C-locant preferred)
    ("CC1CN(C)CCN1",        "2,N'-dimethylpiperazine"),
])
def test_phase388_mixed_locant_sort(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
