"""Phase 898: organosilicon/germanium/tin thiols (silanethiol etc.) --
found via the same fresh probe sweep as Phase 896/897, directly prompted
by checking whether the -OH gap fixed for these elements in Phase 886/894
also existed for -SH.

C[Si](C)(C)S (trimethylsilanethiol) was named "trimethylsilane" -- the
whole -SH group silently dropped, same as if it were the plain
tetraorgano compound. Same bug for germanium and tin.

Root cause: _detect_silicon_groups / _detect_germanium_tin_groups only
ever checked for O-H (silanol/germanol/stannanol, Phase 231/379/886) --
there was no S-H branch at all, so Si/Ge/Sn-SH fell through to the plain
organo-element name. Fixed by mirroring the existing O-H detection with an
S-H variant (new silanethiol_org/germanethiol_org/stannanethiol_org group
types) and a new _name_element_thiol_family namer mirroring the existing
_name_element_ol_family helper (Phase 895/896 lesson applied directly this
time: thiol-family suffixes are all consonant-starting so no elision logic
is needed, simpler than the ol-family case). Verified via OPSIN.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C[Si](C)(C)S",  "trimethylsilanethiol"),
    ("C[Ge](C)(C)S",  "trimethylgermanethiol"),
    ("C[Sn](C)(C)S",  "trimethylstannanethiol"),
    ("C[Si](C)(S)S",  "dimethylsilanedithiol"),
    # regression: the -OH family (Phase 231/379/886) unchanged
    ("C[Si](C)(C)O", "trimethylsilanol"),
    ("C[Ge](C)(C)O", "trimethylgermanol"),
    # regression: plain organo-element (no O/S) unchanged
    ("C[Si](C)(C)C", "tetramethylsilane"),
    ("C[Ge](C)(C)C", "tetramethylgermane"),
    # regression: halide handling (Phase 863/886) unaffected
    ("C[Si](C)(C)Cl", "chloro(trimethyl)silane"),
])
def test_phase898_element_thiols(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase898_not_confused_with_plain_silane():
    result = smiles_to_iupac("C[Si](C)(C)S")
    assert result != smiles_to_iupac("C[Si](C)(C)C")
