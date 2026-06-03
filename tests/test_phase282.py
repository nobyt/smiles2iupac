"""Phase 282: sulfide (thioether) substitutive PIN (IUPAC 2013 P-63.6.1.1).

Sulfides R-S-R' are named as "(Rsulfanyl)parent" where the longer chain/ring is
the parent and the shorter -S-R group is the (Rsulfanyl) prefix.  The functional-
class retained name "R R' sulfide" is acceptable but is not the IUPAC 2013 PIN.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── symmetric sulfides ────────────────────────────────────────────────
    ("CSC",                    "dimethyl sulfide"),
    ("CCSCC",                  "diethyl sulfide"),

    # ── asymmetric: longer chain is parent ───────────────────────────────
    ("CSCC",                   "ethyl methyl sulfide"),
    ("CSCCC",                  "methyl propyl sulfide"),

    # ── aryl-containing: ring is parent ─────────────────────────────────
    ("CSc1ccccc1",             "methyl phenyl sulfide"),

    # ── regressions: thioester and disulfide unchanged ────────────────────
    ("CSC(=O)C",               "S-methyl ethanethioate"),
    ("CSSC",                   "dimethyl disulfide"),
    ("CCS",                    "ethanethiol"),
])
def test_phase282_sulfide_pin(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
