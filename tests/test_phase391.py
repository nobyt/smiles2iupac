"""Phase 391: Thiolactam naming – exocyclic =S on N-heterocycle ring C (IUPAC 2013 P-65.1.2.4).

Before this fix, S=C1CCCCN1 was named 'piperidine-2-thiol' (treating =S as
-SH substituent) instead of 'piperidine-2-thione' (exocyclic C=S suffix).

The '-thione' suffix starts with 't' (consonant), so the terminal 'e' of the
ring base name is NOT elided (contrast: '-one' elides the 'e').
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # piperidine-based thiolactam
    ("S=C1CCCCN1",        "piperidine-2-thione"),
    # pyrrolidine-based thiolactam (no 'e' elision: pyrrolidine not pyrrolidin)
    ("S=C1CCCN1",         "pyrrolidine-2-thione"),
    # N-substituted thiolactam
    ("CN1CCCC1=S",        "N-methylpyrrolidine-2-thione"),
    ("CN1CCCCC1=S",       "N-methylpiperidine-2-thione"),
    # regression: lactam (-one) still works (elides 'e' because -one starts with 'o')
    ("O=C1CCCCN1",        "piperidin-2-one"),
    ("O=C1CCCN1",         "pyrrolidin-2-one"),
    # regression: all-carbon ring thioketone unchanged
    ("S=C1CCCC1",         "cyclopentanethione"),
    ("S=C1CCCCC1",        "cyclohexanethione"),
])
def test_phase391_thiolactam(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
