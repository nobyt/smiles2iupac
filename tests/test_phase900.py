"""Phase 900: cyanate/thiocyanate anions colliding with the neutral acids.

Investigated at the user's request after Phase 899 -- flagged as a
possible bug during that session's probe sweep but not yet verified.

[S-]C#N (the thiocyanate anion, SCN-) was named "thiocyanic acid" --
identical to the neutral acid HSCN (SC#N). Same bug for [O-]C#N
(cyanate anion) vs cyanic acid.

Root cause: _name_cyanate/_name_thiocyanate already distinguished the
ester form (R-O-C#N / R-S-C#N -> "cyanato{alkane}"/"thiocyanato{alkane}")
from the "free acid" form, but the free-acid branch returned "cyanic
acid"/"thiocyanic acid" unconditionally whenever there was no carbon
substituent on the O/S -- never checking whether that O/S carried an H
(true neutral acid) or a negative formal charge (the anion). Fixed by
checking formal_charge == -1 and returning "cyanate"/"thiocyanate" in that
case. Verified via OPSIN.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("[S-]C#N", "thiocyanate"),
    ("SC#N",    "thiocyanic acid"),
    ("[O-]C#N", "cyanate"),
    ("OC#N",    "cyanic acid"),
    # regression: ester forms unchanged
    ("CSC#N", "thiocyanatomethane"),
    ("COC#N", "cyanatomethane"),
])
def test_phase900_cyanate_thiocyanate_anions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase900_thiocyanate_not_confused_with_acid():
    assert smiles_to_iupac("[S-]C#N") != smiles_to_iupac("SC#N")


def test_phase900_cyanate_not_confused_with_acid():
    assert smiles_to_iupac("[O-]C#N") != smiles_to_iupac("OC#N")
