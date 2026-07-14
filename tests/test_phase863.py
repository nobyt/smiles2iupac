"""Phase 863: halogen substituents on organoboranes.

R_nBHal_{3-n} (e.g. methyldifluoroborane) previously dropped the halogens
entirely -- CB(F)F was named "methylborane". Silicon already handled the
analogous case (C[Si](Cl)(Cl)Cl -> "trichloro(methyl)silane"); this mirrors
_name_organic_silane's halo+alkyl handling onto _name_organic_borane:

    CB(F)F   -> difluoro(methyl)borane
    CB(Cl)Cl -> dichloro(methyl)borane

Halogen prefixes and (parenthesised) alkyl names are collected and sorted
alphabetically, exactly as for silane. Amino substituents on boron
(aminoboranes) are a separate class not handled by the silane reference
either, so they are left for a follow-up.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # dihalo(alkyl)boranes
    ("CB(F)F",       "difluoro(methyl)borane"),
    ("CB(Cl)Cl",     "dichloro(methyl)borane"),
    ("CB(Br)Br",     "dibromo(methyl)borane"),
    ("CB(I)I",       "diiodo(methyl)borane"),
    ("CCB(Cl)Cl",    "dichloro(ethyl)borane"),
    # mixed halides
    ("CB(F)Cl",      "chlorofluoro(methyl)borane"),
    # halo + two alkyls
    ("CB(C)Cl",      "chloro(dimethyl)borane"),
    # regression: plain organoboranes and boron oxyacids unchanged
    ("CB",           "methylborane"),
    ("CCB",          "ethylborane"),
    ("CB(C)C",       "trimethylborane"),
    ("CB(O)O",       "methylboronic acid"),
    ("CB(C)O",       "dimethylborinic acid"),
    # regression: silicon reference unchanged
    ("C[Si](Cl)(Cl)Cl", "trichloro(methyl)silane"),
    ("C[Si](C)(C)C",    "tetramethylsilane"),
])
def test_phase863_haloboranes(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
