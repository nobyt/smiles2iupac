"""Phase 249: halosilane naming (IUPAC 2013 P-68.3).

Organohalosilanes: halogen substituents cited alphabetically with alkyl
substituents as prefixes to silane.
  chloro(methyl)silane, dichloro(ethyl)silane, etc.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # monohalogenosilanes (1 halogen + 1 alkyl)
    ("Cl[SiH2]C",          "chloro(methyl)silane"),
    ("Cl[SiH2]CC",         "chloro(ethyl)silane"),
    ("F[SiH2]C",           "fluoro(methyl)silane"),
    ("Br[SiH2]C",          "bromo(methyl)silane"),
    # dihalo + alkyl
    ("[SiH](Cl)(Cl)C",     "dichloro(methyl)silane"),
    ("[SiH](Cl)(Cl)CC",    "dichloro(ethyl)silane"),
    # trihalo + alkyl
    ("[Si](Cl)(Cl)(Cl)C",  "trichloro(methyl)silane"),
    # mono halo + multiple alkyls
    ("Cl[Si](C)(C)C",      "chloro(trimethyl)silane"),
    # regression: pure alkylsilane unchanged
    ("[SiH3]C",            "methylsilane"),
    ("[SiH2](C)C",         "dimethylsilane"),
    # regression: halo-only silane (via retained names) unchanged
    ("[SiH2](Cl)Cl",       "dichlorosilane"),
    ("[SiH3]Cl",           "chlorosilane"),
])
def test_phase249_halosilane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
