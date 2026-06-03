"""Phase 375: Sulfenyl halide and sulfenate ester naming (IUPAC 2013 P-65.4).

Sulfenyl halides R-S-X (X = Cl, Br, F, I): named '{stem}anesulfenyl {halide}'.
Sulfenate esters R-S-O-R': named '{alkyl} {stem}anesulfenate'.

Previously both patterns fell through all S-detection branches (which required
either S=O or S-OH) and were named as thioethers or ethers.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Sulfenyl chlorides
    ("CSCl",           "methanesulfenyl chloride"),
    ("CCSCl",          "ethanesulfenyl chloride"),
    ("CCCSCl",         "propane-1-sulfenyl chloride"),
    # Sulfenyl bromide/fluoride
    ("CSBr",           "methanesulfenyl bromide"),
    ("CSF",            "methanesulfenyl fluoride"),
    # Sulfenate esters
    ("CSOC",           "methyl methanesulfenate"),
    ("CSOCC",          "ethyl methanesulfenate"),
    ("CCSOC",          "methyl ethanesulfenate"),
    ("CCSOCC",         "ethyl ethanesulfenate"),
    # Regressions: sulfenic acid unchanged
    ("CSO",            "methanesulfenic acid"),
    ("CCSO",           "ethanesulfenic acid"),
    # Regressions: sulfide unchanged
    ("CSC",            "dimethyl sulfide"),
    ("CSCC",           "ethyl methyl sulfide"),
])
def test_phase375_sulfenyl_halide_sulfenate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
