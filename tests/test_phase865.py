"""Phase 865: thio-substituted sulfonic / sulfinic acids.

Sulfur analogue of the Phase 864 thiophosphonic acids. Sulfur-substituted
sulfonic acid (R-S(=O)2-OH) and sulfinic acid (R-S(=O)-OH) were mis-named or
garbled: CS(=S)(=O)O was wrongly "methanesulfinic acid" (dropping the =S and
demoting sulfonic->sulfinic), and CS(=O)(=O)S garbled to "sulfonylmethane".

Naming: the number of S on the central sulfur (=S plus -SH) selects
thio/dithio/trithio; the number of doubly-bonded chalcogens selects the level
(2 -> sulfono, 1 -> sulfino). NO positional prefix, by the same tautomer
principle as the phosphorus acids -- e.g. R-S(=O)(=S)-OH and R-S(=O)2-SH
interconvert too fast to isolate, so IUPAC gives them one name:

  CS(=S)(=O)O / CS(=O)(=O)S -> methanesulfonothioic acid
  CS(=S)(=S)O / CS(=S)(=O)S -> methanesulfonodithioic acid
  CS(=S)(=S)S               -> methanesulfonotrithioic acid
  CS(=S)O / CS(=O)S         -> methanesulfinothioic acid
  CS(=S)S                   -> methanesulfinodithioic acid

The carboxylic thioic S-acid/O-acid convention (Phase 149) is deliberately
NOT merged -- it is a separate retained IUPAC convention.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # sulfonothioic (2 double chalcogens): tautomer pairs share one name
    ("CS(=S)(=O)O",         "methanesulfonothioic acid"),
    ("CS(=O)(=O)S",         "methanesulfonothioic acid"),
    ("CCS(=S)(=O)O",        "ethanesulfonothioic acid"),
    # sulfonodithioic (2 S)
    ("CS(=S)(=S)O",         "methanesulfonodithioic acid"),
    ("CS(=S)(=O)S",         "methanesulfonodithioic acid"),
    # sulfonotrithioic (3 S)
    ("CS(=S)(=S)S",         "methanesulfonotrithioic acid"),
    # sulfinothioic (1 double chalcogen)
    ("CS(=S)O",             "methanesulfinothioic acid"),
    ("CS(=O)S",             "methanesulfinothioic acid"),
    # sulfinodithioic (2 S)
    ("CS(=S)S",             "methanesulfinodithioic acid"),
    # aromatic
    ("c1ccccc1S(=S)(=O)O",  "benzenesulfonothioic acid"),
    # regression: oxy acids and other S groups unchanged
    ("CS(=O)(=O)O",         "methanesulfonic acid"),
    ("CCCS(=O)(=O)O",       "propane-1-sulfonic acid"),
    ("CS(=O)O",             "methanesulfinic acid"),
    ("CSO",                 "methanesulfenic acid"),
    ("CS(=O)(=O)Cl",        "methanesulfonyl chloride"),
    ("CS(=O)(=O)N",         "methanesulfonamide"),
    ("CS(=O)(=O)OC",        "methyl methanesulfonate"),
])
def test_phase865_thiosulfur_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
