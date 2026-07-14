"""Phase 864: thio-substituted phosphonic / phosphinic acids.

Sulfur-substituted analogues of phosphonic acid (R-P(=O)(OH)2) and phosphinic
acid (R2-P(=O)(OH)) were mis-named or garbled: CP(=S)(O)O was wrongly called
"methylphosphonous acid" (dropping the =S and demoting P(V) to P(III)).

Naming (IUPAC): the number of sulfur atoms on P (=S plus -SH) selects
thio/dithio/trithio; NO positional prefix is used, because e.g.
R-P(=S)(OH)2 and R-P(=O)(OH)(SH) are rapidly interconverting tautomers that
cannot be isolated or distinguished, so IUPAC assigns them a single name:

  1 S  -> {alkyl(s)}phospho{no|ni}thioic acid
  2 S  -> ...phospho{no|ni}dithioic acid
  3 S  -> ...phosphonotrithioic acid   (c=1 only)
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phosphonothioic (c=1): tautomer pairs share one name (no O,O-/O,S- tag)
    ("CP(=S)(O)O",   "methylphosphonothioic acid"),
    ("CP(=O)(S)O",   "methylphosphonothioic acid"),
    ("CCP(=S)(O)O",  "ethylphosphonothioic acid"),
    # phosphonodithioic (2 S)
    ("CP(=O)(S)S",   "methylphosphonodithioic acid"),
    ("CP(=S)(S)O",   "methylphosphonodithioic acid"),
    # phosphonotrithioic (3 S)
    ("CP(=S)(S)S",   "methylphosphonotrithioic acid"),
    # phosphinothioic (c=2)
    ("CP(C)(=S)O",   "dimethylphosphinothioic acid"),
    ("CP(C)(=O)S",   "dimethylphosphinothioic acid"),
    # phosphinodithioic (2 S, c=2)
    ("CP(C)(=S)S",   "dimethylphosphinodithioic acid"),
    # regression: oxy acids unchanged
    ("CP(=O)(O)O",   "methylphosphonic acid"),
    ("CP(C)(=O)O",   "dimethylphosphinic acid"),
    ("CP(O)O",       "methylphosphonous acid"),
    ("CP(C)O",       "dimethylphosphinous acid"),
])
def test_phase864_thiophosphorus_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
