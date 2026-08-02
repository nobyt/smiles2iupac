"""Phase 883: amino (-NH2) substituents on organoboranes/organosilanes.

CB(N)N was named "methylborane" -- dropping both -NH2 groups entirely,
colliding with the plain compound CB (also "methylborane"). Same info-loss
bug class as Phase 863 (which fixed the analogous halogen case,
CB(F)F -> difluoro(methyl)borane) but amine needed a DIFFERENT fix: it's
not a substituent prefix like halogens, it's the SUFFIX "-amine" (mirroring
how amine is always a suffix on carbon parents, e.g. methanediamine).

Verified via OPSIN parse-back (bundled opsin-cli-2.9.0-jar-with-dependencies.jar):
  "methylboranediamine"  -> CB(N)N       (matches, confirms the target)
  "methylsilanetriamine" -> C[Si](N)(N)N (matches, confirms the target)
  "aminomethylborane"    -> NCB          (WRONG molecule: amino on the
                                           methyl carbon, not on boron --
                                           rules out the naive prefix guess)

Implementation: new _simple_amino_neighbors (plain -NH2 only, N with no
other heavy neighbors) + _amino_suffix_word (elision: borane+amine ->
boranamine for n=1; no elision for di-/tri- since they start with a
consonant) + _name_borane_silane_with_halo_amino, a generalization of
Phase 863's halogen-mixing helper that also appends the amine suffix.
Mixed halogen+amino (e.g. CB(F)N) keeps halogen/alkyl as parenthesized
prefixes (same as Phase 863) with amine still as the suffix.

Scoped to simple, unsubstituted -NH2 only (matches this tool's established
carbon-required scope elsewhere -- a carbon-free aminoborane/aminosilane
like plain NB or N[SiH2]N is out of scope, same as the rest of the codebase).
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # the reported bug and its Si mirror
    ("CB(N)N",       "methylboranediamine"),
    ("C[Si](N)(N)N", "methylsilanetriamine"),
    # single amino (elision: borane+amine -> boranamine)
    ("CBN",          "methylboranamine"),
    # multiple alkyls + diamino
    ("C[Si](C)(N)N", "dimethylsilanediamine"),
    # mixed halogen + amino
    ("CB(F)N",       "fluoro(methyl)boranamine"),
    # regression: Phase 863/pre-existing halide-only cases unchanged
    ("CB(F)F",             "difluoro(methyl)borane"),
    ("C[Si](Cl)(Cl)Cl",    "trichloro(methyl)silane"),
    # regression: plain (no halogen/amino) organoborane/silane unchanged
    ("CB",           "methylborane"),
    ("C[Si](C)C",    "trimethylsilane"),
])
def test_phase883_amino_borane_silane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase883_not_confused_with_plain_methylborane():
    assert smiles_to_iupac("CB(N)N") != smiles_to_iupac("CB")
    assert smiles_to_iupac("C[Si](N)(N)N") != smiles_to_iupac("C[Si](C)C")
