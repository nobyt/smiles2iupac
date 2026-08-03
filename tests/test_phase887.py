"""Phase 887: sulfonylureas -- fix a severe collision where the whole
sulfonyl group was dropped and the compound collapsed to bare "urea".

Found via a fresh probe sweep. CS(=O)(=O)NC(=O)N (a sulfonylurea -- the core
scaffold of an entire real drug class, e.g. glibenclamide-type
antidiabetics) was named "urea", identical to the actual parent compound
NC(=O)N. All 3 heavy atoms of the methanesulfonyl group were silently
dropped.

Root cause (two bugs, same underlying gap): the urea-detection helpers only
ever checked whether a urea nitrogen carried a CARBON substituent (to defer
to the N-alkylurea namer) -- they never considered other heteroatom
substituents like a sulfonyl group.

1. _is_urea's "both N are plain -NH2" guard excluded N-alkyl/N-aryl (checked
   for a carbon neighbor) but not N-sulfonyl -- so the sulfonylurea's
   substituted nitrogen looked "unsubstituted" to it, and the whole molecule
   was confidently reported as bare "urea".
2. Once that was fixed (so _is_urea correctly says no), the compound fell
   through to _name_substituted_urea_if_match, whose get_subs() also only
   ever looked for carbon substituents on each urea nitrogen -- so the
   sulfonyl substituent was invisible there too, and the function returned
   None (treating the sulfonyl-bearing N as if it had no substituent),
   letting the molecule fall through to an unrelated, nonsensical
   "1-aminoformamide" fallback further down the pipeline.

Fixed both: _is_urea's guard now requires every non-carbonyl neighbor of
each flanking N to be literally H (not just "not carbon"); get_subs() in
_name_substituted_urea_if_match now also recognizes a R-SO2- substituent on
the nitrogen and names it "{alkane}sulfonyl", feeding it through the same
N-/N'- prefix-building logic already used for alkyl substituents. All new
name forms verified via OPSIN parse-back.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # the reported bug and its ethyl mirror
    ("CS(=O)(=O)NC(=O)N",  "N-methanesulfonylurea"),
    ("CCS(=O)(=O)NC(=O)N", "N-ethanesulfonylurea"),
    # mixed: sulfonyl on one N, alkyl on the other
    ("CS(=O)(=O)NC(=O)NC", "N-methanesulfonyl-N'-methylurea"),
    # regression: plain urea and N-alkylureas (Phase 49/83) unchanged
    ("NC(=O)N",     "urea"),
    ("CNC(=O)N",    "N-methylurea"),
    ("CNC(=O)NC",   "N,N'-dimethylurea"),
    # regression: thiourea (separate detector, checked earlier) unaffected
    ("NC(=S)N",     "thiourea"),
])
def test_phase887_sulfonylurea(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase887_not_confused_with_plain_urea():
    assert smiles_to_iupac("CS(=O)(=O)NC(=O)N") != smiles_to_iupac("NC(=O)N")
