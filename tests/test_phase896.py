"""Phase 896: chalcogenolate anions (thiolate/selenolate/tellurolate) --
the S/Se/Te mirror of Phase 895's alkoxide fix.

Found via a fresh probe sweep prompted directly by the alkoxide fix:
CC[S-] (ethanethiolate) was named "sulfanylethane" -- the negative charge
silently dropped, reporting the anion as identical to the neutral thioether
substituent name. Same bug for Se and Te.

Root cause: same shape as Phase 895's alkoxide gap -- there was no
"thiolate"/"selenolate"/"tellurolate" functional-group detection at all,
so the deprotonated chalcogen fell through to a generic substituent-naming
path. Fixed by mirroring the existing thiol/selenol/tellurol detectors
(same elif chain, keyed on formal_charge == -1 and no H instead of an
H neighbor) and adding FunctionalGroupSpec entries with suffix
"thiolate"/"selenolate"/"tellurolate". Also added the required explicit
per-suffix branches in name_assembler.py's _build_name_body (per the
Phase 895 lesson: FunctionalGroupSpec.chain_template is dead code and a
missing branch silently garbles the name), mirroring the existing
thiol/selenol/tellurol branches' locant/elision rules exactly (all
consonant-starting, so "ane" is never elided). Verified via OPSIN.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CC[S-]",       "ethanethiolate"),
    ("CC[Se-]",      "ethaneselenolate"),
    ("CC[Te-]",      "ethanetellurolate"),
    ("c1ccccc1[S-]", "benzenethiolate"),
    ("CCC[S-]",      "propane-1-thiolate"),
    # regression: neutral thiol/selenol unchanged
    ("CCS",     "ethanethiol"),
    ("CC[SeH]", "ethaneselenol"),
    ("CCCS",    "propane-1-thiol"),
])
def test_phase896_chalcogenolate_anions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase896_thiolate_not_confused_with_thioether():
    result = smiles_to_iupac("CC[S-]")
    assert "sulfanylethane" not in result
