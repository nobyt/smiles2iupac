"""Phase 894: silicon/germanium heterocycles (silolane, germolane, etc.) --
the Si/Ge mirror of Phase 888's boron-ring fix.

Found via a fresh probe sweep after Phase 893, prompted directly by having
just fixed the analogous boron case: C[Si]1(C)CCCC1 (a simple
silacyclopentane, 1,1-dimethylsilolane) was named "1,1-dimethylcyclopentane"
-- silicon silently read as an ordinary ring carbon, same bug class as the
cyclic boronic esters fixed in Phase 888. Same bug for germanium.

Root cause: _HETERO_SYMBOLS (ring_handler.py) didn't include "Si"/"Ge", so
those rings fell straight into the carbocycle path before ever reaching
heterocycle_handler.py's naming tables.

Fixed by adding "Si"/"Ge" to _HETERO_SYMBOLS, adding them to the
heteroatom-seniority _PRIORITY table (below B, matching real IUPAC
replacement-nomenclature order: ...Si > Ge > Sn > Pb > B), and adding
Hantzsch-Widman name entries for sizes 3-8 (silirane/siletane/silolane/
silinane/silepane/silocane and the germ- equivalents) -- all individually
verified via OPSIN parse-back.

DELIBERATELY NOT extended to multi-heteroatom rings combining Si/Ge with
O/N/S (e.g. an attempted "dioxasilolane" / mixed Si+O larger-ring a-
nomenclature): the underlying _canonical_sig tie-break (which picks
whichever ring-traversal direction sorts alphabetically smaller as a tuple
of element-symbol strings) happens to align with correct IUPAC numbering
for "B" (since "B" < "C" alphabetically) but NOT for "Si"/"Ge" (since
"Si"/"Ge" > "C" alphabetically) -- and a first attempt at the
multi-heteroatom case produced a name OPSIN could not parse
("5,5-dimethyl-1-oxa-5-silocane" -- "Unable to assign all locants"),
so it was reverted rather than landing an unverified, plausibly-wrong
naming path. The single-heteroatom Hantzsch-Widman path doesn't hit this
ambiguity (there's no second heteroatom position to get the direction
wrong relative to), so it stayed in scope.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C[Si]1(C)CCCC1",    "1,1-dimethylsilolane"),
    ("C[Ge]1(C)CCCC1",    "1,1-dimethylgermolane"),
    ("[SiH2]1CCCC1",      "silolane"),
    ("C[Si]1(C)CCCCC1",   "1,1-dimethylsilinane"),
    # regression: plain carbocycles unchanged
    ("C1CCCC1",  "cyclopentane"),
    ("C1CCCCC1", "cyclohexane"),
    # regression: existing O-only heterocycles unchanged
    ("C1COCCO1", "1,4-dioxane"),
    # regression: boron rings (Phase 888) unaffected by the Si/Ge additions
    ("CB1OCCO1", "2-methyl-1,3,2-dioxaborolane"),
])
def test_phase894_silicon_germanium_heterocycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase894_not_misnamed_as_carbocycle():
    result = smiles_to_iupac("C[Si]1(C)CCCC1")
    assert "cyclopentane" not in result
