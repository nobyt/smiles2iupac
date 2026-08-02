"""Phase 885: sulfonimidamide + fix a same-session regression in sulfilimine/
sulfoximine (Phase 884) that silently collapsed sulfondiimide into them.

Found via a fresh probe sweep after Phase 884 landed. Two bugs, same family:

1. NEW GROUP: sulfonimidamide C-S(=O)(=NR)-NR'2, the =NR-extended analog of
   sulfinamide. CS(=O)(=N)N was wrongly "methanesulfinamide" -- a DIFFERENT
   real compound (only 1 O, no imide N). Root cause: the existing
   sulfinamide branch only checked `len(n_neighbors) >= 1`, never the N's
   bond order, so it silently absorbed the =NH as if it were just another
   ordinary amide-N neighbor. Fixed by inserting a sulfonimidamide branch
   (exactly 1 double-bonded N + exactly 1 single-bonded N) before
   sulfinamide. Verified via OPSIN parse-back: "methanesulfonimidamide" ->
   CS(=N)(N)=O, matching. Amide N takes the "N-" locant, imide N takes
   "N'-" (also OPSIN-confirmed).

2. REGRESSION IN PHASE 884 ITSELF: the sulfoximine/sulfilimine detection
   added this session used `any(...)` to check for a double-bonded N,
   which also (wrongly) matches when there are TWO double-bonded N's --
   i.e. a sulfondiimide C-S(=NR)2-C. CS(=N)(=N)C was silently reported as
   "dimethyl sulfilimine" (a real but DIFFERENT compound with only one
   =NR), dropping the second imine entirely. OPSIN confirms the correct
   name for this structure needs a distinct systematic parent
   ("dimethyl-lambda6-sulfanediimine") that this codebase's existing
   "sulf-" functional-class naming style doesn't cover -- rather than
   half-implement lambda-convention naming (used nowhere else in this
   codebase), the detection guards were tightened to require EXACTLY one
   double-bonded N, so the two-imine case now honestly falls through
   (unmodeled) instead of confidently claiming the wrong compound.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # sulfonimidamide: base + both N-locant variants
    ("CS(=O)(=N)N",   "methanesulfonimidamide"),
    ("CS(=O)(=N)NC",  "N-methylmethanesulfonimidamide"),
    ("CS(=O)(=NC)N",  "N'-methylmethanesulfonimidamide"),
    # regression: plain sulfinamide/sulfonamide (no imide N) unchanged
    ("CS(=O)N",       "methanesulfinamide"),
    ("CS(=O)(=O)N",   "methanesulfonamide"),
    # regression: sulfilimine/sulfoximine (single =N, Phase 884) unchanged
    ("C[S](C)=N",       "dimethyl sulfilimine"),
    ("C[S](C)(=O)=N",   "dimethyl sulfoximine"),
])
def test_phase885_sulfonimidamide_and_regressions(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase885_sulfonimidamide_not_confused_with_sulfinamide():
    assert smiles_to_iupac("CS(=O)(=N)N") != smiles_to_iupac("CS(=O)N")


def test_phase885_sulfondiimide_not_misclaimed_as_sulfilimine():
    # a sulfondiimide (2 imide N on S) must not be silently reported as the
    # single-imine sulfilimine -- honest non-claim is safer than a wrong name
    result = smiles_to_iupac("CS(=N)(=N)C")
    assert result != "dimethyl sulfilimine"
