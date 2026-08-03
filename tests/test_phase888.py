"""Phase 888: cyclic boronic esters (1,3,2-dioxaborolane / -dioxaborinane).

Found via the fresh probe sweep (same session as 884-887). This was a
structurally WRONG name, not just info loss: CB1OCCO1 (2-methyl-1,3,2-
dioxaborolane, a simple cyclic boronic ester) was named "methylcyclopentane"
-- boron was silently read as an ordinary ring carbon. The pinacol boronate
ester family (4,4,5,5-tetramethyl-2-aryl-1,3,2-dioxaborolanes) is one of
the most common reagent classes in modern synthetic chemistry (Suzuki
coupling), so this was a significant real-world gap, not an exotic edge
case.

Root cause (two layers, both needed fixing):
1. _HETERO_SYMBOLS (ring_handler.py) -- the set that decides whether a ring
   is even treated as heterocyclic at all -- didn't include "B", so any
   boron-containing ring fell straight into the carbocycle naming path,
   which has no way to represent boron and just miscounted it as carbon.
2. heterocycle_handler.py's ring-composition pattern table (_RETAINED_NAMES,
   keyed by canonical heteroatom-signature tuples) had no boron entries.
   Also added "B" to the _PRIORITY heteroatom-seniority table (lowest
   priority, matching real IUPAC replacement-nomenclature seniority: O
   ranks above B, so O correctly becomes the ring's locant-1 starting
   point) -- though the existing .get(sig, 99) fallback already gave B the
   right (lowest) priority by accident, adding it explicitly documents the
   real seniority position for any future ring combining B with N/P/Se/Te.

Both new ring patterns (5-membered O,B,O,C,C and 6-membered O,B,O,C,C,C)
and every generated name verified via OPSIN parse-back before landing.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # bare parent rings
    ("B1OCCO1",   "1,3,2-dioxaborolane"),
    ("B1OCCCO1",  "1,3,2-dioxaborinane"),
    # substituted at boron (the common case -- boronic ester of an alcohol)
    ("CB1OCCO1",  "2-methyl-1,3,2-dioxaborolane"),
    ("CB1OCCCO1", "2-methyl-1,3,2-dioxaborinane"),
    # the real-world case: pinacol boronate ester of phenylboronic acid
    ("CC1(C)OB(OC1(C)C)c1ccccc1", "4,4,5,5-tetramethyl-2-phenyl-1,3,2-dioxaborolane"),
    # regression: plain carbocycles of the same ring size unchanged
    ("C1CCCC1",  "cyclopentane"),
    ("C1CCCCC1", "cyclohexane"),
    # regression: existing O,O-heterocycles (no boron) unchanged
    ("C1COCO1",  "1,3-dioxolane"),
    # regression: non-cyclic boron compounds unchanged
    ("CB(O)O",   "methylboronic acid"),
    ("CB(O)C",   "dimethylborinic acid"),
])
def test_phase888_dioxaborolane_dioxaborinane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase888_not_misnamed_as_carbocycle():
    result = smiles_to_iupac("CB1OCCO1")
    assert "cyclopentane" not in result
