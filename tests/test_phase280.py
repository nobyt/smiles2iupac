"""Phase 280: methylidene substituents on acyclic chains (IUPAC 2013 P-31.1.8.4).

An exo C=C on a chain is named with the "-ylidene" suffix when the functional-group
principal chain is forced to bypass the double bond (e.g. by carboxyl or aldehyde
groups selecting a longer saturated path).  The preferred name for =CH2 is
"methylidene" (not "methylene") on acyclic parents.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── carboxyl-anchored chains ──────────────────────────────────────────
    # itaconic acid: COOH-CH2-C(=CH2)-COOH
    ("OC(=O)CC(=C)C(=O)O",  "2-methylidenebutanedioic acid"),
    ("OC(=O)C(=C)CC(=O)O",  "2-methylidenebutanedioic acid"),
    # 2-methylidenepentanoic acid
    ("OC(=O)C(=C)CCC",      "2-methylidenepentanoic acid"),

    # ── aldehyde-anchored chains ──────────────────────────────────────────
    # CHO-C(=CH2)-CCC  →  2-methylidenepentanal
    ("O=CC(=C)CCC",         "2-methylidenepentanal"),

    # ── regressions: plain alkene, not alkylidene ─────────────────────────
    # chain_finder pulls C=C into the main chain → ordinary alkene name
    ("CC(=C)CC",            "2-methylbut-1-ene"),
    ("C=C(C)CC",            "2-methylbut-1-ene"),
])
def test_phase280_chain_alkylidene(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
