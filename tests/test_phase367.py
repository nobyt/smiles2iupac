"""Phase 367: N-locant prefix for non-aromatic ring nitrogen substituents.

Previously _collect_hetero_substituents returned numeric locants for all ring
atoms, so N-substituted saturated heterocycles were named with '1-' instead of
'N-' (e.g. '1-methylpiperidine' instead of 'N-methylpiperidine').

IUPAC 2013 P-31.1.3.4: when a non-aromatic ring nitrogen bears a substituent
not part of the ring name, use 'N' (or "N'" for the second N, etc.) as the
locant.  Aromatic N-heterocycles continue to use numeric locants.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Single-N saturated rings
    ("C1CCCCN1C",          "N-methylpiperidine"),
    ("C1CCCN1C",           "N-methylpyrrolidine"),
    ("C1CN(C)CCO1",        "N-methylmorpholine"),
    ("C1CCCCN1CC",         "N-ethylpiperidine"),
    # Lactams
    ("O=C1CCCCN1C",        "N-methylpiperidin-2-one"),
    ("O=C1CCCN1C",         "N-methylpyrrolidin-2-one"),
    ("O=C1CCCCN1CC",       "N-ethylpiperidin-2-one"),
    # 8-membered azacyclooctane
    ("CN1CCCCCCC1",        "N-methylazocane"),
    # Piperazine (two N atoms): N and N'
    ("CN1CCNCC1",          "N-methylpiperazine"),
    ("C1CN(C)CCN1C",       "N,N'-dimethylpiperazine"),
    # C-substituents still use numeric locants
    ("CC1CCCCN1",          "2-methylpiperidine"),
    # Aromatic N-heterocycles unchanged (numeric locants)
    ("Cn1cccc1",           "1-methylpyrrole"),
    ("Cn1cccn1",           "1-methylpyrazole"),
    # Regressions: unsubstituted rings unchanged
    ("C1CCCCN1",           "piperidine"),
    ("C1CCCN1",            "pyrrolidine"),
    ("C1CNCCO1",           "morpholine"),
    ("O=C1NC(=O)CC1",      "pyrrolidine-2,5-dione"),
])
def test_phase367_n_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
