"""Phase 272: partially unsaturated monocyclic heterocycles (IUPAC 2013 P-31.1.3).

5- and 6-membered rings with one heteroatom and at least one ring double bond
that are not aromatic — named as "X,Y-dihydro[parent]" or "X,Y,Z,W-tetrahydro[parent]".

Parent selection (IUPAC P-31.1.3.4):
  5-membered: furan, thiophene, pyrrole
  6-membered: pyridine, pyran, thiopyran
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ── 5-membered O (furan parent) ──────────────────────────────────────
    # sp3 at 2,3: O-CH2-CH2-CH=CH (db at 4-5)
    ("C1=CCCO1",      "2,3-dihydrofuran"),
    ("C1CC=CO1",      "2,3-dihydrofuran"),
    # sp3 at 2,5: O-CH2-CH=CH-CH2 (db at 3-4)
    ("C1C=CCO1",      "2,5-dihydrofuran"),

    # ── 5-membered S (thiophene parent) ──────────────────────────────────
    # sp3 at 2,3 (db at 4-5)
    ("C1=CCCS1",      "2,3-dihydrothiophene"),
    ("C1CC=CS1",      "2,3-dihydrothiophene"),
    # sp3 at 2,5 (db at 3-4)
    ("S1CC=CC1",      "2,5-dihydrothiophene"),
    ("C1C=CCS1",      "2,5-dihydrothiophene"),

    # ── 6-membered N (pyridine parent, 2 dbs remain) ─────────────────────
    # 1,2-dihydropyridine: N(sp3,H) at 1, C(sp3) at 2, dbs at 3-4 and 5-6
    ("C1=CC=CCN1",    "1,2-dihydropyridine"),

    # ── 6-membered N (pyridine parent, 1 db remains) ─────────────────────
    # 1,2,3,4-tetrahydropyridine: sp3 at 1(N),2,3,4; db at 5-6
    ("C1=CCCCN1",     "1,2,3,4-tetrahydropyridine"),
    # 1,2,3,6-tetrahydropyridine: sp3 at 1(N),2,3,6; db at 4-5
    ("N1CCC=CC1",     "1,2,3,6-tetrahydropyridine"),

    # ── regressions: fully saturated → Hantzsch-Widman / retained names ──
    ("C1CCCO1",       "oxolane"),
    ("C1CCCS1",       "thiolane"),
    ("C1CCNCC1",      "piperidine"),
    ("C1CCCCN1",      "piperidine"),

    # ── regressions: aromatic → retained names ────────────────────────────
    ("c1ccoc1",       "furan"),
    ("c1ccsc1",       "thiophene"),
    ("c1ccncc1",      "pyridine"),
])
def test_phase272_partial_unsat_heterocycles(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
