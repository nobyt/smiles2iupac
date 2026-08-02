"""Phase 886: halide/hydroxide/amino substituents on organogermanes,
organostannanes, organobismuthanes, organostibanes, and organoplumbanes.

Found via the fresh probe sweep in the same session as Phase 884/885. Same
info-loss bug class as Phase 863 (boron halides) and Phase 883 (B/Si
amines): the germane_org/stannane_org/bismuthane_org/stibane_org/
plumbane_org namers only ever looked at C-neighbors of the central atom, so
ANY other substituent (halogen, -OH, -NH2) was silently dropped and
different real compounds collapsed to the same bare name:

  C[Ge](Cl)(Cl)Cl -> "methylgermane"  (same as plain CG, drops 3 Cl)
  C[Sn](O)(O)O     -> "methylstannane" (drops 3 OH)
  C[Bi](N)N        -> "methylbismuthane" (drops 2 NH2)
  C[Pb](Cl)(Cl)Cl  -> "methylplumbane" (drops 3 Cl)

Verified via OPSIN parse-back for every new compound family before
implementing.

Implementation:
- Halide + amino: the germane/stannane/bismuthane/stibane/plumbane namers
  now delegate to the SAME _name_borane_silane_with_halo_amino helper
  written for Phase 863/883 -- it was already fully generic (parameterized
  on the base word and the central atom), so no new halogen/amine logic was
  needed, just swapping which function each element's namer calls.
- Hydroxyl: new germanol_org/stannanol_org/bismuthanol_org/stibanol_org/
  plumbanol_org group types mirror the pre-existing silanol_org exactly
  (R_nEl(OH)_{m-n} -> element-ol/-diol/-triol/-tetraol), sharing a newly
  extracted _name_element_ol_family helper (silanol_org's own namer was
  refactored to call it too, behavior-preserving).
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # halides
    ("C[Ge](Cl)(Cl)Cl", "trichloro(methyl)germane"),
    ("C[Sn](Cl)(Cl)Cl", "trichloro(methyl)stannane"),
    ("C[Pb](Cl)(Cl)Cl", "trichloro(methyl)plumbane"),
    ("C[Bi](Cl)Cl",     "dichloro(methyl)bismuthane"),
    ("C[Sb](Cl)Cl",     "dichloro(methyl)stibane"),
    # amines
    ("C[Ge](N)(N)N", "methylgermanetriamine"),
    ("C[Sn](N)(N)N", "methylstannanetriamine"),
    ("C[Pb](N)(N)N", "methylplumbanetriamine"),
    ("C[Bi](N)N",    "methylbismuthanediamine"),
    ("C[Sb](N)N",    "methylstibanediamine"),
    # hydroxyls (ol/diol/triol family)
    ("C[Ge](O)(O)O", "methylgermanetriol"),
    ("C[Sn](O)(O)O", "methylstannanetriol"),
    ("C[Pb](O)(O)O", "methylplumbanetriol"),
    ("C[Bi](O)O",    "methylbismuthanediol"),
    ("C[Sb](O)O",    "methylstibanediol"),
    ("C[Ge](O)C",    "dimethylgermanol"),
    # regression: plain organo-El (no halide/OH/amino) unchanged
    ("C[Ge](C)(C)C", "tetramethylgermane"),
    ("C[Sn](C)(C)C", "tetramethylstannane"),
    ("C[Pb](C)(C)C", "tetramethylplumbane"),
    ("C[Bi](C)C",    "trimethylbismuthane"),
    ("C[Sb](C)C",    "trimethylstibane"),
    # regression: Bi/Sb oxide/sulfide (Phase 858) unchanged
    ("C[Bi](C)(C)=O", "trimethylbismuthane oxide"),
    ("C[Bi](C)(C)=S", "trimethylbismuthane sulfide"),
    ("C[Sb](C)(C)=O", "trimethylstibane oxide"),
    ("C[Sb](C)(C)=S", "trimethylstibane sulfide"),
    # regression: silanol_org (Phase 231/379) unchanged after the
    # _name_element_ol_family extraction
    ("C[Si](C)(C)O", "trimethylsilanol"),
    ("C[Si](C)(O)O", "dimethylsilanediol"),
])
def test_phase886_ge_sn_bi_sb_pb_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase886_not_confused_with_plain_methylgermane():
    assert smiles_to_iupac("C[Ge](Cl)(Cl)Cl") != smiles_to_iupac("C[Ge](C)(C)C")
    assert smiles_to_iupac("C[Sn](O)(O)O") != smiles_to_iupac("C[Sn](C)(C)C")
