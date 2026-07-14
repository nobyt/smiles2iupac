"""Phase 244: organobismuth/stibane/plumbane hydrides and Se/Te oxyacids (IUPAC 2013 P-68.1, P-65.3.3).

Group 15 organometallics:
  R_nBiH_{3-n}  → {alkyl}bismuthane
  R_nSbH_{3-n}  → {alkyl}stibane
  R_nPbH_{4-n}  → {alkyl}plumbane

Selenium/tellurium oxyacids:
  R-Se-OH           → {stem}aneselenenic acid
  R-Se(=O)-OH       → {stem}aneseleninic acid
  R-Se(=O)2-OH      → {stem}aneselenonic acid
  (Te analogues follow the same pattern)

Phase 857: R3As=O (arsane oxide, arsenic analog of phosphine oxide).
Previously matched no branch in the arsenic detection chain at all,
producing a garbled fallback name instead of a real arsane-oxide name.

Phase 858: R3As=S, R3Bi=O/=S, R3Sb=O/=S -- same missing-branch bug found
by probing the rest of the group 15/16 element family. R3Pb=O/=S have no
valid analog (would exceed lead's normal valence) so are not covered.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # organobismuthanes
    ("C[BiH2]",          "methylbismuthane"),
    ("CC[BiH2]",         "ethylbismuthane"),
    ("CC[BiH]CC",        "diethylbismuthane"),
    # organostibanes
    ("C[SbH2]",          "methylstibane"),
    ("CC[SbH2]",         "ethylstibane"),
    # organoplumbanes
    ("C[PbH3]",          "methylplumbane"),
    ("CC[PbH3]",         "ethylplumbane"),
    # arsane oxides/sulfides (Phase 857/858)
    ("C[As](C)(C)=O",    "trimethylarsane oxide"),
    ("CC[As](CC)(CC)=O", "triethylarsane oxide"),
    ("C[As](C)(C)=S",    "trimethylarsane sulfide"),
    ("CC[As](CC)(CC)=S", "triethylarsane sulfide"),
    # regression: plain arsane unchanged
    ("C[As](C)C",        "trimethylarsane"),
    # bismuthane oxides/sulfides (Phase 858)
    ("C[Bi](C)(C)=O",    "trimethylbismuthane oxide"),
    ("CC[Bi](CC)(CC)=O", "triethylbismuthane oxide"),
    ("C[Bi](C)(C)=S",    "trimethylbismuthane sulfide"),
    # regression: plain bismuthane (4 substituents) unchanged
    ("C[Bi](C)(C)C",     "tetramethylbismuthane"),
    # stibane oxides/sulfides (Phase 858)
    ("C[Sb](C)(C)=O",    "trimethylstibane oxide"),
    ("CC[Sb](CC)(CC)=O", "triethylstibane oxide"),
    ("C[Sb](C)(C)=S",    "trimethylstibane sulfide"),
    # regression: plain stibane unchanged
    ("C[Sb](C)C",        "trimethylstibane"),
    # selenium oxyacids
    ("C[Se]O",           "methaneselenenic acid"),
    ("CC[Se]O",          "ethaneselenenic acid"),
    ("C[Se](=O)O",       "methaneseleninic acid"),
    ("CC[Se](=O)O",      "ethaneseleninic acid"),
    ("CCC[Se](=O)O",     "propane-1-seleninic acid"),
    ("C[Se](=O)(=O)O",   "methaneselenonic acid"),
    # tellurium oxyacids
    ("C[Te]O",           "methanetellurenic acid"),
    ("C[Te](=O)O",       "methanetellurinic acid"),
    ("C[Te](=O)(=O)O",   "methanetelluronic acid"),
    # regression: germane/stannane unchanged
    ("C[GeH3]",          "methylgermane"),
    ("C[SnH3]",          "methylstannane"),
    # regression: selenol/selenide unchanged
    ("C[SeH]",           "methaneselenol"),
    ("C[Se]C",           "methylselanylmethane"),
])
def test_phase244_group15_and_chalcogen_oxyacids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
