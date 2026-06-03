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
    ("C[Se]C",           "dimethyl selenide"),
])
def test_phase244_group15_and_chalcogen_oxyacids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
