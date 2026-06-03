"""Phase 175: isocyanate/isothiocyanate substitutive PIN (IUPAC 2013 P-65.3.1).

R-N=C=O → isocyanato{alkane} (PIN; "{alkyl} isocyanate" is retained acceptable):
  CN=C=O   → isocyanatomethane
  CCN=C=O  → isocyanatoethane
  aromatic: isocyanatobenzene (substitutive, already correct)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 脂肪族: 置換命名 PIN
    ("CN=C=O",           "isocyanatomethane"),
    ("CCN=C=O",          "isocyanatoethane"),
    ("CCCN=C=O",         "isocyanatopropane"),
    ("CCCCN=C=O",        "isocyanatobutane"),
    ("CC(C)N=C=O",       "2-isocyanatopropane"),
    # 芳香族: 置換命名 (既存動作維持)
    ("c1ccccc1N=C=O",    "isocyanatobenzene"),
    # 回帰: イソチオシアネート
    ("CN=C=S",           "isothiocyanatomethane"),
    ("CCN=C=S",          "isothiocyanatoethane"),
])
def test_phase175_isocyanate_substitutive(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
