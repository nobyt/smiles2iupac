"""Phase 162/175: isocyanate/isothiocyanate substitutive PIN (IUPAC 2013 P-65.3.1).

PIN uses substitutive prefix "isocyanato-" / "isothiocyanato-" on parent alkane.
"{alkyl} isocyanate" / "{alkyl} isothiocyanate" are retained acceptable but not PIN.
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isocyanate 脂肪族 (置換命名 PIN)
    ("CN=C=O",          "isocyanatomethane"),
    ("CCN=C=O",         "isocyanatoethane"),
    ("CCCN=C=O",        "isocyanatopropane"),
    ("CCCCN=C=O",       "isocyanatobutane"),
    ("CC(C)N=C=O",      "2-isocyanatopropane"),
    # isocyanate 芳香族 (置換命名)
    ("O=C=Nc1ccccc1",   "isocyanatobenzene"),
    ("O=C=Nc1ccc(C)cc1","1-isocyanato-4-methylbenzene"),
    # isothiocyanate 脂肪族 (置換命名 PIN)
    ("CN=C=S",          "isothiocyanatomethane"),
    ("CCN=C=S",         "isothiocyanatoethane"),
    ("CCCN=C=S",        "isothiocyanatopropane"),
    ("CC(C)N=C=S",      "2-isothiocyanatopropane"),
    ("CCCCN=C=S",       "isothiocyanatobutane"),
    # isothiocyanate 芳香族
    ("S=C=Nc1ccccc1",   "isothiocyanatobenzene"),
    # 回帰: アミド・チオアミドは変わらない
    ("CC(=O)N",         "acetamide"),
    ("CC(=S)N",         "ethanethioamide"),
])
def test_phase162_isocyanate_substitutive(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
