"""Phase 162/175: isocyanate/isothiocyanate functional-class PIN (IUPAC 2013 P-65.3.1).

PIN uses functional-class name "{alkyl} isocyanate" / "{alkyl} isothiocyanate".
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # isocyanate 脂肪族
    ("CN=C=O",          "methyl isocyanate"),
    ("CCN=C=O",         "ethyl isocyanate"),
    ("CCCN=C=O",        "propyl isocyanate"),
    ("CCCCN=C=O",       "butyl isocyanate"),
    ("CC(C)N=C=O",      "propan-2-yl isocyanate"),
    # isocyanate 芳香族
    ("O=C=Nc1ccccc1",   "phenyl isocyanate"),
    ("O=C=Nc1ccc(C)cc1","4-methylphenyl isocyanate"),
    # isothiocyanate 脂肪族
    ("CN=C=S",          "methyl isothiocyanate"),
    ("CCN=C=S",         "ethyl isothiocyanate"),
    ("CCCN=C=S",        "propyl isothiocyanate"),
    ("CC(C)N=C=S",      "propan-2-yl isothiocyanate"),
    ("CCCCN=C=S",       "butyl isothiocyanate"),
    # isothiocyanate 芳香族
    ("S=C=Nc1ccccc1",   "phenyl isothiocyanate"),
    # 回帰: アミド・チオアミドは変わらない
    ("CC(=O)N",         "acetamide"),
    ("CC(=S)N",         "ethanethioamide"),
])
def test_phase162_isocyanate_functional_class(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
