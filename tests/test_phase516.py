"""Phase 516: isocyanate/isothiocyanate functional-class PIN (IUPAC 2013 P-65.3.1)

R-N=C=O → "{R} isocyanate"; R-N=C=S → "{R} isothiocyanate".
Previous implementation used substitutive prefix form "isocyanato{alkane}".
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # aliphatic isocyanates
    ("CN=C=O",          "methyl isocyanate"),
    ("CCN=C=O",         "ethyl isocyanate"),
    ("CCCN=C=O",        "propyl isocyanate"),
    ("CC(C)N=C=O",      "propan-2-yl isocyanate"),
    # aromatic isocyanates
    ("O=C=Nc1ccccc1",   "phenyl isocyanate"),
    ("O=C=Nc1ccc(C)cc1","4-methylphenyl isocyanate"),
    # aliphatic isothiocyanates
    ("CN=C=S",          "methyl isothiocyanate"),
    ("CCN=C=S",         "ethyl isothiocyanate"),
    ("CC(C)N=C=S",      "propan-2-yl isothiocyanate"),
    # aromatic isothiocyanates
    ("S=C=Nc1ccccc1",   "phenyl isothiocyanate"),
    # E/Z unsaturated R group
    ("C/C=C/CN=C=O",    "(2E)-but-2-en-1-yl isocyanate"),
    # regression: diisocyanate uses substitutive prefix (unchanged)
    ("O=C=NCCN=C=O",    "1,2-diisocyanatoethane"),
    # regression: isocyanic acid unaffected
    ("N=C=O",           "isocyanic acid"),
])
def test_phase516_isocyanate_functional_class(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
