"""Phase 384: prop-1-ene / prop-1-yne require locant '1' (IUPAC 2013 P-31.1.2.2);
acrylic acid / methacrylic acid are not IUPAC 2013 preferred names.

IUPAC 2013 preferred names always cite locants for multiple bonds,
even when there is no ambiguity: propene → prop-1-ene, propyne → prop-1-yne.

Also, acrylic acid and methacrylic acid are not retained preferred IUPAC names
(not in Table 65.1 of IUPAC 2013); their preferred names are prop-2-enoic acid
and 2-methylprop-2-enoic acid respectively.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # propene: locant 1 required
    ("C=CC",              "prop-1-ene"),
    ("CC=C",              "prop-1-ene"),
    # propyne: locant 1 required
    ("C#CC",              "prop-1-yne"),
    ("CC#C",              "prop-1-yne"),
    # longer chains still correct (regression)
    ("C=CCC",             "but-1-ene"),
    ("C#CCC",             "but-1-yne"),
    ("CC=CC",             "but-2-ene"),
    # ethene/ethyne: no locant (2-carbon chain; omitting is IUPAC rule for 2C)
    ("C=C",               "ethene"),
    ("C#C",               "ethyne"),
    # acrylic acid: preferred IUPAC name is systematic
    ("C=CC(=O)O",         "prop-2-enoic acid"),
    # methacrylic acid: preferred IUPAC name is systematic
    ("C=C(C)C(=O)O",      "2-methylprop-2-enoic acid"),
])
def test_phase384_propene_propyne_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
