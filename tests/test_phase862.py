"""Phase 862: tellurocarboxylic acids -- the tellurium analogues of the thioic
(Phase 149) and selenoic (Phase 861) acids, completing the chalcogen set.

    C(=O)-TeH   -> {stem}anetelluroic Te-acid
    C(=Te)-OH   -> {stem}anetelluroic O-acid
    C(=Te)-TeH  -> {stem}aneditelluroic acid

Before this phase these fell through to garbled generic names such as
"1-oxoethanetellurol". Implementation mirrors Phase 861 exactly: detection in
the shared thiocarboxylic block plus three FunctionalGroupSpecs, reusing the
table-driven _name_thioic_acid namer so chain/unsaturation/aromatic paths work
for free. The codebase already mirrors tellurium across the chalcogen groups
(tellurol, telluramide, telluroxide, tellurone, telluronic acid), so this keeps
the carboxylic-acid family symmetric with sulfur and selenium.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # telluroic Te-acid  C(=O)-TeH
    ("C(=O)[TeH]",          "methanetelluroic Te-acid"),
    ("CC(=O)[TeH]",         "ethanetelluroic Te-acid"),
    ("CCC(=O)[TeH]",        "propanetelluroic Te-acid"),
    # telluroic O-acid  C(=Te)-OH
    ("CC(=[Te])O",          "ethanetelluroic O-acid"),
    ("CCC(=[Te])O",         "propanetelluroic O-acid"),
    # ditelluroic acid  C(=Te)-TeH
    ("CC(=[Te])[TeH]",      "ethaneditelluroic acid"),
    # unsaturated (multiple-bond path)
    ("C=CC(=O)[TeH]",       "prop-2-enetelluroic Te-acid"),
    # aromatic carbo-telluroic
    ("[TeH]C(=O)c1ccccc1",  "benzenecarbotelluroic Te-acid"),
    ("OC(=[Te])c1ccccc1",   "benzenecarbotelluroic O-acid"),
    # regression: thioic + selenoic unchanged
    ("CC(=O)S",             "ethanethioic S-acid"),
    ("CC(=S)S",             "ethanedithioic acid"),
    ("CC(=O)[SeH]",         "ethaneselenoic Se-acid"),
    ("CC(=[Se])[SeH]",      "ethanediselenoic acid"),
])
def test_phase862_tellurocarboxylic_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
