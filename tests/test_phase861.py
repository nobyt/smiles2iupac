"""Phase 861: selenocarboxylic acids -- the selenium analogues of the thioic
acids (Phase 149).

The thioic-acid machinery already names:
    C(=O)-SH   -> {stem}anethioic S-acid
    C(=S)-OH   -> {stem}anethioic O-acid
    C(=S)-SH   -> {stem}anedithioic acid
but the selenium analogues fell through to a generic (garbled) name such as
"1-oxoethaneselenol". This phase mirrors the pattern to selenium:
    C(=O)-SeH  -> {stem}aneselenoic Se-acid
    C(=Se)-OH  -> {stem}aneselenoic O-acid
    C(=Se)-SeH -> {stem}anediselenoic acid

Implementation reuses the existing _name_thioic_acid namer (chain finding,
multiple-bond, aromatic "benzenecarbo..." handling) by table-driving the acid
word on group_type, so aromatic and unsaturated cases work for free.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # aliphatic selenoic Se-acid  C(=O)-SeH
    ("C(=O)[SeH]",          "methaneselenoic Se-acid"),
    ("CC(=O)[SeH]",         "ethaneselenoic Se-acid"),
    ("CCC(=O)[SeH]",        "propaneselenoic Se-acid"),
    # selenoic O-acid  C(=Se)-OH
    ("CC(=[Se])O",          "ethaneselenoic O-acid"),
    ("CCC(=[Se])O",         "propaneselenoic O-acid"),
    # diselenoic acid  C(=Se)-SeH
    ("CC(=[Se])[SeH]",      "ethanediselenoic acid"),
    ("CCC(=[Se])[SeH]",     "propanediselenoic acid"),
    # unsaturated (exercises the multiple-bond path)
    ("C=CC(=O)[SeH]",       "prop-2-eneselenoic Se-acid"),
    # aromatic carbo-selenoic (benzene ring directly attached)
    ("[SeH]C(=O)c1ccccc1",  "benzenecarboselenoic Se-acid"),
    ("OC(=[Se])c1ccccc1",   "benzenecarboselenoic O-acid"),
    # regression: thioic acids unchanged
    ("CC(=O)S",             "ethanethioic S-acid"),
    ("CC(=S)O",             "ethanethioic O-acid"),
    ("CC(=S)S",             "ethanedithioic acid"),
    ("SC(=O)c1ccccc1",      "benzenecarbothioic S-acid"),
    ("C=CC(=O)S",           "prop-2-enethioic S-acid"),
])
def test_phase861_selenocarboxylic_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
