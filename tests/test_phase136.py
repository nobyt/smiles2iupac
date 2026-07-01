"""Phase 136: 二塩基酸・ヒドロキシ酸 保留名 → IUPAC 2013 PINs (P-65.1.1.4, P-65.1.3)

oxalic acid, malonic acid, adipic acid (retained PINs),
butanedioic/pentanedioic/heptanedioic/octanedioic/nonanedioic/decanedioic acid (was succinic..sebacic),
maleic acid, (E)-but-2-enedioic acid (was fumaric acid),
(E/Z)-3-phenylprop-2-enoic acid (was cinnamic acid),
phthalic acid, benzene-1,3-dicarboxylic acid (was isophthalic acid), terephthalic acid,
2-hydroxypropanoic acid (was lactic acid, L/D/racemic),
2-hydroxybutanedioic acid (was malic acid, L/D/racemic),
2,3-dihydroxybutanedioic acid (was tartaric acid, L/D/meso/racemic), citric acid
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # aliphatic diacids (oxalic/malonic/adipic: retained PINs; others: systematic PINs)
    ("OC(=O)C(=O)O",         "oxalic acid"),
    ("OC(=O)CC(=O)O",        "malonic acid"),
    ("OC(=O)CCC(=O)O",       "butanedioic acid"),
    ("OC(=O)CCCC(=O)O",      "pentanedioic acid"),
    ("OC(=O)CCCCC(=O)O",     "adipic acid"),
    ("OC(=O)CCCCCC(=O)O",    "heptanedioic acid"),
    ("OC(=O)CCCCCCC(=O)O",   "octanedioic acid"),
    ("OC(=O)CCCCCCCC(=O)O",  "nonanedioic acid"),
    ("OC(=O)CCCCCCCCC(=O)O", "decanedioic acid"),
    # unsaturated diacids
    ("OC(=O)/C=C\\C(=O)O",   "maleic acid"),
    ("OC(=O)/C=C/C(=O)O",    "(E)-but-2-enedioic acid"),
    # cinnamic acid → PIN: (E/Z)-3-phenylprop-2-enoic acid
    ("OC(=O)/C=C/c1ccccc1",  "(E)-3-phenylprop-2-enoic acid"),
    ("OC(=O)/C=C\\c1ccccc1", "(Z)-3-phenylprop-2-enoic acid"),
    # aromatic diacids
    ("OC(=O)c1ccccc1C(=O)O",       "phthalic acid"),
    ("OC(=O)c1cccc(C(=O)O)c1",     "benzene-1,3-dicarboxylic acid"),
    ("OC(=O)c1ccc(C(=O)O)cc1",     "terephthalic acid"),
    # hydroxy acids → PINs (Phase 735)
    ("C[C@@H](O)C(=O)O",    "(2R)-2-hydroxypropanoic acid"),
    ("C[C@H](O)C(=O)O",     "(2S)-2-hydroxypropanoic acid"),
    ("CC(O)C(=O)O",          "2-hydroxypropanoic acid"),
    ("OC(=O)[C@@H](O)CC(=O)O", "(2S)-2-hydroxybutanedioic acid"),
    ("OC(=O)[C@H](O)CC(=O)O",  "(2R)-2-hydroxybutanedioic acid"),
    ("OC(=O)C(O)CC(=O)O",      "2-hydroxybutanedioic acid"),
    ("OC(=O)[C@@H](O)[C@H](O)C(=O)O",  "(2S,3S)-2,3-dihydroxybutanedioic acid"),
    ("OC(=O)[C@H](O)[C@@H](O)C(=O)O",  "(2R,3R)-2,3-dihydroxybutanedioic acid"),
    ("OC(=O)[C@@H](O)[C@@H](O)C(=O)O", "(2R,3S)-2,3-dihydroxybutanedioic acid"),
    ("OC(=O)C(O)C(O)C(=O)O",            "2,3-dihydroxybutanedioic acid"),
    # citric acid
    ("OC(=O)CC(O)(CC(=O)O)C(=O)O", "citric acid"),
    # 回帰: glycine → PIN 2-aminoacetic acid (Phase 735)
    ("NCC(=O)O",        "2-aminoacetic acid"),
    ("C[C@H](N)C(=O)O", "(2S)-2-aminopropanoic acid"),
    # 回帰: simple acids unchanged
    ("CC(=O)O",  "acetic acid"),
    ("CCCC(=O)O","butanoic acid"),
    # 回帰: 2-ハロゲン dioic acids are NOT retained names (Phase 76)
    ("ClC(CC(=O)O)C(=O)O", "2-chlorobutanedioic acid"),
])
def test_phase136_acid_retained(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
