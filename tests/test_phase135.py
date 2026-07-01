"""Phase 135: α-アミノ酸 PIN (IUPAC 2013 P-102.4)

L-/D-アミノ酸の PIN は (R/S) 系統名。
glycine, alanine, valine, leucine, isoleucine, proline,
phenylalanine, tryptophan, methionine, serine, threonine, cysteine,
tyrosine, asparagine, glutamine, lysine, arginine, histidine,
aspartic acid, glutamic acid (L/D および無立体中心形)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # glycine → PIN: 2-aminoacetic acid (Phase 735)
    ("NCC(=O)O", "2-aminoacetic acid"),
    # alanine
    ("N[C@@H](C)C(=O)O", "(2S)-2-aminopropanoic acid"),
    ("N[C@H](C)C(=O)O", "(2R)-2-aminopropanoic acid"),
    ("NC(C)C(=O)O", "alanine"),
    # valine
    ("N[C@@H](C(C)C)C(=O)O", "(2S)-2-amino-3-methylbutanoic acid"),
    ("N[C@H](C(C)C)C(=O)O", "(2R)-2-amino-3-methylbutanoic acid"),
    # leucine
    ("N[C@@H](CC(C)C)C(=O)O", "(2S)-2-amino-4-methylpentanoic acid"),
    # isoleucine
    ("N[C@@H]([C@@H](C)CC)C(=O)O", "(2S,3S)-2-amino-3-methylpentanoic acid"),
    # proline
    ("OC(=O)[C@@H]1CCCN1", "(2S)-pyrrolidine-2-carboxylic acid"),
    ("OC(=O)[C@H]1CCCN1", "(2R)-pyrrolidine-2-carboxylic acid"),
    # phenylalanine
    ("N[C@@H](Cc1ccccc1)C(=O)O", "(2S)-2-amino-3-phenylpropanoic acid"),
    ("N[C@H](Cc1ccccc1)C(=O)O", "(2R)-2-amino-3-phenylpropanoic acid"),
    # tryptophan
    ("N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O", "(2S)-2-amino-3-(1H-indol-3-yl)propanoic acid"),
    # methionine
    ("N[C@@H](CCSC)C(=O)O", "(2S)-2-amino-4-methylsulfanylbutanoic acid"),
    # serine
    ("N[C@@H](CO)C(=O)O", "(2S)-2-amino-3-hydroxypropanoic acid"),
    ("N[C@H](CO)C(=O)O", "(2R)-2-amino-3-hydroxypropanoic acid"),
    # threonine
    ("N[C@@H]([C@@H](O)C)C(=O)O", "(2S,3S)-2-amino-3-hydroxybutanoic acid"),
    # cysteine
    ("N[C@@H](CS)C(=O)O", "(2R)-2-amino-3-sulfanylpropanoic acid"),
    # tyrosine
    ("N[C@@H](Cc1ccc(O)cc1)C(=O)O", "(2S)-2-amino-3-(4-hydroxyphenyl)propanoic acid"),
    # asparagine
    ("N[C@@H](CC(=O)N)C(=O)O", "(2S)-2,4-diamino-4-oxobutanoic acid"),
    # glutamine
    ("N[C@@H](CCC(=O)N)C(=O)O", "(2S)-2,5-diamino-5-oxopentanoic acid"),
    # lysine
    ("N[C@@H](CCCCN)C(=O)O", "(2S)-2,6-diaminohexanoic acid"),
    ("N[C@H](CCCCN)C(=O)O", "(2R)-2,6-diaminohexanoic acid"),
    # arginine
    ("N[C@@H](CCCNC(=N)N)C(=O)O", "(2S)-2-amino-5-(diaminomethylideneamino)pentanoic acid"),
    # histidine
    ("N[C@@H](Cc1cnc[nH]1)C(=O)O", "(2S)-2-amino-3-(1H-imidazol-5-yl)propanoic acid"),
    # aspartic acid
    ("N[C@@H](CC(=O)O)C(=O)O", "(2S)-2-aminobutanedioic acid"),
    ("N[C@H](CC(=O)O)C(=O)O", "(2R)-2-aminobutanedioic acid"),
    ("NC(CC(=O)O)C(=O)O", "2-aminobutanedioic acid"),
    # glutamic acid → PIN: 2-aminopentanedioic acid (Phase 735)
    ("N[C@@H](CCC(=O)O)C(=O)O", "(2S)-2-aminopentanedioic acid"),
    ("N[C@H](CCC(=O)O)C(=O)O", "(2R)-2-aminopentanedioic acid"),
    ("NC(CCC(=O)O)C(=O)O", "2-aminopentanedioic acid"),
    # 回帰: non-amino-acid compounds unchanged
    ("CC(=O)O", "acetic acid"),
    ("CCCC(N)=O", "butanamide"),
    ("c1ccccc1", "benzene"),
])
def test_phase135_amino_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
