"""Phase 135: α-アミノ酸保留名 (IUPAC 2013 P-12.1)

glycine, alanine, valine, leucine, isoleucine, proline,
phenylalanine, tryptophan, methionine, serine, threonine, cysteine,
tyrosine, asparagine, glutamine, lysine, arginine, histidine,
aspartic acid, glutamic acid (L/D および無立体中心形)
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # glycine (no stereo center)
    ("NCC(=O)O", "glycine"),
    # alanine
    ("N[C@@H](C)C(=O)O", "L-alanine"),
    ("N[C@H](C)C(=O)O", "D-alanine"),
    ("NC(C)C(=O)O", "alanine"),
    # valine
    ("N[C@@H](C(C)C)C(=O)O", "L-valine"),
    ("N[C@H](C(C)C)C(=O)O", "D-valine"),
    # leucine
    ("N[C@@H](CC(C)C)C(=O)O", "L-leucine"),
    # isoleucine
    ("N[C@@H]([C@@H](C)CC)C(=O)O", "L-isoleucine"),
    # proline
    ("OC(=O)[C@@H]1CCCN1", "L-proline"),
    ("OC(=O)[C@H]1CCCN1", "D-proline"),
    # phenylalanine
    ("N[C@@H](Cc1ccccc1)C(=O)O", "L-phenylalanine"),
    ("N[C@H](Cc1ccccc1)C(=O)O", "D-phenylalanine"),
    # tryptophan
    ("N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O", "L-tryptophan"),
    # methionine
    ("N[C@@H](CCSC)C(=O)O", "L-methionine"),
    # serine
    ("N[C@@H](CO)C(=O)O", "L-serine"),
    ("N[C@H](CO)C(=O)O", "D-serine"),
    # threonine
    ("N[C@@H]([C@@H](O)C)C(=O)O", "L-threonine"),
    # cysteine
    ("N[C@@H](CS)C(=O)O", "L-cysteine"),
    # tyrosine
    ("N[C@@H](Cc1ccc(O)cc1)C(=O)O", "L-tyrosine"),
    # asparagine
    ("N[C@@H](CC(=O)N)C(=O)O", "L-asparagine"),
    # glutamine
    ("N[C@@H](CCC(=O)N)C(=O)O", "L-glutamine"),
    # lysine
    ("N[C@@H](CCCCN)C(=O)O", "L-lysine"),
    ("N[C@H](CCCCN)C(=O)O", "D-lysine"),
    # arginine
    ("N[C@@H](CCCNC(=N)N)C(=O)O", "L-arginine"),
    # histidine
    ("N[C@@H](Cc1cnc[nH]1)C(=O)O", "L-histidine"),
    # aspartic acid
    ("N[C@@H](CC(=O)O)C(=O)O", "L-aspartic acid"),
    ("N[C@H](CC(=O)O)C(=O)O", "D-aspartic acid"),
    ("NC(CC(=O)O)C(=O)O", "aspartic acid"),
    # glutamic acid
    ("N[C@@H](CCC(=O)O)C(=O)O", "L-glutamic acid"),
    ("N[C@H](CCC(=O)O)C(=O)O", "D-glutamic acid"),
    ("NC(CCC(=O)O)C(=O)O", "glutamic acid"),
    # 回帰: non-amino-acid compounds unchanged
    ("CC(=O)O", "acetic acid"),
    ("CCCC(N)=O", "butanamide"),
    ("c1ccccc1", "benzene"),
])
def test_phase135_amino_acids(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
