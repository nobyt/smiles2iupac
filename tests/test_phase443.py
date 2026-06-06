"""Phase 443: pyrazolo/imidazo/pyrrolo-quinoline and pyrazolo-quinoxaline retained names
(IUPAC 2013 P-31.1.3).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-pyrazolo[3,4-b]quinoline
    ("c1ccc2nc3[nH]ncc3cc2c1",  "1H-pyrazolo[3,4-b]quinoline"),
    # pyrazolo[1,5-a]quinoline
    ("c1ccc2nc3ccnn3cc2c1",     "pyrazolo[1,5-a]quinoline"),
    # pyrazolo[1,5-a]quinoxaline
    ("c1ccc2nn3nccc3nc2c1",     "pyrazolo[1,5-a]quinoxaline"),
    # 1H-imidazo[4,5-b]quinoline
    ("c1ccc2nc3[nH]cnc3cc2c1",  "1H-imidazo[4,5-b]quinoline"),
    # 3H-imidazo[4,5-b]quinoline (tautomer)
    ("c1ccc2nc3nc[nH]c3cc2c1",  "3H-imidazo[4,5-b]quinoline"),
    # 1H-pyrrolo[3,2-b]quinoline
    ("c1ccc2nc3cc[nH]c3cc2c1",  "1H-pyrrolo[3,2-b]quinoline"),
    # regressions
    ("c1ccc2nc3[nH]ncc3cc2c1",  "1H-pyrazolo[3,4-b]quinoline"),
    ("c1ccc2nc3ccnn3cc2c1",     "pyrazolo[1,5-a]quinoline"),
    ("c1cnc2ncccc2c1",          "1,8-naphthyridine"),
    ("c1ccc2nc3cc[nH]c3cc2c1",  "1H-pyrrolo[3,2-b]quinoline"),
])
def test_phase443_quinoline_fused_hetero(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
