"""Phase 523: Heteroaromatic ring dicarboxylic acid naming (IUPAC 2013 P-65.1.1.4).

Furan, thiophene, pyrrole, and pyridine rings bearing two COOH groups
must be named as <ring>-X,Y-dicarboxylic acid, not as cycloalkane analogues.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # furan dicarboxylic acids
    ("OC(=O)c1ccc(C(=O)O)o1",   "furan-2,5-dicarboxylic acid"),
    ("OC(=O)c1cc(C(=O)O)co1",   "furan-2,4-dicarboxylic acid"),
    # thiophene dicarboxylic acids
    ("OC(=O)c1ccc(C(=O)O)s1",   "thiophene-2,5-dicarboxylic acid"),
    ("OC(=O)c1cc(C(=O)O)cs1",   "thiophene-2,4-dicarboxylic acid"),
    # pyrrole dicarboxylic acids (1H indicated hydrogen)
    ("OC(=O)c1ccc(C(=O)O)[nH]1",  "1H-pyrrole-2,5-dicarboxylic acid"),
    ("OC(=O)c1[nH]cc(C(=O)O)c1",  "1H-pyrrole-2,4-dicarboxylic acid"),
    ("OC(=O)c1cc[nH]c1C(=O)O",    "1H-pyrrole-2,3-dicarboxylic acid"),
    # pyridine dicarboxylic acids
    ("OC(=O)c1cccc(C(=O)O)n1",  "pyridine-2,6-dicarboxylic acid"),
    ("OC(=O)c1ccc(C(=O)O)cn1",  "pyridine-2,5-dicarboxylic acid"),
    ("OC(=O)c1ccnc(C(=O)O)c1",  "pyridine-2,4-dicarboxylic acid"),
    ("OC(=O)c1cncc(C(=O)O)c1",  "pyridine-3,5-dicarboxylic acid"),
    # regression: monocarboxylic acids unchanged
    ("OC(=O)c1ccco1",    "furan-2-carboxylic acid"),
    ("OC(=O)c1cccs1",    "thiophene-2-carboxylic acid"),
    ("OC(=O)c1ccc[nH]1", "1H-pyrrole-2-carboxylic acid"),
    ("OC(=O)c1ccccn1",   "pyridine-2-carboxylic acid"),
    # regression: benzene uses retained IUPAC 2013 names
    ("OC(=O)c1ccccc1C(=O)O",    "phthalic acid"),
    ("OC(=O)c1cccc(C(=O)O)c1",  "isophthalic acid"),
    ("OC(=O)c1ccc(C(=O)O)cc1",  "terephthalic acid"),
    # regression: cycloalkane dicarboxylic acid unchanged
    ("OC(=O)C1CCCCC1C(=O)O",    "cyclohexane-1,2-dicarboxylic acid"),
])
def test_phase523_heteroaromatic_diacid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
