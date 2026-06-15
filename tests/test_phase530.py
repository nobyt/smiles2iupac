"""Phase 530: Heteroaromatic carbothioic acid naming (IUPAC 2013 P-65.1.6).

Rings bearing -C(=O)SH or -C(=S)OH were returning "methanethioic S/O-acid"
because _name_thioic_acid only had a benzene hardcoded path, and
_collect_acid_chain yielded a 1-carbon chain for ring-adjacent carbonyls.
Now correctly named as "<ring>-carbothioic S/O-acid".
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Heteroaromatic carbothioic S-acids
    ("SC(=O)c1ccccn1",   "pyridine-2-carbothioic S-acid"),
    ("SC(=O)c1cccnc1",   "pyridine-3-carbothioic S-acid"),
    ("SC(=O)c1ccncc1",   "pyridine-4-carbothioic S-acid"),
    ("SC(=O)c1cccs1",    "thiophene-2-carbothioic S-acid"),
    ("SC(=O)c1ccco1",    "furan-2-carbothioic S-acid"),
    ("SC(=O)c1ccc[nH]1", "1H-pyrrole-2-carbothioic S-acid"),
    # Heteroaromatic carbothioic O-acids
    ("OC(=S)c1ccccn1",   "pyridine-2-carbothioic O-acid"),
    ("OC(=S)c1cccs1",    "thiophene-2-carbothioic O-acid"),
    # Regression: benzene unchanged
    ("SC(=O)c1ccccc1",   "benzenecarbothioic S-acid"),
    ("OC(=S)c1ccccc1",   "benzenecarbothioic O-acid"),
    # Regression: aliphatic unchanged
    ("SC(=O)C",          "ethanethioic S-acid"),
    ("OC(=S)C",          "ethanethioic O-acid"),
])
def test_phase530_heteroaromatic_carbothioic_acid(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
