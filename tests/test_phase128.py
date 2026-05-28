"""Phase 128: 保留名 guanidine・cyanic/thiocyanic/isocyanic acid (IUPAC 2013 P-66.4.1.3)"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # guanidine: aminomethanimidamide → guanidine (IUPAC 2013 PIN)
    ("NC(=N)N", "guanidine"),
    # isocyanic acid: H-N=C=O
    ("N=C=O", "isocyanic acid"),
    # cyanic acid: HO-C≡N
    ("OC#N", "cyanic acid"),
    ("N#CO", "cyanic acid"),
    # thiocyanic acid: HS-C≡N
    ("SC#N", "thiocyanic acid"),
    # 回帰: amidine unchanged
    ("CC(=N)N", "ethanimidamide"),
    # 回帰: isocyanate ester unchanged
    ("CN=C=O", "methyl isocyanate"),
    # 回帰: thiocyanate ester unchanged
    ("CSC#N", "methyl thiocyanate"),
    # 回帰: cyanate ester unchanged
    ("COC#N", "methyl cyanate"),
    # 回帰: thiol unchanged
    ("CCS", "ethanethiol"),
])
def test_phase128_retained_names(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
