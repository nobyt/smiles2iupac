"""Phase 847: IUPAC 2013 PINs for methylated purine-2,6-diones (caffeine, theophylline, theobromine).

P-14.7.1: indicated H and dihydro prefixes are omitted when the N at that locant bears a substituent.
For purine-2,6-dione ring system: N1, N3, N7 may each bear a methyl; drop the corresponding
dihydro/indicated-H tokens from the retained ring base name.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Caffeine: N1, N3, N7 all methylated — drop 3,7-dihydro and 1H- entirely
    ("Cn1c(=O)c2c(ncn2C)n(C)c1=O",  "1,3,7-trimethylpurine-2,6-dione"),
    # Theophylline: N1, N3 methylated; N7-H remains → 7H- prefix form
    ("Cn1c(=O)c2[nH]cnc2n(C)c1=O",  "1,3-dimethyl-7H-purine-2,6-dione"),
    # Theobromine: N3, N7 methylated; N1-H remains (lactam NH, implicit in 2,6-dione)
    ("Cn1cnc2c1c(=O)[nH]c(=O)n2C",  "3,7-dimethylpurine-2,6-dione"),
    # Xanthine (unsubstituted): 3,7-dihydro unchanged
    ("O=c1[nH]c(=O)c2[nH]cnc2[nH]1", "3,7-dihydropurine-2,6-dione"),
])
def test_phase847(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
