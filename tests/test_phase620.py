"""Phase 620: non-aromatic N=C fused heterocycles — 1H-pyrrolo[3,4-c]pyridine,
1H-pyrrolo[3,4-d]pyridazine, 3H-imidazo[4,5-c]pyridazine with methyl substituents."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-pyrrolo[3,4-c]pyridine (1H at sp3 C1; positions 1,3,4,6,7)
    ("C1=NCc2ccncc21",    "1H-pyrrolo[3,4-c]pyridine"),
    ("CC1N=Cc2cnccc21",   "1-methyl-1H-pyrrolo[3,4-c]pyridine"),
    ("CC1=NCc2ccncc21",   "3-methyl-1H-pyrrolo[3,4-c]pyridine"),
    ("Cc1nccc2c1C=NC2",   "4-methyl-1H-pyrrolo[3,4-c]pyridine"),
    ("Cc1cc2c(cn1)C=NC2", "6-methyl-1H-pyrrolo[3,4-c]pyridine"),
    ("Cc1cncc2c1CN=C2",   "7-methyl-1H-pyrrolo[3,4-c]pyridine"),
    # 1H-pyrrolo[3,4-d]pyridazine (1H at sp3 C1; positions 1,4,5,7)
    ("C1=NC=C2CN=NC=C12",      "1H-pyrrolo[3,4-d]pyridazine"),
    ("CC1N=NC=C2C=NC=C21",     "1-methyl-1H-pyrrolo[3,4-d]pyridazine"),
    ("CC1=C2C=NC=C2CN=N1",     "4-methyl-1H-pyrrolo[3,4-d]pyridazine"),
    ("CC1=NC=C2CN=NC=C21",     "5-methyl-1H-pyrrolo[3,4-d]pyridazine"),
    ("CC1=C2CN=NC=C2C=N1",     "7-methyl-1H-pyrrolo[3,4-d]pyridazine"),
    # 3H-imidazo[4,5-c]pyridazine (3H at sp3 C3; positions 3,4,6)
    ("C1=NC2=CCN=NC2=N1",  "3H-imidazo[4,5-c]pyridazine"),
    ("CC1C=C2N=CN=C2N=N1", "3-methyl-3H-imidazo[4,5-c]pyridazine"),
    ("CC1=C2N=CN=C2N=NC1", "4-methyl-3H-imidazo[4,5-c]pyridazine"),
    ("CC1=NC2=CCN=NC2=N1", "6-methyl-3H-imidazo[4,5-c]pyridazine"),
])
def test_phase620_nonaromatic_nc_fused(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
