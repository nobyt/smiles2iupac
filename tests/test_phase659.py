"""Phase 659: 2H-chromene and 4H-chromene methyl derivatives."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 2H-chromene (2H-1-benzopyran): O(1)-C2H2(sp3)-C3=C4-C4a-C(5-8)-C8a-O1
    # (parent in Phase 134; O1 not methylable)
    ("C1=Cc2ccccc2OC1",    "2H-chromene"),
    ("C1=Cc2ccccc2OC(C)1", "2-methyl-2H-chromene"),
    ("C(C)1=Cc2ccccc2OC1", "3-methyl-2H-chromene"),
    ("C1=C(C)c2ccccc2OC1", "4-methyl-2H-chromene"),
    ("C1=Cc2c(C)cccc2OC1", "5-methyl-2H-chromene"),
    ("C1=Cc2cc(C)ccc2OC1", "6-methyl-2H-chromene"),
    ("C1=Cc2ccc(C)cc2OC1", "7-methyl-2H-chromene"),
    ("C1=Cc2cccc(C)c2OC1", "8-methyl-2H-chromene"),
    # 4H-chromene (4H-1-benzopyran): O(1)-C2=C3-C4H2(sp3)-C4a-C(5-8)-C8a-O1
    # (parent in Phase 134; O1 not methylable)
    ("C1=COc2ccccc2C1",    "4H-chromene"),
    ("C1=C(C)Oc2ccccc2C1", "2-methyl-4H-chromene"),
    ("C(C)1=COc2ccccc2C1", "3-methyl-4H-chromene"),
    ("C1=COc2ccccc2C(C)1", "4-methyl-4H-chromene"),
    ("C1=COc2cccc(C)c2C1", "5-methyl-4H-chromene"),
    ("C1=COc2ccc(C)cc2C1", "6-methyl-4H-chromene"),
    ("C1=COc2cc(C)ccc2C1", "7-methyl-4H-chromene"),
    ("C1=COc2c(C)cccc2C1", "8-methyl-4H-chromene"),
])
def test_phase659(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
