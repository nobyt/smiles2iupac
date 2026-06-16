"""Phase 560: Substituted phenazine and 1,10-phenanthroline naming.
Phenazine has D2h symmetry → only 2 unique positions (locants 1 and 2).
1,10-phenanthroline has C2 symmetry → 4 unique positions (locants 2–5).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phenazine
    ("c1ccc2nc3ccccc3nc2c1",    "phenazine"),
    ("Cc1cccc2nc3ccccc3nc12",   "1-methylphenazine"),
    ("Cc1ccc2nc3ccccc3nc2c1",   "2-methylphenazine"),
    # 1,10-phenanthroline (C2 symmetric: 2≡9, 3≡8, 4≡7, 5≡6)
    ("c1cnc2c(c1)ccc1cccnc12",  "1,10-phenanthroline"),
    ("Cc1ccc2ccc3cccnc3c2n1",   "2-methyl-1,10-phenanthroline"),
    ("Cc1cnc2c(ccc3cccnc32)c1", "3-methyl-1,10-phenanthroline"),
    ("Cc1ccnc2c1ccc1cccnc12",   "4-methyl-1,10-phenanthroline"),
    ("Cc1cc2cccnc2c2ncccc12",   "5-methyl-1,10-phenanthroline"),
])
def test_phase560_phenazine_110phen(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
