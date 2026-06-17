"""Phase 552: Locant maps for quinoxaline, quinazoline, phthalazine, cinnoline.
These bicyclic diazine ring systems share the same 10-atom locant map structure
as quinoline/isoquinoline, enabling correct substituent locant assignment.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # quinoxaline (1,4-benzodiazine)
    ("c1ccc2nccnc2c1",         "quinoxaline"),
    ("c1ccc2nc(C)cnc2c1",     "2-methylquinoxaline"),
    ("c1ccc2ncc(C)nc2c1",     "3-methylquinoxaline"),
    ("c1ccc2nccnc2c1",         "quinoxaline"),  # sanity
    # quinazoline (1,3-benzodiazine)
    ("c1ccc2ncncc2c1",         "quinazoline"),
    ("c1ccc2nc(C)ncc2c1",     "2-methylquinazoline"),
    ("c1ccc2ncnc(C)c2c1",     "4-methylquinazoline"),
    # phthalazine (2,3-benzodiazine)
    ("c1ccc2cnncc2c1",         "phthalazine"),
    ("c1ccc2c(C)nncc2c1",     "1-methylphthalazine"),
    ("c1ccc2cnnc(C)c2c1",     "4-methylphthalazine"),
    # cinnoline (1,2-benzodiazine)
    ("c1ccc2nnccc2c1",         "cinnoline"),
    ("c1ccc2nnc(C)cc2c1",     "3-methylcinnoline"),
    ("c1ccc2nnccc2c1",         "cinnoline"),  # sanity
])
def test_phase552_bicyclic_diazines(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
