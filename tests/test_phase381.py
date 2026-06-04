"""Phase 381: Fix asymmetric 1,2-disubstituted hydrazine locant format.

'1,2-methyl-2-propylhydrazine' → '1-methyl-2-propylhydrazine'
The '1,2-' prefix was wrong; each substituent takes its own locant.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # asymmetric 1,2-disubstituted hydrazines
    ("CNNCCC",   "1-methyl-2-propylhydrazine"),
    ("CNNCC",    "1-ethyl-2-methylhydrazine"),
    ("CNNCCCC",  "1-butyl-2-methylhydrazine"),
    ("CCNNCCC",  "1-ethyl-2-propylhydrazine"),
    # symmetric cases must still work
    ("CNNC",     "1,2-dimethylhydrazine"),
    ("CCNNCC",   "1,2-diethylhydrazine"),
    # mono-substituted still correct
    ("CNN",      "methylhydrazine"),
    ("NN",       "hydrazine"),
    # 1,1-disubstituted still correct
    ("CN(N)C",   "1,1-dimethylhydrazine"),
])
def test_phase381_hydrazine_locant(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
