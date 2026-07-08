"""Phase 850: N-methyl lactams of phthalazine/quinazoline/cinnoline/quinoxaline → -one suffix.

Extends Phase 849 oxo→-one conversion to four additional fused N-heterocycles.
When N is substituted (methyl etc.), the (nH) indicated hydrogen is omitted per P-14.7.1.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Phthalazinone: keto at C1, N-2 methyl → no (2H)
    ("Cn1ncc2ccccc2c1=O",         "2-methylphthalazin-1-one"),
    # Phthalazinone: NH form (already worked via hydroxy path)
    ("O=c1[nH]ncc2ccccc12",       "phthalazin-1(2H)-one"),
    # Quinazolinone: keto at C4, N-3 methyl → no (3H)
    ("Cn1cnc2ccccc2c1=O",         "3-methylquinazolin-4-one"),
    # Quinazolinone: NH form
    ("O=c1[nH]cnc2ccccc12",       "quinazolin-4(3H)-one"),
    # Cinnolinone: keto at C3, N-2 methyl → no (2H)
    ("Cn1nc2ccccc2cc1=O",         "2-methylcinnolin-3-one"),
    # Quinoxalinone: keto at C2, N-1 methyl → no (1H)
    ("Cn1c(=O)cnc2ccccc21",       "1-methylquinoxalin-2-one"),
    # Quinoxalinone: NH form (already worked)
    ("O=c1cnc2ccccc2[nH]1",       "quinoxalin-2(1H)-one"),
])
def test_phase850(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
