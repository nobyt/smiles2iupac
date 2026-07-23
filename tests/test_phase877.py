"""Phase 877: biphenyl as a substituent -> [1,1'-biphenyl]-N-yl (stop ring drop).

When a biphenyl was a substituent on a chain/principal group (e.g. an aryl
ketone), the second ring was dropped entirely:
CC(=O)c1ccc(-c2ccccc2)cc1 was named "1-phenylethanone" -- identical to
acetophenone (CC(=O)c1ccccc1). Two different compounds collapsed to one name.

Root cause: the aryl-substituent BFS in _name_carbon_substituent traversed the
inter-ring single bond, collected 12 aromatic carbons, and (len != 6) fell back
to bare "phenyl". Fix: detect the biphenyl-as-substituent (_find_benzene_
biphenyl) and emit the [1,1'-biphenyl]-N-yl PIN, with the attachment carbon
numbered relative to the inter-ring carbon (para 4 / meta 3 / ortho 2). Only
pure (unsubstituted-ring) biphenyl substituents use this path.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # biphenylyl on a ketone / chain -> PIN substituent, no dropped ring
    ("CC(=O)c1ccc(-c2ccccc2)cc1",  "1-([1,1'-biphenyl]-4-yl)ethan-1-one"),
    ("CCC(=O)c1ccc(-c2ccccc2)cc1", "1-([1,1'-biphenyl]-4-yl)propan-1-one"),
    ("CC(=O)c1cccc(-c2ccccc2)c1",  "1-([1,1'-biphenyl]-3-yl)ethan-1-one"),
    ("CC(=O)c1ccccc1-c1ccccc1",    "1-([1,1'-biphenyl]-2-yl)ethan-1-one"),
    ("OCc1ccc(-c2ccccc2)cc1",      "([1,1'-biphenyl]-4-yl)methanol"),
    # the flagged defect: must NOT equal acetophenone
    # (checked explicitly below)

    # regression: acetophenone and simple substituted variants unchanged
    ("CC(=O)c1ccccc1",       "1-phenylethanone"),
    ("CC(=O)c1ccc(Cl)cc1",   "1-(4-chlorophenyl)ethan-1-one"),
    ("CC(=O)c1ccc(C)cc1",    "1-(4-methylphenyl)ethan-1-one"),
    ("OCc1ccccc1",           "phenylmethanol"),
    # regression: diphenylmethane (rings not directly bonded) unchanged
    ("C(c1ccccc1)c1ccccc1",  "(phenylmethyl)benzene"),
    ("OC(c1ccccc1)c1ccccc1", "diphenylmethanol"),
    # regression: biphenyl as parent (Phase 874-876) unchanged
    ("c1ccccc1-c1ccccc1",       "1,1'-biphenyl"),
    ("Cc1ccc(-c2ccc(C)cc2)cc1", "4,4'-dimethyl-1,1'-biphenyl"),
    ("OC(=O)c1ccc(-c2ccccc2)cc1", "[1,1'-biphenyl]-4-carboxylic acid"),
    ("c1ccc2ccccc2c1",       "naphthalene"),
])
def test_phase877_biphenylyl_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase877_not_confused_with_acetophenone():
    # 4-acetylbiphenyl and acetophenone must get distinct names.
    assert smiles_to_iupac("CC(=O)c1ccc(-c2ccccc2)cc1") != smiles_to_iupac("CC(=O)c1ccccc1")
