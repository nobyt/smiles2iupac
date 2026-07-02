"""Phase 754: 5-membered ring α-hydroxy/thiol N-heterocycles → lactam/thiolactam (IUPAC 2013).

α-position enol tautomers of imidazole, pyrazole, and oxazole prefer
the N-H lactam/thiolactam form per IUPAC 2013 P-31.1.7.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-imidazol-2-ol → 1H-imidazol-2(3H)-one
    ("Oc1ncc[nH]1",             "1H-imidazol-2(3H)-one"),
    # 1H-pyrazol-3-ol → 1H-pyrazol-3(2H)-one
    ("Oc1cc[nH]n1",             "1H-pyrazol-3(2H)-one"),
    # 1,3-oxazol-2-ol → 1,3-oxazol-2(3H)-one
    ("Oc1ncco1",                "1,3-oxazol-2(3H)-one"),
    # thiol counterparts
    ("Sc1ncc[nH]1",             "1H-imidazol-2(3H)-thione"),
    ("Sc1cc[nH]n1",             "1H-pyrazol-3(2H)-thione"),
    ("Sc1ncco1",                "1,3-oxazol-2(3H)-thione"),
    # Regression: parent rings unaffected
    ("c1cnc[nH]1",              "1H-imidazole"),
    ("c1cn[nH]c1",              "1H-pyrazole"),
    ("c1cnco1",                 "1,3-oxazole"),
    ("c1cncs1",                 "1,3-thiazole"),
    # Regression: benzo-fused 5-membered ring unchanged (Phase 746)
    ("Oc1nc2ccccc2[nH]1",       "1H-benzimidazol-2(3H)-one"),
])
def test_phase754_five_membered_ring_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
