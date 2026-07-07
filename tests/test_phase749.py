"""Phase 749: α-thiol N-heterocycles → preferred thiolactam tautomers (IUPAC 2013).

Parallel to Phases 744–748 for OH→lactam: SH at α-positions adjacent to ring N
prefers the N-H thiolactam (thione) form (IUPAC 2013 P-31.1.7).
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 6-membered monocyclic thiol → thione
    ("Sc1ccccn1",               "pyridin-2(1H)-thione"),
    ("Sc1ccncc1",               "pyridin-4(1H)-thione"),
    ("Sc1ncccn1",               "pyrimidin-2(1H)-thione"),
    ("Sc1ccncn1",               "1H-pyrimidin-6-thione"),
    ("Sc1cnccn1",               "pyrazin-2(1H)-thione"),
    ("Sc1cccnn1",               "1H-pyridazin-6-thione"),
    # benzo-fused bicyclic thiol → thione
    ("Sc1nc2ccccc2cc1",         "quinolin-2(1H)-thione"),
    ("Sc1ccnc2ccccc12",         "quinolin-4(1H)-thione"),
    ("Sc1cnc2ccccc2n1",         "quinoxalin-2(1H)-thione"),
    # 5-membered fused thiol → thione
    ("Sc1nc2ccccc2s1",          "1,3-benzothiazol-2(3H)-thione"),
    ("Sc1nc2ccccc2o1",          "1,3-benzoxazol-2(3H)-thione"),
    ("Sc1nc2ccccc2[nH]1",       "1H-benzimidazol-2(3H)-thione"),
    # Regression: keto (thione) SMILES unchanged
    ("S=C1NC=CC=C1",            "pyridin-2(1H)-thione"),
    ("S=C1NC=CC=N1",            "pyrimidin-2(1H)-thione"),
    ("S=C1NC=NC=C1",            "1H-pyrimidin-6-thione"),
    ("S=C1NC=CN=C1",            "pyrazin-2(1H)-thione"),
    # Regression: parent rings unaffected
    ("c1ccccn1",                "pyridine"),
    ("c1ccncn1",                "pyrimidine"),
    ("c1cnccn1",                "pyrazine"),
])
def test_phase749_azine_thiol_tautomers(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
