"""Phase 221: symmetric perhalogenated ethane → omit locants (IUPAC 2013 P-14.5).

When all 6 hydrogen positions of ethane are replaced by the same halogen,
only one arrangement is possible, so locants are unnecessary.
"""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # hexafluoroethane (6F on ethane, only one possible arrangement)
    ("FC(F)(F)C(F)(F)F",   "hexafluoroethane"),
    # hexachloroethane
    ("ClC(Cl)(Cl)C(Cl)(Cl)Cl", "hexachloroethane"),
    # regression: asymmetric halogenation still needs locants
    ("FC(F)(F)CCl",        "2-chloro-1,1,1-trifluoroethane"),
    # regression: single substituent on ethane still no locant
    ("CCF",                "fluoroethane"),
    # regression: 4F on ethane still needs locants (ambiguous: 1,1,1,2 vs 1,1,2,2)
    ("FC(F)(F)CF",         "1,1,1,2-tetrafluoroethane"),
])
def test_phase221_perhalogenated_ethane(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
