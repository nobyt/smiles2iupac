"""Phase 320: cycloalkane diester, diselenoamide, dicarboximidamide (IUPAC 2013 P-65.1).

Extends Phase 268/308/315 handlers to support non-aromatic (cycloalkane) rings,
and adds benzene/cycloalkane dicarboximidamide (diamidine) naming.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # diester cycloalkane
    ("COC(=O)C1CCC(C(=O)OC)CC1",   "dimethyl cyclohexane-1,4-dicarboxylate"),
    ("CCOC(=O)C1CCC(C(=O)OCC)CC1", "diethyl cyclohexane-1,4-dicarboxylate"),
    ("COC(=O)C1CCCC1C(=O)OC",      "dimethyl cyclopentane-1,2-dicarboxylate"),
    # diselenoamide cycloalkane
    ("NC(=[Se])C1CCC(C(=[Se])N)CC1", "cyclohexane-1,4-dicarboselenoamide"),
    ("NC(=[Se])C1CCCC1C(=[Se])N",    "cyclopentane-1,2-dicarboselenoamide"),
    # dicarboximidamide on ring
    ("NC(=N)C1CCC(C(=N)N)CC1",     "cyclohexane-1,4-dicarboximidamide"),
    ("NC(=N)c1ccc(C(=N)N)cc1",     "benzene-1,4-dicarboximidamide"),
    # regressions: chain diester unchanged
    ("CCOC(=O)CC(=O)OCC",          "diethyl malonate"),
    # regressions: benzene diester unchanged
    ("COC(=O)c1ccc(C(=O)OC)cc1",   "dimethyl benzene-1,4-dicarboxylate"),
    # regressions: chain diselenoamide unchanged
    ("NC(=[Se])C(=[Se])N",         "ethanediselenoamide"),
    # regressions: single amidine unchanged
    ("NC(=N)C",                    "ethanimidamide"),
    ("NC(=N)c1ccccc1",             "benzenecarboximidamide"),
])
def test_phase320_cycloalkane_digroups(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
