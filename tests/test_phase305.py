"""Phase 305: imidate ester naming fix (IUPAC 2013 P-65.1.2.3).

Ester of imidic acid RC(=NH)OH → "{alkyl} {stem}animidate" (not "animidate").
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C(=N)OC",     "methyl methanimidate"),
    ("CC(=N)OC",    "methyl ethanimidate"),
    ("CC(=N)OCC",   "ethyl ethanimidate"),
    ("CCC(=N)OC",   "methyl propanimidate"),
    # regressions: imidic acid unchanged
    ("CC(=N)O",     "ethanimidic acid"),
    # regressions: amine/amide unchanged
    ("CC(=O)N",     "acetamide"),
])
def test_phase305_imidate_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
