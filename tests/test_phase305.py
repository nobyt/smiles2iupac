"""Phase 305: imidoate ester naming fix (IUPAC 2013 P-65.1.2.3).

Ester of imidic acid RC(=NH)OH → "{alkyl} {stem}animidoate" (not "animidate").
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("C(=N)OC",     "methyl methanimidoate"),
    ("CC(=N)OC",    "methyl ethanimidoate"),
    ("CC(=N)OCC",   "ethyl ethanimidoate"),
    ("CCC(=N)OC",   "methyl propanimidoate"),
    # regressions: imidic acid unchanged
    ("CC(=N)O",     "ethanimidic acid"),
    # regressions: amine/amide unchanged
    ("CC(=O)N",     "acetamide"),
])
def test_phase305_imidoate_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
