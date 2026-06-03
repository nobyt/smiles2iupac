"""Phase 373: Imidate ester suffix corrected to -imidate (IUPAC 2013 P-65.1.1.2).

The acid RC(=NH)OH is 'ethanimidic acid' (suffix -imidic acid).
Esters are named by replacing '-ic acid' with '-ate':
  ethanimidic acid → ethanimidate (not ethanimidate).

Previously the code produced 'ethanimidate' which is inconsistent with
the acid name 'ethanimidic acid' that was already correctly implemented.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Basic imidates
    ("C(=N)OC",          "methyl methanimidate"),
    ("CC(=N)OC",         "methyl ethanimidate"),
    ("CC(=N)OCC",        "ethyl ethanimidate"),
    ("CCC(=N)OC",        "methyl propanimidate"),
    ("CCC(=N)OCC",       "ethyl propanimidate"),
    # N-substituted imidates
    ("COC(=NC)C",        "methyl N-methylethanimidate"),
    ("CCOC(=NC)CC",      "ethyl N-methylpropanimidate"),
    # E/Z imidates
    ("C/C=C/C(=N)OC",    "methyl (2E)-but-2-enimidate"),
    # Regressions: acid form unchanged
    ("CC(=N)O",          "ethanimidic acid"),
    ("CCC(=N)O",         "propanimidic acid"),
    # Regressions: amide unchanged
    ("CC(=O)N",          "acetamide"),
])
def test_phase373_imidate_ester_suffix(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
