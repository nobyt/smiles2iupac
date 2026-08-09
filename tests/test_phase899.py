"""Phase 899: azanide / phosphanide anions -- the N/P mirror of Phase 895/896's
alkoxide/chalcogenolate fixes, and the anionic counterparts of the existing
ammonium/phosphanium cations.

CC[NH-] (the conjugate base of ethylamine) was named "aminoethane" -- the
negative charge silently dropped, reporting the anion as identical to the
neutral amine. CC[P-]C was named "ethylmethylphosphane" -- same bug for the
phosphorus analog.

Root cause: same shape as every anion gap fixed this session -- no
"azanide"/"phosphanide" functional-group detection existed, so the
deprotonated N/P fell through to a generic path (the amine/phosphane
branches require an exact H/C count that a charged, one-fewer-bond anion
doesn't match, and nothing else caught the case).

Fixed by adding detection guarded on formal_charge == -1 (placed before the
primary/secondary/tertiary amine branches for N; a single new elif before
the phosphane branch for P), and namers reusing the EXISTING
_name_by_c_substituents helper with suffix "azanide"/"phosphanide" -- these
follow the same substituent-prefix-on-element-parent pattern as
ammonium/phosphanium (Phase 146/518), not the carbon-chain-suffix pattern,
so no dead chain_template pitfall to worry about here. Verified via OPSIN.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CC[NH-]",  "ethylazanide"),
    ("C[NH-]",   "methylazanide"),
    ("CC[P-]C",  "ethylmethylphosphanide"),
    ("CC[PH-]",  "ethylphosphanide"),
    # regression: neutral amine/phosphane unchanged
    ("CCN",   "ethanamine"),
    ("CCPC",  "ethylmethylphosphane"),
    # regression: ammonium/phosphanium cations unchanged
    ("[NH4+]",        "ammonium"),
    ("C[P+](C)(C)C",  "tetramethylphosphanium"),
])
def test_phase899_azanide_phosphanide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase899_azanide_not_confused_with_neutral_amine():
    result = smiles_to_iupac("CC[NH-]")
    assert result != smiles_to_iupac("CCN")
