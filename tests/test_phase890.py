"""Phase 890: phosphoramidate esters (RO)_n-P(=O)(NH2)(OH)_{2-n}.

Found via the same fresh probe sweep as Phase 889. COP(=O)(N)OC (dimethyl
phosphoramidate -- relevant to ProTide-style prodrug chemistry, e.g. the
Sofosbuvir phosphoramidate class) was named "dimethyl phosphate", silently
dropping the amino group entirely.

Root cause: the phosphate_ester detection branch (`len(o_double) == 1 and
len(c_neighbors) == 0`) never checked for an N directly on phosphorus --
it only ever considered O-ester/O-H combinations. The phosphoramidic_acid
family (Phase 868) already handled the free-acid form (P-NH2 + P-OH only,
no P-O-R esters), but the ester form (P-NH2 + at least one P-O-R) fell
through past all of those acid-specific branches (which require the O's to
be exactly the right OH count) and into the generic phosphate_ester branch,
which has no N-awareness at all.

Fixed with a new phosphoramidate_ester group type/branch (inserted right
after phosphonamidic_acid, before the fallthrough to phosphate_ester) and
a namer mirroring the existing phosphate_ester ester-listing logic
(alphabetized alkyl esters, "hydrogen" for any remaining free -OH) plus an
N-/N,N- prefix for substituents on the amide nitrogen. Verified via OPSIN
parse-back / RDKit canonical-SMILES round-trip.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("COP(=O)(N)OC",   "dimethyl phosphoramidate"),
    ("COP(=O)(NC)OC",  "dimethyl N-methylphosphoramidate"),
    ("COP(=O)(N)O",    "methyl hydrogen phosphoramidate"),
    # regression: plain phosphate esters (no N) unchanged
    ("COP(=O)(O)OC",   "dimethyl hydrogen phosphate"),
    ("COP(=O)(OC)OC",  "trimethyl phosphate"),
    # regression: phosphoramidic acid (free acid, Phase 868) unchanged
    ("CNP(=O)(O)O",    "N-methylphosphoramidic acid"),
])
def test_phase890_phosphoramidate_ester(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase890_not_confused_with_plain_phosphate():
    assert smiles_to_iupac("COP(=O)(N)OC") != smiles_to_iupac("COP(=O)(OC)OC")
