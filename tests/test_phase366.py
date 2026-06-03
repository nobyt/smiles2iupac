"""Phase 366: Chain substituents in N-substituted amide names.

_name_secondary_tertiary_amide previously omitted branch substituents on the
acid chain, yielding e.g. 'N-methylpropanamide' for CC(C)C(=O)NC instead of
the correct '2-methyl-N-methylpropanamide'.

Fixed by calling collect_substituents on the acid chain (excluding the N and
O atoms) and merging with N-substituents using alphabetical ordering with
numeric locants preceding letter (N-) locants for same-name groups.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Secondary amide with chain branch
    ("CC(C)C(=O)NC",    "2-methyl-N-methylpropanamide"),
    ("CC(C)C(=O)NCC",   "N-ethyl-2-methylpropanamide"),
    ("CCC(C)C(=O)NC",   "2-methyl-N-methylbutanamide"),
    # Tertiary amide with chain branch: same-name groups → numeric before N
    ("CC(C)C(=O)N(C)C", "2-methyl-N,N-dimethylpropanamide"),
    # Regressions: no chain branch → unchanged
    ("CC(=O)NC",        "N-methylacetamide"),
    ("CC(=O)N(C)C",     "N,N-dimethylacetamide"),
    ("CCCC(=O)NC",      "N-methylbutanamide"),
    ("CC(C)C(=O)N",     "2-methylpropanamide"),
])
def test_phase366_amide_chain_subs(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
