"""Phase 904: acid anhydride and amidine namers dropped chain substituents --
a fresh instance of the same bug class fixed for sulfur/Se/Te chain acids
in Phase 902/903, found by probing other standalone (non-chain-finder-
pipeline) namers with a trifluoromethyl-substituted test molecule.

Trifluoroacetic anhydride (TFAA, `FC(F)(F)C(=O)OC(=O)C(F)(F)F`), one of the
most common reagents in organic synthesis, was named plain "acetic
anhydride" -- all 6 fluorines across both acyl groups dropped, colliding
with the real, unrelated acetic anhydride. Root cause:
`_name_anhydride`'s `_acid_stem_name` helper picks the acid word purely
from chain length (`CHAIN_PREFIX`/`_ACID_RETAINED` lookup) and never
collected substituents on either acyl chain.

Trifluoroacetamidine (`FC(F)(F)C(=N)N`) had the same bug: was named plain
"ethanimidamide", colliding with the real acetamidine. Root cause:
`_name_substituted_amidine` (which unconditionally handles ALL amidines,
substituted or not, per `__init__.py`'s `_name_acyclic` dispatch) builds
its chain via `_collect_acid_chain` but only ever looked at N-substituents,
never chain C-substituents.

Both fixed using the `_pivot_chain_sub_prefix` helper introduced in
Phase 903 (anhydride needed its own locant-omission variant, matching the
existing behavior of hardcoded retained names like "trifluoroacetic acid":
omit locants when chain length <=2 and only one substituent NAME is
present, since it's unambiguous). OPSIN-verified.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # acid anhydride: substituents on either/both acyl chains
    ("FC(F)(F)C(=O)OC(=O)C(F)(F)F", "trifluoroacetic anhydride"),
    ("ClCC(=O)OC(=O)CCl", "chloroacetic anhydride"),
    ("FC(F)(F)C(=O)OC(=O)C", "acetic trifluoroacetic anhydride"),
    # amidine: substituents on the acid-side chain
    ("FC(F)(F)C(=N)N", "2,2,2-trifluoroethanimidamide"),
    ("ClCC(=N)N", "2-chloroethanimidamide"),
    # regression: plain (unsubstituted) forms unchanged
    ("CC(=O)OC(=O)C", "acetic anhydride"),
    ("CCC(=O)OC(=O)CC", "propanoic anhydride"),
    ("CC(=O)OC(=O)CC", "acetic propanoic anhydride"),
    ("CC(=N)N", "ethanimidamide"),
    ("CC(=N)NC", "N-methylethanimidamide"),
    ("CC(=NC)N", "N'-methylethanimidamide"),
    ("CNC(=N)N", "N-methylguanidine"),
    ("CCCC(=N)N", "butanimidamide"),
])
def test_phase904_anhydride_amidine_chain_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase904_tfaa_not_confused_with_acetic_anhydride():
    assert (smiles_to_iupac("FC(F)(F)C(=O)OC(=O)C(F)(F)F")
            != smiles_to_iupac("CC(=O)OC(=O)C"))


def test_phase904_trifluoroacetamidine_not_confused_with_acetamidine():
    assert smiles_to_iupac("FC(F)(F)C(=N)N") != smiles_to_iupac("CC(=N)N")
