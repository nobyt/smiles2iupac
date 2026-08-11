"""Phase 906: the linear-chain substituent scanner in _name_carbon_substituent
(substituent.py) had no ketone ("oxo") detection at all, unlike its
halogen/hydroxy/amino/sulfanyl siblings sitting right next to it.

`CC(=O)OCC(=O)C` (2-oxopropyl acetate, i.e. acetoxyacetone) was named
plain "propyl acetate" -- the entire ketone on the O-alkyl substituent
chain silently vanished, colliding with the real, unrelated propyl
acetate. This was flagged as a new, separate, unrelated bug at the end of
Phase 905 (confirmed pre-existing via `git stash`, not introduced by any
of the 902-905 fixes) and is now closed.

Root cause: the chain-scanning block in `_name_carbon_substituent` that
builds names like "2-hydroxyethyl"/"2-chloroethyl"/"2-aminoethyl" checks
halogen, hydroxy, amino, and sulfanyl substituents on non-root chain
carbons in sequence -- but never checked for an exocyclic C=O (ketone).
A root atom's OWN C=O is handled separately by the earlier acyl branch
(formyl/acetyl/propanoyl), so this gap only affected ketones at
non-root chain positions. Since aromatic-ketone substituents
(`CC(=O)OCC(=O)c1ccccc1`) route through a completely different branched-
substituent code path that already had correct oxo handling, only the
*aliphatic* ketone-in-chain case was broken -- which is why it slipped
past the earlier probe sweeps in this session (they mostly used CF3,
not a chain ketone).

Fixed by adding an "oxo" detection block mirroring the existing
halogen/hydroxy/amino/sulfanyl blocks exactly (same locant-grouping /
multiplier logic). Since `_name_carbon_substituent` is the same shared
low-level primitive fixed for acyl groups in Phase 905, this single fix
propagates to every substituent-chain call site (ester O-alkyl groups,
N-substituents on amides/amines, etc.), not just esters.

OPSIN-verified.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CC(=O)OCC(=O)C", "2-oxopropyl acetate"),
    ("CC(=O)OCCC(=O)C", "3-oxobutyl acetate"),
    ("CC(C)C(=O)OCC(=O)C", "2-oxopropyl 2-methylpropanoate"),
    ("CC(=O)NCC(=O)C", "N-(2-oxopropyl)acetamide"),
    # regression: aromatic ketone (different code path) unaffected
    ("CC(=O)OCC(=O)c1ccccc1", "2-oxo-2-phenylethyl acetate"),
    # regression: other chain-substituent types unaffected
    ("CC(=O)OCCO", "2-hydroxyethyl acetate"),
    ("CC(=O)OCCCl", "2-chloroethyl acetate"),
    ("CC(=O)OCC", "ethyl acetate"),
    ("CC(=O)OCCC", "propyl acetate"),
    ("CC(=O)Nc1ccccc1", "N-phenylacetamide"),
])
def test_phase906_chain_ketone_substituent_not_dropped(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase906_oxopropyl_acetate_not_confused_with_propyl_acetate():
    assert (smiles_to_iupac("CC(=O)OCC(=O)C")
            != smiles_to_iupac("CC(=O)OCCC"))
