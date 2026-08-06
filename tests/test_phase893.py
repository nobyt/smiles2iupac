"""Phase 893: carboxylate anions dropped ALL chain substituents.

Found via the same fresh probe sweep as Phase 889-892 (originally noticed
as a multi-component salt bug: CC(N)C(=O)[O-].[NH4+], the ammonium salt of
alanine's conjugate base, was named "ammonium propanoate" -- losing the
2-amino group and colliding with plain ammonium propanoate).

Investigating showed this wasn't amino-specific or salt-specific at all:
_name_carboxylate (the namer for any deprotonated -C(=O)[O-]) never
collected substituents on the acid chain at ALL, unlike the neutral
carboxylic-acid path (which goes through the general chain-finder pipeline
and calls collect_substituents normally). ANY substituent -- amino,
hydroxy, halogen, anything -- silently vanished the moment the acid was
deprotonated:
  CC(N)C(=O)[O-]  -> "propanoate"   (drops 2-amino)
  CC(O)C(=O)[O-]  -> "propanoate"   (drops 2-hydroxy, same as unsubstituted!)
  CC(Cl)C(=O)[O-] -> "propanoate"   (drops 2-chloro)
All three different real anions collapsed to the exact same wrong name.

Fixed by calling the same collect_substituents helper used elsewhere in
this codebase (e.g. the amide namer) on the acid chain, excluding the
carboxylate oxygens, and prefixing the results the same way ("{loc}-{name}"
per substituent, joined by "-", concatenated directly onto the acid stem --
no separating hyphen before the stem, matching ordinary substituent-prefix-
to-stem concatenation). Verified via OPSIN / RDKit canonical-SMILES
round-trip.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CC(N)C(=O)[O-]",  "2-aminopropanoate"),
    ("CC(O)C(=O)[O-]",  "2-hydroxypropanoate"),
    ("CC(Cl)C(=O)[O-]", "2-chloropropanoate"),
    ("NCC(=O)[O-]",     "2-aminoacetate"),
    ("ClCC(=O)[O-]",    "2-chloroacetate"),
    # the original multi-component-salt bug report
    ("CC(N)C(=O)[O-].[NH4+]", "ammonium 2-aminopropanoate"),
    # regression: unsubstituted carboxylates unchanged
    ("CC(=O)[O-]",    "acetate"),
    ("C(=O)[O-]",     "formate"),
    ("CC=CC(=O)[O-]", "but-2-enoate"),
])
def test_phase893_carboxylate_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase893_distinct_anions_not_collapsed():
    amino = smiles_to_iupac("CC(N)C(=O)[O-]")
    hydroxy = smiles_to_iupac("CC(O)C(=O)[O-]")
    chloro = smiles_to_iupac("CC(Cl)C(=O)[O-]")
    plain = smiles_to_iupac("CCC(=O)[O-]")
    assert len({amino, hydroxy, chloro, plain}) == 4
