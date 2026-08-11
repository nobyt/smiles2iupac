"""Phase 905: imidate ester, thioester family, hydrazide, and N-acyl-amide
(imide) namers -- plus the generic acyl-substituent namer used EVERYWHERE
(N-acyl, O-acyl, ring-acyl substituents) -- all had the same
chain-substituent-dropping bug as Phase 902-904, found by continuing the
same CF3-substituted-input probe into ester/hydrazide/imide territory.

- `_name_imidate_ester`, `_name_thioester`, `_name_o_thioester`,
  `_name_s_dithioate_ester`, `_name_hydrazide`: each built its acid chain
  independently and never called collect_substituents, so e.g.
  `FC(F)(F)C(=O)NN` (trifluoroacetohydrazide) collapsed to plain
  "acetohydrazide". Fixed with `_pivot_chain_sub_prefix` (from Phase 903),
  now also excluding the carbonyl's own =O/=S via a new
  `_carbonyl_dbl_bonded` helper -- the first attempt at this fix
  regressed `CSC(=O)C` (S-methyl ethanethioate) to
  "S-methyl 1-oxoethanethioate" because the carbonyl oxygen, not excluded,
  got picked up as a fake "oxo" substituent; caught by testing the
  unsubstituted case immediately after the substituted one.

- `_name_carbon_substituent`'s acyl branch (substituent.py) -- the
  function that names ANY acyl group used as a substituent anywhere in
  the codebase (N-acyl amides/imides, O-acyl esters, acyl groups on
  rings) -- had the identical gap: picked "acetyl"/"propanoyl" purely
  from carbon count, ignoring substituents on those carbons. This made
  `CC(=O)NC(=O)C(F)(F)F` (N-(trifluoroacetyl)acetamide) collapse to
  "N-acetylacetamide", silently losing an entire trifluoroacetyl group.
  This is probably the highest-leverage single fix in this bug-class
  sweep since it's shared by every acyl-as-substituent call site.

- `_name_secondary_tertiary_amide`'s own separate chain-substituent
  collection (added correctly, unlike the others above) used a naive
  string-join that didn't group repeated identical substituents with a
  multiplier, so `FC(F)(F)C(=O)NC(=O)C` produced the malformed
  "N-acetyl-2-fluoro-2-fluoro-2-fluoroacetamide" instead of
  "N-acetyl-2,2,2-trifluoroacetamide" -- confirms the latent
  naive-join bug flagged (but not touched) in the Phase 902 memory note
  about `_name_carboxylate` can recur anywhere the same shortcut is used.

OPSIN-verified.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # imidate ester
    ("FC(F)(F)C(=N)OC", "methyl 2,2,2-trifluoroethanimidate"),
    # thioester family
    ("FC(F)(F)C(=O)SC", "S-methyl 2,2,2-trifluoroethanethioate"),
    ("FC(F)(F)C(=S)OC", "O-methyl 2,2,2-trifluoroethanethioate"),
    ("FC(F)(F)C(=S)SC", "S-methyl 2,2,2-trifluoroethanedithioate"),
    # hydrazide
    ("FC(F)(F)C(=O)NN", "trifluoroacetohydrazide"),
    ("ClCC(=O)NN", "chloroacetohydrazide"),
    # N-acyl amide / imide: substituent on either or both acyl chains
    ("CC(=O)NC(=O)C(F)(F)F", "N-(trifluoroacetyl)acetamide"),
    ("FC(F)(F)C(=O)NC(=O)C", "N-acetyl-2,2,2-trifluoroacetamide"),
    ("FC(F)(F)C(=O)NC(=O)C(F)(F)F", "2,2,2-trifluoro-N-(trifluoroacetyl)acetamide"),
    # acyl-as-substituent used elsewhere (ester O-alkyl side, N-substituent)
    ("FC(F)(F)C(=O)Oc1ccccc1", "phenyl 2,2,2-trifluoroacetate"),
    ("CNC(=O)C(F)(F)F", "2,2,2-trifluoro-N-methylacetamide"),
    # regression: plain (unsubstituted) forms unchanged
    ("CC(=N)OCC", "ethyl ethanimidate"),
    ("CSC(=O)C", "S-methyl ethanethioate"),
    ("COC(=S)C", "O-methyl ethanethioate"),
    ("CSC(=S)C", "S-methyl ethanedithioate"),
    ("CC(=O)NN", "acetohydrazide"),
    ("CC(=O)NNC", "N'-methylacetohydrazide"),
    ("CC(=O)NC(=O)C", "N-acetylacetamide"),
    ("CC(=O)Oc1ccccc1", "phenyl acetate"),
    ("CNC(=O)C", "N-methylacetamide"),
])
def test_phase905_acid_ester_hydrazide_imide_chain_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase905_trifluoroacetohydrazide_not_confused_with_acetohydrazide():
    assert (smiles_to_iupac("FC(F)(F)C(=O)NN")
            != smiles_to_iupac("CC(=O)NN"))


def test_phase905_n_trifluoroacetylacetamide_not_confused_with_n_acetylacetamide():
    assert (smiles_to_iupac("CC(=O)NC(=O)C(F)(F)F")
            != smiles_to_iupac("CC(=O)NC(=O)C"))
