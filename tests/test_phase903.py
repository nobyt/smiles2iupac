"""Phase 903: the rest of the sulfur/selenium/tellurium/nitroso chain-namer
family had the SAME chain-substituent-dropping bug fixed for the core
sulfonic-acid family in Phase 902.

After fixing `_name_sulfur_chain_acid`/`_name_sulfonate_anion`/etc in
Phase 902, a systematic check of every other function using the same
`_chain_through_pivot` chain-building helper (grep for the call site) found
~14 more namers with the identical bug -- each builds its carbon chain
independently and never runs it through `collect_substituents`:

  _name_sulfonate_sulfinate_ester, _name_sulfonyl_azide,
  _name_sulfonohydrazide, _name_sulfinylhydrazide, _name_sulfenamide,
  _name_sulfinamide, _name_sulfonimidamide, _name_sulfenic_acid,
  _name_sulfenyl_halide, _name_sulfenate_ester, _name_sulfinyl_chloride,
  _name_chalcogen_oxyacid (selenonic/seleninic/selenenic/telluronic/
  tellurinic/tellurenic acid -- 6 group types via one shared function),
  _name_nitroso, _name_di_xisocyanate, the sulfonylurea N-substituent
  helper inside _name_substituted_urea_if_match, and
  _name_o_substituted_oxime.

Extracted a shared `_pivot_chain_sub_prefix(graph, chain, excluded_atoms,
get_atom)` helper (mirrors assemble_name's own _build_prefix + the
chain-length-1 locant-omission rule) and applied it to all of them,
including refactoring the 6 functions fixed inline in Phase 902 to use the
same helper instead of duplicating the collect_substituents/_build_prefix
boilerplate. `_name_nitroso` needed different handling since "nitroso"
itself is a substituent-style prefix (not a suffix) -- it's now merged
into the same alphabetized substituent list rather than prepended
separately, so trifluoro + nitroso group correctly into
"trifluoronitrosomethane". Verified via OPSIN where OPSIN's own grammar
supports the construct (it doesn't for the O-substituted-oxime family --
a pre-existing OPSIN limitation unrelated to this fix, confirmed by
checking that OPSIN also rejects the untouched control name
"O-methylpropan-2-one oxime").
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # sulfinyl / sulfenyl chloride
    ("FC(F)(F)S(=O)Cl", "trifluoromethanesulfinyl chloride"),
    ("FC(F)(F)SCl", "trifluoromethanesulfenyl chloride"),
    # sulfinamide / sulfenamide / sulfonimidamide
    ("FC(F)(F)S(=O)N", "trifluoromethanesulfinamide"),
    ("FC(F)(F)SN", "trifluoromethanesulfenamide"),
    # selenonic/seleninic/selenenic + telluronic/tellurinic/tellurenic acid
    ("FC(F)(F)[Se](=O)(=O)O", "trifluoromethaneselenonic acid"),
    ("FC(F)(F)[Se](=O)O", "trifluoromethaneseleninic acid"),
    ("FC(F)(F)[Se]O", "trifluoromethaneselenenic acid"),
    ("FC(F)(F)[Te](=O)(=O)O", "trifluoromethanetelluronic acid"),
    ("FC(F)(F)[Te](=O)O", "trifluoromethanetellurinic acid"),
    ("FC(F)(F)[Te]O", "trifluoromethanetellurenic acid"),
    # sulfonohydrazide / sulfinylhydrazide / sulfonyl azide
    ("FC(F)(F)S(=O)(=O)NN", "trifluoromethanesulfonohydrazide"),
    ("FC(F)(F)S(=O)(=O)N=[N+]=[N-]", "trifluoromethanesulfonyl azide"),
    # sulfenate / sulfonate / sulfinate esters
    ("FC(F)(F)SOC", "methyl trifluoromethanesulfenate"),
    ("FC(F)(F)S(=O)OC", "methyl trifluoromethanesulfinate"),
    ("FC(F)(F)S(=O)(=O)OC", "methyl trifluoromethanesulfonate"),
    # nitroso: substituent merges into the same alphabetized prefix list
    ("FC(F)(F)N=O", "trifluoronitrosomethane"),
    # diisocyanate/diisothiocyanate chain substituents
    ("O=C=NC(Cl)CN=C=O", "1-chloro-1,2-diisocyanatoethane"),
    # sulfonylurea N-substituent (chain substituents on the sulfonyl side)
    ("FC(F)(F)S(=O)(=O)NC(=O)N", "N-trifluoromethanesulfonylurea"),
    # O-substituted oxime
    ("CC(=NOC)CCl", "O-methyl-3-chloropropan-2-one oxime"),
    ("ClCC=NOC", "O-methyl-2-chloroethanal oxime"),
    # regression: plain (unsubstituted) forms unchanged
    ("CS(=O)Cl", "methanesulfinyl chloride"),
    ("CSCl", "methanesulfenyl chloride"),
    ("CS(=O)N", "methanesulfinamide"),
    ("CSN", "methanesulfenamide"),
    ("C[Se](=O)(=O)O", "methaneselenonic acid"),
    ("C[Te](=O)O", "methanetellurinic acid"),
    ("CS(=O)(=O)NN", "methanesulfonohydrazide"),
    ("CS(=O)(=O)N=[N+]=[N-]", "methanesulfonyl azide"),
    ("CSOC", "methyl methanesulfenate"),
    ("CN=O", "nitrosomethane"),
    ("CCN=O", "nitrosoethane"),
    ("CCCN=O", "1-nitrosopropane"),
    ("O=C=NCCN=C=O", "1,2-diisocyanatoethane"),
    ("S=C=NCCN=C=S", "1,2-diisothiocyanatoethane"),
    ("CS(=O)(=O)NC(=O)N", "N-methanesulfonylurea"),
    ("CC(=NOC)C", "O-methylpropan-2-one oxime"),
    ("CC=NOC", "O-methylethanal oxime"),
])
def test_phase903_chain_pivot_namers_keep_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase903_triflyl_chloride_not_confused_with_thionyl_variant():
    assert (smiles_to_iupac("FC(F)(F)S(=O)Cl")
            != smiles_to_iupac("CS(=O)Cl"))


def test_phase903_trifluoronitrosomethane_not_confused_with_nitrosomethane():
    assert smiles_to_iupac("FC(F)(F)N=O") != smiles_to_iupac("CN=O")


def test_phase903_chloro_diisocyanatoethane_not_confused_with_plain():
    assert (smiles_to_iupac("O=C=NC(Cl)CN=C=O")
            != smiles_to_iupac("O=C=NCCN=C=O"))
