"""Phase 902: sulfonic/sulfonate/sulfonamide/sulfonyl-halide/disulfonic-acid
namers silently dropped every chain substituent (halo/hydroxy/...).

Found via a fresh probe sweep: trifluoromethanesulfonic acid (triflic acid,
one of the strongest known Brønsted acids and a ubiquitous reagent as its
anhydride/esters) was named "methanesulfonic acid" -- all 3 fluorines
dropped, colliding with the real, unrelated methanesulfonic acid:
  FC(F)(F)S(=O)(=O)O  -> wrongly "methanesulfonic acid"
  CS(=O)(=O)O          -> "methanesulfonic acid"      (same output!)

Root cause: `_name_sulfur_chain_acid` (backing sulfonic_acid/sulfinic_acid/
the Phase865 thio variants), `_name_sulfonate_anion`, `_name_sulfonamide`'s
acid-chain branch, `_name_sulfonyl_chloride`, `_name_disulfonic_acid` and
`_name_disulfonamide` were all written as standalone functions that build
the chain via `_chain_through_pivot` but never call `collect_substituents`
on it -- the same "namer bypasses the general chain+substituent pipeline"
bug class as Phase 893's carboxylate-anion fix. Fixed by adding a
collect_substituents + _build_prefix call (mirroring the main pipeline's
own `assemble_name`, including its chain-length-1 locant-omission rule) to
each of the 6 functions above. OPSIN-verified for several of these
(trifluoromethanesulfonic acid, 1,1-difluoroethanesulfonic acid,
2-hydroxyethanesulfonic acid, 1-chloroethanesulfonic acid,
trifluoromethanesulfonate).
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # sulfonic acid: chain substituents no longer dropped
    ("FC(F)(F)S(=O)(=O)O", "trifluoromethanesulfonic acid"),
    ("CC(F)(F)S(=O)(=O)O", "1,1-difluoroethanesulfonic acid"),
    ("ClCS(=O)(=O)O", "chloromethanesulfonic acid"),
    ("ClCCS(=O)(=O)O", "2-chloroethanesulfonic acid"),
    ("CC(Cl)S(=O)(=O)O", "1-chloroethanesulfonic acid"),
    ("OCCS(=O)(=O)O", "2-hydroxyethanesulfonic acid"),
    # sulfonate anion
    ("FC(F)(F)S(=O)(=O)[O-]", "trifluoromethanesulfonate"),
    # sulfonamide (acid-chain side substituents)
    ("FC(F)(F)S(=O)(=O)N", "trifluoromethanesulfonamide"),
    ("FC(F)(F)S(=O)(=O)NC", "N-methyltrifluoromethanesulfonamide"),
    # sulfonyl chloride
    ("FC(F)(F)S(=O)(=O)Cl", "trifluoromethanesulfonyl chloride"),
    ("CC(Cl)S(=O)(=O)Cl", "1-chloroethanesulfonyl chloride"),
    # disulfonic acid / disulfonamide
    ("OS(=O)(=O)C(Cl)CS(=O)(=O)O", "1-chloroethane-1,2-disulfonic acid"),
    ("NS(=O)(=O)C(Cl)CS(=O)(=O)N", "1-chloroethane-1,2-disulfonamide"),
    # regression: plain (unsubstituted) forms unchanged
    ("CS(=O)(=O)O", "methanesulfonic acid"),
    ("CS(=O)(=O)[O-]", "methanesulfonate"),
    ("CS(=O)(=O)N", "methanesulfonamide"),
    ("CS(=O)(=O)Cl", "methanesulfonyl chloride"),
    ("CCCS(=O)(=O)O", "propane-1-sulfonic acid"),
    ("OS(=O)(=O)CCS(=O)(=O)O", "ethane-1,2-disulfonic acid"),
])
def test_phase902_sulfur_acid_chain_substituents(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase902_triflic_acid_not_confused_with_methanesulfonic_acid():
    assert (smiles_to_iupac("FC(F)(F)S(=O)(=O)O")
            != smiles_to_iupac("CS(=O)(=O)O"))


def test_phase902_triflyl_chloride_not_confused_with_mesyl_chloride():
    assert (smiles_to_iupac("FC(F)(F)S(=O)(=O)Cl")
            != smiles_to_iupac("CS(=O)(=O)Cl"))
