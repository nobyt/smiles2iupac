"""Phase 908: a systematic sweep for the 902-907 chain-substituent-dropping
bug shape, this time driven by a script that greps every `_collect_acid_chain`
call site in group_namers.py and flags any surrounding function that never
calls collect_substituents (directly or via a known wrapper like
_pivot_chain_sub_prefix). This turned up 13 more standalone namers with the
identical gap, all confirmed real with a CF2/CF3-substituted probe molecule
and OPSIN-verified:

- _name_diester / _name_dicarboxylate: diester/dicarboxylate-dianion chains
  never collected substituents, e.g. dimethyl 2,2-difluoromalonate ->
  "dimethyl malonate" (drops both F).
- _name_carboxylate: a *different*, older bug (Phase893's own fix) -- it DID
  collect substituents but joined them with a naive `f"{loc}-{name}"` string
  join instead of `_build_prefix`'s multiplier-grouping, so repeated
  substituents didn't merge: trifluoroacetate -> "2-fluoro-2-fluoro-2-
  fluoroacetate" instead of "2,2,2-trifluoroacetate". This was flagged as a
  known latent bug in the Phase902 memory note and finally fixed here.
- _name_thioic_acid / _name_imidic_acid: same gap on the acid chain.
- _name_diacid_halide (aliphatic chain path; the aromatic-ring path was
  already fixed in Phase895): needed the halide atoms themselves excluded
  too, not just pgrp.atom_indices, or the acid chain's own halogens got
  double-counted as substituents.
- _name_peroxyacid / _name_acyl_azide: same gap on the acid chain.
- _name_thioamide / _name_selenoamide / _name_telluramide: same gap, PLUS
  the acid chain's own C=S/C=Se/C=Te got picked up as a spurious
  "sulfanyl"/"selanyl"/"tellanyl" substituent once collect_substituents was
  added, since these three built their own ad-hoc exclusion set instead of
  reusing pgrp.atom_indices like the sibling dithioamide/diselenoamide
  functions already did -- needed `_carbonyl_dbl_bonded` (generalized here
  to include Te, not just O/S/Se) added explicitly.
- _name_dithioamide / _name_diselenoamide: same gap on the acid chain
  between the two C=S/C=Se carbons.
- _name_acyl_peroxide: dropped substituents on BOTH acyl chains at once,
  e.g. bis(trifluoroacetyl) peroxide collapsed to plain "diethanoyl
  peroxide", colliding with the real unrelated diacetyl peroxide.

Two more bugs found separately by grepping for other naive same-function DFS
chain walkers (the same shape that made `_alkylidene_name` wrong before this
session's earlier Phase907 fix):

- _name_nitrone: only counted chain length via a manual DFS, e.g. a
  CF3-substituted nitrone collided with the plain unsubstituted one (both ->
  "N-methylethanimine N-oxide"). Rewritten to reuse find_principal_chain
  (anchored via the unrelated "aldehyde" spec purely for its anchor_c1=True
  chain-orientation behavior -- nitrone always fixes the imine carbon at
  locant 1) + collect_substituents.
- _name_azo_compound: only ever compared chain LENGTH between the two
  R-N=N-R' sides (via another manual DFS) and built the name from bare
  CHAIN_PREFIX stems -- so it only handled the symmetric case, and even then
  any substituent silently vanished: bis(trifluoromethyl)diazene (a real,
  well-known fluorine-chemistry reagent) collided with plain dimethyldiazene
  (both sides "1 carbon"). Rewritten to use the shared _name_carbon_substituent
  builder per side, which also incidentally adds previously-unsupported
  asymmetric R-N=N-R' naming (e.g. ethylmethyldiazene) for free.
- _name_diazo_compound / _name_diazonium: same chain-length-only DFS shape;
  _name_diazonium's aliphatic branch even had a code comment admitting it
  ("Pick longest chain (simplified: just follow first)"). 2-diazo-1,1,1-
  trifluoroethane collapsed to "diazoethane"; 2,2,2-trifluoroethyldiazonium
  collapsed to "ethanediazonium".

All OPSIN-verified.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("COC(=O)C(F)(F)C(=O)OC", "dimethyl 2,2-difluoromalonate"),
    ("CCOC(=O)CC(=O)OCC", "diethyl malonate"),
])
def test_phase908_diester_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("[O-]C(=O)C(F)(F)C(=O)[O-]", "2,2-difluoromalonate"),
    ("[O-]C(=O)CC(=O)[O-]", "malonate"),
    ("[O-]C(=O)CCCCC(=O)[O-]", "adipate"),
])
def test_phase908_dicarboxylate_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("[O-]C(=O)C(F)(F)F", "2,2,2-trifluoroacetate"),
    ("CC(N)C(=O)[O-]", "2-aminopropanoate"),
    ("[O-]C(=O)C", "acetate"),
    ("[O-]C(=O)CO", "2-hydroxyacetate"),
])
def test_phase908_carboxylate_repeated_substituent_multiplier(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("FC(F)(F)C(=O)S", "trifluoroethanethioic S-acid"),
    ("CC(=O)S", "ethanethioic S-acid"),
    ("CC(=S)O", "ethanethioic O-acid"),
    ("CC(=S)S", "ethanedithioic acid"),
])
def test_phase908_thioic_acid_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("FC(F)(F)C(=N)O", "trifluoroethanimidic acid"),
    ("CC(=N)O", "ethanimidic acid"),
    ("CC(=NC)O", "N-methylethanimidic acid"),
])
def test_phase908_imidic_acid_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("ClC(=O)C(F)(F)C(=O)Cl", "2,2-difluoropropanedioyl dichloride"),
    ("FC(=O)C(F)(F)C(=O)Cl", "2,2-difluoropropanedioyl chloride fluoride"),
    ("ClC(=O)CC(=O)Cl", "propanedioyl dichloride"),
    ("ClC(=O)C(=O)Cl", "ethanedioyl dichloride"),
    ("O=C(Cl)c1ccccc1C(=O)Cl", "benzene-1,2-dicarbonyl dichloride"),
])
def test_phase908_diacid_halide_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("FC(F)(F)C(=O)OO", "trifluoroethaneperoxoic acid"),
    ("CC(=O)OO", "ethaneperoxoic acid"),
    ("c1ccccc1C(=O)OO", "benzeneperoxoic acid"),
])
def test_phase908_peroxyacid_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("FC(F)(F)C(=O)N=[N+]=[N-]", "trifluoroacetyl azide"),
    ("CC(=O)N=[N+]=[N-]", "acetyl azide"),
    ("O=CN=[N+]=[N-]", "formyl azide"),
    ("CCCC(=O)N=[N+]=[N-]", "butanoyl azide"),
])
def test_phase908_acyl_azide_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("FC(F)(F)C(=S)N", "trifluoroethanethioamide"),
    ("CC(=S)N", "ethanethioamide"),
    ("CC(=S)NC", "N-methylethanethioamide"),
    ("c1ccccc1C(=S)N", "benzothioamide"),
    ("c1ccncc1C(=S)N", "pyridine-3-carbothioamide"),
])
def test_phase908_thioamide_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("FC(F)(F)C(=[Se])N", "trifluoroethaneselenoamide"),
    ("CC(=[Se])N", "ethaneselenoamide"),
])
def test_phase908_selenoamide_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("FC(F)(F)C(=[Te])N", "trifluoroethanetelluramide"),
    ("CC(=[Te])N", "ethanetelluramide"),
])
def test_phase908_telluramide_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("NC(=S)C(F)(F)C(=S)N", "2,2-difluoropropanedithioamide"),
    ("NC(=S)C(=S)N", "ethanedithioamide"),
    ("NC(=S)C1CCC(C(=S)N)CC1", "cyclohexane-1,4-dicarbothioamide"),
])
def test_phase908_dithioamide_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("NC(=[Se])C(F)(F)C(=[Se])N", "2,2-difluoropropanediselenoamide"),
    ("NC(=[Se])C(=[Se])N", "ethanediselenoamide"),
])
def test_phase908_diselenoamide_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("FC(F)(F)C(=O)OOC(=O)C(F)(F)F", "ditrifluoroethanoyl peroxide"),
    ("FC(F)(F)C(=O)OOC(=O)C", "ethanoyl trifluoroethanoyl peroxide"),
    ("CC(=O)OOC(=O)C", "diethanoyl peroxide"),
    ("CC(=O)OOC(=O)CC", "ethanoyl propanoyl peroxide"),
])
def test_phase908_acyl_peroxide_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("C[N+](=CC(F)(F)F)[O-]", "N-methyl-2,2,2-trifluoroethanimine N-oxide"),
    ("C[N+](=CC)[O-]", "N-methylethanimine N-oxide"),
    ("C[N+](=C)[O-]", "N-methylmethanimine N-oxide"),
])
def test_phase908_nitrone_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("FC(F)(F)N=NC(F)(F)F", "bis(trifluoromethyl)diazene"),
    ("CN=NC", "dimethyldiazene"),
    ("CCN=NCC", "diethyldiazene"),
    ("CN=NCC", "ethylmethyldiazene"),
    ("c1ccccc1N=Nc1ccccc1", "diphenyldiazene"),
])
def test_phase908_azo_compound_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("FC(F)(F)C=[N+]=[N-]", "1-diazo-2,2,2-trifluoroethane"),
    ("C=[N+]=[N-]", "diazomethane"),
    ("CC=[N+]=[N-]", "diazoethane"),
    ("CCC=[N+]=[N-]", "diazopropane"),
])
def test_phase908_diazo_compound_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("FC(F)(F)C[N+]#N", "2,2,2-trifluoroethanediazonium"),
    ("C[N+]#N", "methanediazonium"),
    ("CC[N+]#N", "ethanediazonium"),
    ("c1ccccc1[N+]#N", "benzenediazonium"),
])
def test_phase908_diazonium_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("Cc1ccc(cc1)[N+]#N", "4-methylbenzenediazonium"),
    ("FC(F)(F)c1ccccc1[N+]#N", "2-(trifluoromethyl)benzenediazonium"),
    ("Clc1ccccc1[N+]#N", "2-chlorobenzenediazonium"),
    ("c1ccccc1[N+]#N", "benzenediazonium"),
])
def test_phase908_arenediazonium_ring_substituent(smiles, expected):
    """Previously ANY ring substituent on an arenediazonium salt fell to a
    nonsensical generic fallback (e.g. "1-(N)-4-methylbenzene", dropping
    the diazonium group's charge entirely) since only a bare, fully
    unsubstituted benzene ring was matched. Fixed by reusing the shared
    `_aryl_sulfonyl_prefix` ring-substituent builder."""
    assert smiles_to_iupac(smiles) == expected
