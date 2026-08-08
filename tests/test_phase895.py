"""Phase 895: two bugs found in a fresh probe sweep after Phase 894.

1. Aromatic diacid halides (phthaloyl chloride and isomers) --
   O=C(Cl)c1ccccc1C(=O)Cl was named "methanedioyl dichloride", a
   1-carbon compound name that completely dropped the benzene ring.
   _name_diacid_halide only ever handled the linear-chain case (oxalyl
   chloride, malonyl chloride, ...) via _collect_acid_chain, which can't
   traverse out through an aromatic ring, so it collapsed to a
   1-carbon "chain". Fixed by mirroring the existing ring-aware logic
   already used for the analogous diCARBOXYLIC ACID case (_name_dioic_acid,
   Phase 311/523): detect when both acyl-halide carbons are substituents
   on the SAME benzene ring and emit "benzene-X,Y-dicarbonyl {halide}".
   Verified via OPSIN parse-back.

2. Alkoxide anions (the anion half of metal alkoxides like sodium
   ethoxide) -- CC[O-] was named "oxyethane", a nonsensical string with
   no IUPAC meaning. There was no "alkoxide" functional group detection
   at all; the deprotonated oxygen fell through to a generic substituent-
   naming path. Added detection (mirrors the existing "alcohol" detector,
   but for O with formal charge -1 and no H) and a FunctionalGroupSpec
   with suffix "olate". Discovered along the way that FunctionalGroupSpec's
   chain_template/chain_template_mb fields are actually DEAD CODE -- the
   real suffix-to-name logic is a giant per-suffix `if suffix == "X":`
   dispatch in name_assembler.py, and any suffix without its own explicit
   branch there falls through to a raw "{stem}{suffix}" concatenation with
   no locant or "ane"-linker handling at all (this is exactly what
   produced "etholate" instead of "ethanolate" on the first attempt). Added
   an explicit "olate" branch mirroring the existing "ol" branch's locant-
   omission rules. Verified via OPSIN, including sodium ethoxide itself.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # aromatic diacid halides
    ("O=C(Cl)c1ccccc1C(=O)Cl",   "benzene-1,2-dicarbonyl dichloride"),
    ("O=C(Cl)c1ccc(cc1)C(=O)Cl", "benzene-1,4-dicarbonyl dichloride"),
    ("O=C(Cl)c1cccc(c1)C(=O)Cl", "benzene-1,3-dicarbonyl dichloride"),
    # regression: linear-chain diacid halides (Phase pre-existing) unchanged
    ("ClC(=O)C(=O)Cl",  "ethanedioyl dichloride"),
    ("ClC(=O)CC(=O)Cl", "propanedioyl dichloride"),
    ("O=C(Cl)c1ccccc1", "benzoyl chloride"),
    # alkoxide anions
    ("CC[O-]",         "ethanolate"),
    ("C[O-]",          "methanolate"),
    ("CC(C)(C)[O-]",   "2-methylpropan-2-olate"),
    ("CCCC[O-]",       "butan-1-olate"),
    ("c1ccccc1[O-]",   "phenolate"),
    ("CC[O-].[Na+]",   "sodium ethanolate"),
    ("C[O-].[K+]",     "potassium methanolate"),
    # regression: neutral alcohols and carboxylates unchanged
    ("CCO",          "ethanol"),
    ("CC(=O)[O-]",   "acetate"),
])
def test_phase895_diacid_halide_and_alkoxide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase895_phthaloyl_chloride_not_dropped_to_1_carbon():
    result = smiles_to_iupac("O=C(Cl)c1ccccc1C(=O)Cl")
    assert "methane" not in result


def test_phase895_alkoxide_not_nonsensical():
    result = smiles_to_iupac("CC[O-]")
    assert "oxyethane" not in result
