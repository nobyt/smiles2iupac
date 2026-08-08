"""Phase 897: sulfenamides (R-S-NR'2) -- found in the same fresh probe
sweep as Phase 896's chalcogenolate anions.

CSNC (N-methyl-S-methyl sulfenamide, CH3-S-NH-CH3) was named
"sulfanylmethane" -- effectively "methanethiol", with the N-methyl group on
the sulfenamide nitrogen dropped entirely. There was no sulfenamide
functional-group detection at all (unlike the closely-related sulfenic
acid C-S-OH and sulfenate ester C-S-O-R, Phase 166/375, which already had
detectors); the C-S-N pattern fell straight through to a generic thioether
name.

Added detection mirroring the existing sulfenic_acid/sulfenate_ester elif
branches (S bonded to exactly 1 C and at least 1 N, no O/halogen), and a
dedicated namer (_name_sulfenamide, following the same "{stem}anesulfen-
amide" + N-/N,N- substituent-prefix pattern already used for sulfinamide/
sulfonamide) registered directly in PGRP_DISPATCH -- avoiding the Phase 895
chain_template pitfall entirely, since dual chain+N-substituent naming
needs a dedicated function regardless. Verified via OPSIN.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("CSNC",   "N-methylmethanesulfenamide"),
    ("CSN",    "methanesulfenamide"),
    ("CCSN",   "ethanesulfenamide"),
    ("CSNCC",  "N-ethylmethanesulfenamide"),
    # regression: sulfenic acid / sulfenate ester / sulfenyl halide unchanged
    ("CSO",    "methanesulfenic acid"),
    ("CSOC",   "methyl methanesulfenate"),
    ("CSCl",   "methanesulfenyl chloride"),
    # regression: plain thioether unchanged
    ("CSC",    "(methylsulfanyl)methane"),
])
def test_phase897_sulfenamide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase897_not_confused_with_methanethiol():
    result = smiles_to_iupac("CSNC")
    assert "sulfanylmethane" not in result
