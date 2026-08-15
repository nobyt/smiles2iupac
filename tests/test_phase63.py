"""Phase 63: ジアミン / トリアミン命名

Phase 925 update: the "トリアミン" case below (NCCCN(CCCN)CCC, a
tris-(3-aminopropyl)amine-type molecule with a TERTIARY central N) was
originally expected to produce "propane-1,3-triamine" -- confirmed via
OPSIN this was ALREADY an invalid, unparseable name predating this
session ("Mismatch between locant and multiplier counts (2 and 3)"), a
latent pre-existing bug from whenever diamine/triamine merging was first
implemented. Phase925 fixed the root cause (aggregate_groups() no longer
merges amine instances into diamine/triamine when any instance is
secondary/tertiary, since their extra N-alkyl substituents can't be
represented as plain merged-suffix locants) -- this correctly routes the
molecule through the existing single-amine + N-substituent-prefix
pathway instead, giving a valid, OPSIN+RDKit round-trip verified name.
"""

import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # ジアミン (diamine) — 直鎖
    ("NCCN", "ethane-1,2-diamine"),
    ("NCCCN", "propane-1,3-diamine"),
    ("NCCCCN", "butane-1,4-diamine"),
    ("NCCCCCN", "pentane-1,5-diamine"),
    ("NCCCCCCN", "hexane-1,6-diamine"),
    # ジアミン — 分岐鎖
    ("NCC(N)C", "propane-1,2-diamine"),
    ("NCC(N)CC", "butane-1,2-diamine"),
    ("NCC(N)CCC", "pentane-1,2-diamine"),
    # トリアミン (Phase925: 三級 N を含むため diamine/triamine にマージ
    # されず、単独アミン + N-置換基接頭辞として正しく命名される)
    ("NCCCN(CCCN)CCC", "3-[(3-aminopropyl)(propyl)amino]propan-1-amine"),
])
def test_phase63_diamine(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
