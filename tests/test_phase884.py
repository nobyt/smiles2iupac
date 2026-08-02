"""Phase 884: P=N/S=N/As=N imine family (iminophosphorane, sulfilimine,
sulfoximine, arsane imine).

These four groups were flagged (but not implemented) during the Phase 860
session as genuinely new functional groups needing real imine naming logic,
not a mechanical mirror-fix. Before this phase all four silently dropped
the =N entirely:

  C[P](C)(C)=N     -> "trimethylphosphane"   (drops =N)
  C[S](C)=N        -> "(methylsulfanyl)methane" (drops =N, wrong structure --
                       falls to plain thioether naming)
  C[S](C)(=O)=N    -> "dimethyl sulfoxide"   (drops =N, collides with the
                       real dimethyl sulfoxide CS(=O)C)
  C[As](C)(C)=N    -> "trimethylarsane"      (drops =N)

Researched and verified via OPSIN parse-back (bundled
opsin-cli-2.9.0-jar-with-dependencies.jar):
  - Iminophosphorane/arsane imine take the "-imine" SUFFIX on the
    phosphane/arsane parent (phosphan-imine -> phosphanimine, elided;
    arsan-imine -> arsanimine), mirroring how amine on B/Si is a suffix
    (Phase 883) and how organic imines already work in this codebase
    (ethanimine, N-methylethanimine). NOT a "amino"-style substituent
    prefix. Confirmed "trimethylphosphanimine" -> CP(=N)(C)C,
    "N-methyltrimethylphosphanimine" -> CN=P(C)(C)C (both round-trip via
    RDKit canonical SMILES equality).
  - Sulfilimine/sulfoximine use the SAME functional-class two-word format
    as the existing sulfoxide/sulfone namer ("dimethyl sulfoxide" style),
    just with the new type word: "dimethyl sulfilimine" -> CS(=N)C,
    "dimethyl sulfoximine" -> CS(=O)(=N)C. N-substitution gets an "N-name"
    prefix in front (no extra hyphen needed, OPSIN parses
    "N-methyldimethyl sulfilimine" -> CN=S(C)C correctly).

Implementation: new phosphine_imine/arsane_imine functional-group branches
mirror the existing phosphine_oxide/sulfide and arsane_oxide/sulfide
detection exactly (same C/H-only-else guard), just swapped to a
double-bonded N; new sulfilimine/sulfoximine branches inserted before
sulfoxide/sulfone in _detect_sulfur_groups (same priority position, since
sulfoxide's condition would otherwise wrongly match sulfoximine first).
Namers: _name_phosphine_arsane_imine (shared P/As helper) and
_name_sulfilimine_sulfoximine (reuses _dual_c_group_prefix, extracted from
the existing _name_sulfoxide_sulfone).
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # iminophosphorane
    ("C[P](C)(C)=N",  "trimethylphosphanimine"),
    ("C[P](C)(C)=NC", "N-methyltrimethylphosphanimine"),
    # sulfilimine
    ("C[S](C)=N",     "dimethyl sulfilimine"),
    ("C[S](C)=NC",    "N-methyldimethyl sulfilimine"),
    # sulfoximine
    ("C[S](C)(=O)=N",  "dimethyl sulfoximine"),
    ("C[S](C)(=O)=NC", "N-methyldimethyl sulfoximine"),
    # arsane imine
    ("C[As](C)(C)=N",  "trimethylarsanimine"),
    ("C[As](C)(C)=NC", "N-methyltrimethylarsanimine"),
    # regression: the O/S analogs (phosphine/arsane oxide/sulfide, Phase
    # 857-860) and sulfoxide/sulfone (Phase 519) unchanged
    ("C[P](C)(C)=O",  "trimethylphosphane oxide"),
    ("C[P](C)(C)=S",  "trimethylphosphane sulfide"),
    ("C[As](C)(C)=O", "trimethylarsane oxide"),
    ("C[As](C)(C)=S", "trimethylarsane sulfide"),
    ("C[S](C)=O",        "dimethyl sulfoxide"),
    ("C[S](C)(=O)=O",    "dimethyl sulfone"),
    ("CS(=O)CC",          "ethyl methyl sulfoxide"),
])
def test_phase884_imine_family(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase884_sulfoximine_not_confused_with_sulfoxide():
    assert smiles_to_iupac("C[S](C)(=O)=N") != smiles_to_iupac("C[S](C)=O")


def test_phase884_phosphanimine_not_confused_with_phosphane():
    assert smiles_to_iupac("C[P](C)(C)=N") != smiles_to_iupac("C[P](C)(C)")
