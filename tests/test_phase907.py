"""Phase 907: acylhydrazones (R-C(=O)-NH-N=CR'R'', the hydrazide+ketone/
aldehyde condensation product -- a very common motif in medicinal chemistry,
e.g. isoniazid-derived hydrazones) were completely unhandled, plus three
related chain-substituent-dropping bugs found while investigating the
hydrazone/hydrazide family.

Root cause of the main bug: `_is_hydrazide` already excludes the case where
the far nitrogen (N2) carries a C=N double bond, deferring to the
semicarbazone/hydrazone detectors on the *other* carbon (the imine carbon).
But semicarbazone detection additionally requires a urea-style NH2 on the
acyl carbon, so a *plain* acylhydrazone (no urea NH2) matched neither
semicarbazone nor hydrazide. It fell through to plain `_is_amide` at the
acyl carbon (priority=95) while the imine carbon separately registered a
`kethydrazone`/`aldhydrazone` (priority=53) -- amide won the priority race
and the entire hydrazone side silently vanished, e.g.
`CC(=O)N/N=C/c1ccccc1` (acetohydrazide benzaldehyde hydrazone) was named
bare "acetamide".

Fixed by adding a new `_is_acylhydrazone_or_thio` detector (mirrors
`_is_semicarbazone_or_thio` but without requiring the urea NH2), new group
types acylhydrazone/aldacylhydrazone/thioacylhydrazone/aldthioacylhydrazone
(priority=96, beats amide), and a dedicated `_name_acylhydrazone` namer that
builds the OPSIN-verified "N'-{alkylidene}...(thio)hydrazide" form (parent =
hydrazide side, ketone/aldehyde side becomes an N'-alkylidene substituent).

Three more bugs found and fixed along the way (same chain-substituent-
dropping shape as the whole 902-906 arc, just in previously-unswept
namers):

1. `_name_semicarbazone` (semicarbazone/aldsemicarbazone/thiosemicarbazone/
   aldthiosemicarbazone) and `_name_substituted_hydrazone` (N-substituted
   kethydrazone/aldhydrazone) never called collect_substituents on the
   ketone/aldehyde chain -- `FC(F)(F)CC(C)=NNC(N)=O` (a CF3-substituted
   semicarbazone) collapsed to "butan-2-one semicarbazone", dropping the
   CF3 entirely.

2. `_name_thiohydrazide`/`_name_selenohydrazide` had the identical gap on
   the acid chain -- `FC(F)(F)C(=S)NN` collapsed to "ethanethiohydrazide".

3. `_alkylidene_name` (used for N-aryl Schiff bases, e.g. anilines) only
   ever counted chain LENGTH via a naive DFS, never substituents --
   `FC(F)(F)CC(C)=Nc1ccccc1` collapsed to "N-butylideneaniline", dropping
   the CF3. Rewritten to use find_principal_chain + collect_substituents
   like every other chain-based namer, which also gave `_name_acylhydrazone`
   a correct alkylidene-name builder to reuse. Needed a follow-up fix in
   `_name_n_substituted_imine`'s aniline branch to parenthesize the
   alkylidene name when it now contains locants (`_needs_bis_tris`), since
   a bare "N-4,4,4-trifluorobutan-2-ylideneaniline" would be ambiguous.

All OPSIN-verified.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # --- acylhydrazone (new detection + namer) ---
    ("CC(=O)N/N=C/c1ccccc1", "(E)-N'-benzylideneacetohydrazide"),
    ("CC(=O)NN=C", "N'-methylideneacetohydrazide"),
    ("CC(=O)N/N=C(C)/C", "N'-(propan-2-ylidene)acetohydrazide"),
    ("CC(=O)N/N=C(C)/CC(F)(F)F",
     "(E)-N'-(4,4,4-trifluorobutan-2-ylidene)acetohydrazide"),
    ("O=C(NN=Cc1ccccc1)c1ccccc1", "N'-benzylidenebenzohydrazide"),
    ("FC(F)(F)C(=O)N/N=C/c1ccccc1",
     "(E)-N'-benzylidenetrifluoroacetohydrazide"),
    # thio variant
    ("CC(=S)N/N=C/c1ccccc1", "(E)-N'-benzylideneethanethiohydrazide"),
    ("CC(=S)N/N=C(C)/CC(F)(F)F",
     "(E)-N'-(4,4,4-trifluorobutan-2-ylidene)ethanethiohydrazide"),
    # regression: plain hydrazide / plain hydrazone / semicarbazone untouched
    ("CC(=O)NN", "acetohydrazide"),
    ("CC=NN", "ethanal hydrazone"),
    ("CC(=NNC(N)=O)C", "propan-2-one semicarbazone"),
])
def test_phase907_acylhydrazone(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase907_acylhydrazone_not_confused_with_plain_amide():
    assert smiles_to_iupac("CC(=O)N/N=C/c1ccccc1") != "acetamide"


@pytest.mark.parametrize("smiles,expected", [
    ("FC(F)(F)CC(C)=NNC(N)=O", "4,4,4-trifluorobutan-2-one semicarbazone"),
    ("FC(F)(F)CC(C)=NNC(N)=S", "4,4,4-trifluorobutan-2-one thiosemicarbazone"),
    ("FC(F)(F)CC(C)=NNC", "4,4,4-trifluorobutan-2-one N-methylhydrazone"),
    # regressions
    ("CC(=NNC(N)=O)C", "propan-2-one semicarbazone"),
    ("CC=NNC(N)=O", "ethanal semicarbazone"),
    ("CC(=NNC(N)=S)C", "propan-2-one thiosemicarbazone"),
    ("CC=NNC", "ethanal N-methylhydrazone"),
    ("CC(=NNC)C", "propan-2-one N-methylhydrazone"),
])
def test_phase907_semicarbazone_and_n_hydrazone_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("FC(F)(F)C(=S)NN", "trifluoroethanethiohydrazide"),
    ("FC(F)(F)C(=[Se])NN", "trifluoroethaneselenohydrazide"),
    # regressions
    ("CC(=S)NN", "ethanethiohydrazide"),
    ("CC(=[Se])NN", "ethaneselenohydrazide"),
    ("O=C(NN)c1ccccc1", "benzohydrazide"),
])
def test_phase907_thio_seleno_hydrazide_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


@pytest.mark.parametrize("smiles,expected", [
    ("FC(F)(F)CC(C)=Nc1ccccc1", "N-(4,4,4-trifluorobutan-2-ylidene)aniline"),
    # regressions
    ("C=Nc1ccccc1", "N-methylideneaniline"),
    ("CC=Nc1ccccc1", "N-ethylideneaniline"),
    ("CN=Cc1ccccc1", "N-methylphenylmethanimine"),
    ("CN=C(C)C", "N-methylpropan-2-imine"),
])
def test_phase907_alkylidene_name_chain_substituent(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
