"""Phase 901: selenocyanate/tellurocyanate family (Se-C#N / Te-C#N) --
the Se/Te mirror of thiocyanate/cyanate, closing out the gap flagged (but
deferred pending investigation) at the end of the Phase 900 sweep.

Se-C#N / Te-C#N had NO functional-group detection at all -- every form
(anion, neutral acid, ester) fell through to a generic "selanyl/tellanyl"
substituted-nitrile name, which additionally collided the anion with the
neutral acid exactly like the pre-Phase-900 cyanate/thiocyanate bug:
  [Se-]C#N (selenocyanate anion)  -> "selanylmethanenitrile"
  [SeH]C#N (selenocyanic acid)    -> "selanylmethanenitrile"  (same name!)
  C[Se]C#N (methyl ester)         -> "methylselanylmethanenitrile"
                                      (inconsistent style vs thiocyanate's
                                      "thiocyanatomethane")

Fixed by generalizing the existing thiocyanate detection helpers
(_is_chalcocyanate/_get_chalcocyanate_chalcogen, parameterized on the
chalcogen symbol) and adding selenocyanate/tellurocyanate group types with
namers mirroring _name_thiocyanate exactly (ester -> "{chalco}cyanato-
{alkane}", free acid vs anion distinguished by formal_charge, matching the
Phase 900 fix). Also added both new group types to chain_finder.py's
nitrile_types set so the ester form's chain-building excludes the cyano
carbon like the existing cyanate/thiocyanate forms already do. Verified
via OPSIN.
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("[Se-]C#N", "selenocyanate"),
    ("[SeH]C#N", "selenocyanic acid"),
    ("C[Se]C#N", "selenocyanatomethane"),
    ("[Te-]C#N", "tellurocyanate"),
    ("[TeH]C#N", "tellurocyanic acid"),
    ("C[Te]C#N", "tellurocyanatomethane"),
    # regression: S/O analogs (Phase 69/900) unchanged
    ("[S-]C#N", "thiocyanate"),
    ("SC#N",    "thiocyanic acid"),
    ("CSC#N",   "thiocyanatomethane"),
    ("[O-]C#N", "cyanate"),
    ("OC#N",    "cyanic acid"),
    # regression: plain nitrile unchanged
    ("CC#N", "acetonitrile"),
])
def test_phase901_selenocyanate_tellurocyanate(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase901_selenocyanate_anion_not_confused_with_acid():
    assert smiles_to_iupac("[Se-]C#N") != smiles_to_iupac("[SeH]C#N")


def test_phase901_tellurocyanate_anion_not_confused_with_acid():
    assert smiles_to_iupac("[Te-]C#N") != smiles_to_iupac("[TeH]C#N")
