"""Phase 935: hydrate/multi-water multi-component naming produced non-standard
or outright invalid names.

Track B (docs/nomenclature_gap_audit.md) found this via targeted probing of
multi-component (dot-notation) SMILES beyond the original Phase 147 scope
(which only tested charged salt components, e.g. "calcium diacetate").
Adding water as a neutral component alongside a salt exposed that
`_name_multicomponent()`'s generic `_compact()` multiplier logic (designed
for cation/anion names) also ran on neutral fragments, producing
"calcium diacetate diwater" -- non-standard (real hydrate nomenclature uses
"dihydrate", not "diwater"). Testing the same generic logic on two
disconnected neutral acid molecules (no water, no salt) produced
"diacetic acid", which is not even OPSIN-parseable ("diacetic acid" is not
a valid IUPAC name -- OPSIN 2.9.0 returns no result for it), a genuinely
invalid name, not just non-standard.

Fixed by special-casing water (IUPAC P-16.3.3 addition nomenclature):
water fragments are pulled out of the generic neutral-component multiplier
logic and appended as a single "{multiplier}hydrate" suffix at the end of
the assembled name, using "mono" for count=1 (not the bare/empty multiplier
used elsewhere) since "monohydrate" (not "hydrate") is the IUPAC/CAS-
standard term. Verified via OPSIN round-trip for 1/2/3-water cases.

The non-water case ("diacetic acid" for two disconnected identical neutral
organic molecules with no salt/hydrate relationship) is NOT fixed here --
deliberately out of scope: it's a synthetic edge case (zero real-world hits
in tests/pubchem_cache.json's 6,978 entries) and the correct IUPAC
multi-component convention for that scenario isn't as clearly established
as the water case, unlike "dihydrate" which is unambiguous, extremely
common, real-world terminology.
"""
from smiles2iupac import smiles_to_iupac


class TestPhase935HydrateNaming:
    def test_calcium_acetate_dihydrate(self):
        assert (
            smiles_to_iupac("O.O.CC(=O)[O-].[Ca+2].[O-]C(C)=O")
            == "calcium diacetate dihydrate"
        )

    def test_sodium_acetate_monohydrate(self):
        assert smiles_to_iupac("O.CC(=O)[O-].[Na+]") == "sodium acetate monohydrate"

    def test_sodium_acetate_trihydrate(self):
        assert (
            smiles_to_iupac("O.O.O.[Na+].[O-]C(=O)C")
            == "sodium acetate trihydrate"
        )

    def test_sodium_chloride_unaffected(self):
        assert smiles_to_iupac("[Na+].[Cl-]") == "sodium chloride"

    def test_calcium_diacetate_no_water_unaffected(self):
        assert (
            smiles_to_iupac("[Ca+2].[O-]C(=O)C.[O-]C(=O)C")
            == "calcium diacetate"
        )

    def test_plain_acetic_acid_unaffected(self):
        assert smiles_to_iupac("CC(=O)O") == "acetic acid"
