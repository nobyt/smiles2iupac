"""Phase 932: general (non-spiro/bicyclo/PAH) 3+-ring polycyclic systems
silently collapsed to a wrong, incomplete name instead of erroring.

Found while researching Track B ("what IUPAC nomenclature is genuinely
unimplemented, not just buggy") for docs/nomenclature_gap_audit.md.
`find_principal_ring`'s multi-ring fallback (ring_handler.py, "複数環: 最大環
を選択") picks whichever SSSR ring is largest and names ONLY that ring,
discarding every atom in the other ring(s) with no substituent, no error,
nothing -- `smiles_to_iupac("C1CC2CC1C1CCCCC21")` (a genuine, connected,
11-atom/3-ring tricyclic) silently returned plain "cyclohexane" (dropping
5 of 11 ring atoms), and a 15-atom tricyclic ortho-fused system silently
returned "cycloheptane".

This is not a locant-table bug like Phase 930/931 -- it's a real, larger
missing feature (general von Baeyer nomenclature for 3+-ring fused/bridged
systems beyond the specific shapes already handled: spiro, bicyclo[l.m.n],
retained cage names, naphthalene/anthracene/phenanthrene/PAH, biphenyl).
Implementing the correct general algorithm (IUPAC P-23.2.3-5: main ring/
main bridge selection, numbering direction rules for arbitrary N-ring
systems) is a substantial undertaking deliberately out of scope for this
phase. Instead, this phase converts the SILENT wrong answer into a LOUD,
honest `ValueError` -- matching this codebase's established preference
(see the diimine "butane-0,1-diimine" note in project history) for a
failure that's at least visible over one that looks like a valid answer
for the wrong molecule.
"""
import pytest
from smiles2iupac import smiles_to_iupac


class TestPhase932PolycyclicFallbackLoudFailure:
    def test_bridged_ortho_fused_tricyclic_raises(self):
        # Phase 942 implemented 3+-ring von Baeyer nomenclature (P-23.2.3 - P-23.2.5)
        assert smiles_to_iupac("C1CC2CC1C1CCCCC21") == "tricyclo[6.2.1.0²,⁷]undecane"

    def test_ortho_fused_tricyclic_raises(self):
        assert smiles_to_iupac("C1CCC2CCCC3CCCCC3C2C1") == "tricyclo[9.4.0.0²,⁷]pentadecane"

    def test_decalin_still_works(self):
        assert smiles_to_iupac("C1CCC2CCCCC2C1") == "decahydronaphthalene"

    def test_bicyclooctane_still_works(self):
        assert smiles_to_iupac("C1CC2CCC1CC2") == "bicyclo[2.2.2]octane"

    def test_plain_cyclohexane_still_works(self):
        assert smiles_to_iupac("C1CCCCC1") == "cyclohexane"

    def test_naphthalene_still_works(self):
        assert smiles_to_iupac("c1ccc2ccccc2c1") == "naphthalene"

    def test_anthracene_still_works(self):
        assert smiles_to_iupac("c1ccc2cc3ccccc3cc2c1") == "anthracene"

    def test_adamantane_still_works(self):
        assert smiles_to_iupac("C1C2CC3CC1CC(C2)C3") == "adamantane"
