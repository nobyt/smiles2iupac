"""Phase 921: branching polyols (trimethylolpropane, pentaerythritol type)
got confidently WRONG or outright INVALID names.

    CCC(CO)(CO)CO (trimethylolpropane) -> wrongly "butane-1,2,3-triol"
        (butane-1,2,3-triol is a real but DIFFERENT, unrelated compound --
        confirmed structurally distinct via RDKit canonical SMILES; the
        correct name is "2-ethyl-2-(hydroxymethyl)propane-1,3-diol")
    OCC(CO)(CO)CO (pentaerythritol) -> wrongly "propane-1,3-tetraol", an
        outright INVALID name (OPSIN rejects it; a propane chain cannot
        carry 4 -ol suffix locants)

Root cause (three-part chain, all in the acyclic naming pipeline):

1. functional_group.py's aggregate_groups() merges N identical alcohol
   instances into a single diol/triol/tetraol/pentaol/hexaol
   FunctionalGroup purely by COUNT, with no check that all N instances
   can ever share one simple carbon chain. On a branching (neopentyl-type)
   skeleton, a quaternary carbon with 3+ CH2OH arms can put at most 2 of
   those arms on any single simple path (a path can only use 2 of a
   carbon's bonds), so "all N" is structurally impossible whenever N >= 3
   on such a center.

2. chain_finder.py's find_principal_chain() required ALL of the merged
   group's carbons to be a SUBSET of the chosen path
   (required_carbons.issubset(set(p))). When no path could satisfy that
   (the branching case above), it silently discarded the requirement
   entirely and fell back to considering every simple path with no
   preference for how many of the required carbons it contained --
   violating IUPAC P-44.1.1 rule (a) (maximize principal-characteristic-
   group atoms on the chain) and picking an arbitrary chain that could
   contain as few as ZERO of the alcohol carbons.

3. __init__.py's _name_acyclic() then tried to compute suffix locants by
   walking pgrp.atom_indices (all N O/C pairs) against the chosen chain's
   locant_map. Any O/C pair that landed off-chain contributed nothing, so
   for trimethylolpropane only 1 of the 3 alcohols resolved to a locant,
   giving suffix_locants_list of length 1 -- which then hit
   `suffix_locants = ... if len(...) > 1 else None`, collapsing to None.
   name_assembler.py's triol/tetraol/... formatters fall back to
   HARDCODED locants ([1,2,3] for triol, etc.) whenever suffix_locants is
   None, fabricating a locant set with no relation to the real structure.
   Separately, pgrp.atom_indices (all N instances) was passed whole as
   the substituent-exclusion set to collect_substituents, so the 2 (or 3)
   off-chain CH2OH branches were excluded from ordinary substituent
   scanning too -- silently vanishing from the name entirely rather than
   appearing as "hydroxymethyl" substituents.

Fixed in two places:

- chain_finder.py: replaced the issubset-or-give-up-entirely filter with
  an overlap-MAXIMIZING filter (`len(required_carbons & set(p))`, keep
  only paths achieving the max). This is a strict generalization: when a
  full-subset path exists, only those paths achieve the maximum, so
  behavior for every previously-passing molecule is unchanged. When no
  full-subset path exists, the chain that carries the MOST alcohol
  carbons is chosen instead of an arbitrary chain that might carry none
  -- correctly implementing IUPAC P-44.1.1 rule (a) for the
  partial-fit case.

- __init__.py's _name_acyclic(): after the chain is chosen, re-counts how
  many of the merged group's (C, O) instance pairs actually landed on the
  chain. If fewer than the merged count, the FunctionalGroup is rebuilt
  with the correct smaller merge type (triol -> diol, tetraol -> diol,
  etc., down to a bare "alcohol" if only 1 instance fits) and
  atom_indices trimmed to only the on-chain pairs. This both fixes the
  suffix locants (now computed from a truthful, fully-on-chain group) and
  releases the off-chain CH2OH branches from the substituent-exclusion
  set, so they get named through the existing, already-verified generic
  substituent path (substituent.py's hydroxymethyl/2-hydroxyethyl logic,
  previously only exercised for non-principal-group alcohols e.g. on
  cycloalkanes -- Phase 274) as ordinary "(hydroxymethyl)" /
  "bis(hydroxymethyl)" / "tris(hydroxymethyl)" prefixes.

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching. Full pre-existing suite re-run in isolation before this test
file existed, to confirm zero regressions from the two source changes
alone (both are structural generalizations of existing, narrower logic,
not new special-casing, so every already-passing molecule -- including
plain glycerol, ethylene glycol, neopentyl glycol -- keeps its prior,
already-correct chain choice and name unchanged).
"""

from smiles2iupac import smiles_to_iupac


class TestBranchingPolyolDowngrade:
    def test_trimethylolpropane(self):
        assert (
            smiles_to_iupac("CCC(CO)(CO)CO")
            == "2-ethyl-2-(hydroxymethyl)propane-1,3-diol"
        )

    def test_pentaerythritol(self):
        assert (
            smiles_to_iupac("OCC(CO)(CO)CO")
            == "2,2-bis(hydroxymethyl)propane-1,3-diol"
        )

    def test_trimethylolmethane_with_extra_methyl(self):
        assert (
            smiles_to_iupac("OCC(C)(CO)CO")
            == "2-(hydroxymethyl)-2-methylpropane-1,3-diol"
        )

    def test_single_branch_triol_downgrades_to_diol(self):
        assert smiles_to_iupac("OCC(CO)CO") == "2-(hydroxymethyl)propane-1,3-diol"

    def test_two_quaternary_centers_hexaol_to_diol(self):
        assert (
            smiles_to_iupac("OCC(CO)(CO)C(CO)(CO)CO")
            == "2,2,3,3-tetrakis(hydroxymethyl)butane-1,4-diol"
        )


class TestBranchingPolyolRegressions:
    def test_ethylene_glycol(self):
        assert smiles_to_iupac("OCCO") == "ethane-1,2-diol"

    def test_glycerol_unbranched_triol_unaffected(self):
        assert smiles_to_iupac("OCC(O)CO") == "propane-1,2,3-triol"

    def test_neopentyl_glycol_plain_diol_unaffected(self):
        assert smiles_to_iupac("CC(C)(CO)CO") == "2,2-dimethylpropane-1,3-diol"
