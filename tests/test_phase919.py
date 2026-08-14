"""Phase 919: continuing the systematic sweep for "siloed heteroaromatic
ring namer" bugs (same shape as Phase 917/918) found two more instances.

1. _name_dioic_acid's heteroaromatic branch (for rings bearing TWO
   carboxylic acid groups, e.g. pyridinedicarboxylic acids) built a bare
   "{ring}-{loc1},{loc2}-dicarboxylic acid" string and never collected any
   other ring substituent:

       OC(=O)c1cc(ncc1C(=O)O)Cl -> wrongly "pyridine-3,4-dicarboxylic acid"
           (should be "6-chloropyridine-3,4-dicarboxylic acid")

2. _name_secondary_tertiary_amine's heteroaromatic branch (used whenever
   the amine nitrogen has an N-substituent, e.g. N-methyl -- a PRIMARY
   amine on the same ring position works correctly via a different code
   path) had BOTH bugs: no locant-direction flip-check (always used
   whichever ring-traversal direction _match_retained happened to return,
   never trying the reverse to see if it gives the amine itself a lower
   locant) AND no ring-substituent collection at all:

       CNc1ccncc1Cl -> wrongly "N-methylpyridin-4-amine"
           (chlorine dropped; should be
           "3-chloro-N-methylpyridin-4-amine")

Both fixed by collecting substituents via _collect_hetero_substituents
(mirroring the Phase 917/918 fix pattern) and, for the amine case,
picking whichever ring rotation gives the LOWEST locant to the amine's
own ring-attachment carbon (the amine is the suffix/principal group, so
per IUPAC numbering priority it must get the lower locant, with
substituent locants only tie-breaking after that -- confirmed the fix
correctly implements this priority, unlike the untouched PRIMARY-amine
pathway, which was separately found during this investigation to still
get this priority ordering wrong: see project memory for that larger,
deliberately-unfixed finding in the general ring-suffix numbering
pipeline).

All fixed names round-trip verified via OPSIN + RDKit canonical-SMILES
matching.
"""

from smiles2iupac import smiles_to_iupac


class TestDioicAcidHeteroaromaticSubstituents:
    def test_chloropyridine_dicarboxylic_acid(self):
        assert (
            smiles_to_iupac("OC(=O)c1cc(ncc1C(=O)O)Cl")
            == "6-chloropyridine-3,4-dicarboxylic acid"
        )

    def test_plain_pyridine_dicarboxylic_acid_regression(self):
        assert (
            smiles_to_iupac("OC(=O)c1ccncc1C(=O)O")
            == "pyridine-3,4-dicarboxylic acid"
        )

    def test_benzene_dicarboxylic_acid_regression(self):
        assert smiles_to_iupac("OC(=O)c1ccccc1C(=O)O") == "phthalic acid"


class TestSecondaryTertiaryAmineHeteroaromaticSubstituents:
    def test_n_methyl_chloropyridine_amine(self):
        assert (
            smiles_to_iupac("CNc1ccncc1Cl") == "3-chloro-N-methylpyridin-4-amine"
        )

    def test_n_methyl_chloropyridine_amine_other_position(self):
        assert (
            smiles_to_iupac("CNc1ccc(Cl)cn1") == "5-chloro-N-methylpyridin-2-amine"
        )

    def test_plain_n_methylpyridinamine_regression(self):
        assert smiles_to_iupac("CNc1ccccn1") == "N-methylpyridin-2-amine"

    def test_n_methylaniline_regression_non_hetero(self):
        assert smiles_to_iupac("CNc1ccccc1") == "N-methylaniline"

    def test_n_n_dimethylaniline_regression(self):
        assert smiles_to_iupac("CN(C)c1ccccc1") == "N,N-dimethylaniline"
