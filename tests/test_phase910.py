"""Phase 910: N,O-disubstituted hydroxylamines (R-NH-O-R' / R2N-O-R') had no
detection at all — fell through both the N-substituted-only namer (requires
O to be plain -OH) and the O-substituted-only namer (requires N to have no
C neighbors) to a nonsensical generic fallback that misread the N-O bond as
a plain amine and silently dropped the O-alkyl group entirely, e.g.
CONCC(F)(F)F ("N-(2,2,2-trifluoroethyl)-O-methylhydroxylamine") -> wrongly
"2-amino-1,1,1-trifluoroethane".

Also fixes two pre-existing bugs in the sibling _name_n_substituted_hydroxylamine
found via regression testing while building this: identical-name substituents
on N joined with a comma instead of a hyphen (e.g. "N-ethyl,N-methyl..." instead
of "N-ethyl-N-methyl..."), and alphabetization sorting by the raw string
(so "2,2,2-trifluoroethyl" sorted before "methyl" since '2' < 'm' in ASCII)
instead of the substituent word itself.

All forms OPSIN-verified via round-trip SMILES canonicalization.
"""

from smiles2iupac import smiles_to_iupac


class TestNODisubstitutedHydroxylamine:
    def test_simple_n_ethyl_o_methyl(self):
        assert smiles_to_iupac("CONCC") == "N-ethyl-O-methylhydroxylamine"

    def test_symmetric_both_methyl_merges(self):
        # OPSIN-verified merged form (not "N-methyl-O-methylhydroxylamine")
        assert smiles_to_iupac("CONC") == "N,O-dimethylhydroxylamine"

    def test_cf3_on_n_side(self):
        assert (
            smiles_to_iupac("CCONCC(F)(F)F")
            == "O-ethyl-N-2,2,2-trifluoroethylhydroxylamine"
        )

    def test_n_n_dimethyl_o_methyl_all_merge(self):
        assert smiles_to_iupac("CN(C)OC") == "N,N,O-trimethylhydroxylamine"

    def test_n_n_dimethyl_o_cf3ethyl(self):
        assert (
            smiles_to_iupac("CN(C)OCC(F)(F)F")
            == "N,N-dimethyl-O-2,2,2-trifluoroethylhydroxylamine"
        )

    def test_mixed_n_substituents_with_o_methyl_tie(self):
        assert (
            smiles_to_iupac("FC(F)(F)CN(C)OC")
            == "N,O-dimethyl-N-2,2,2-trifluoroethylhydroxylamine"
        )

    def test_weinreb_amide_not_regressed(self):
        # N(OMe)(Me) on a carbonyl carbon must stay an amide, not hydroxylamine
        assert smiles_to_iupac("CC(=O)N(C)OC") == "N-methoxy-N-methylacetamide"
        assert (
            smiles_to_iupac("CC(=O)N(C)OCC(F)(F)F")
            == "N-(2,2,2-trifluoroethoxy)-N-methylacetamide"
        )


class TestHydroxylamineSiblingRegressions:
    def test_n_only_still_works(self):
        assert smiles_to_iupac("CNO") == "N-methylhydroxylamine"

    def test_o_only_still_works(self):
        assert smiles_to_iupac("NOC") == "O-methylhydroxylamine"

    def test_n_n_only_still_works(self):
        assert smiles_to_iupac("CN(O)C") == "N,N-dimethylhydroxylamine"

    def test_asymmetric_n_n_disubstituted_hyphen_not_comma(self):
        # Pre-existing bug: was "N-ethyl,N-methylhydroxylamine" (comma) instead
        # of the correct hyphen-joined, alphabetized form.
        assert smiles_to_iupac("CCN(C)O") == "N-ethyl-N-methylhydroxylamine"

    def test_n_n_disubstituted_alpha_order_ignores_embedded_locants(self):
        # Pre-existing bug: raw-string sort put "2,2,2-trifluoroethyl" before
        # "methyl" (ASCII '2' < 'm'); correct alphabetical order is by the
        # substituent word itself ("methyl" < "trifluoroethyl").
        assert (
            smiles_to_iupac("CN(CC(F)(F)F)O")
            == "N-methyl-N-2,2,2-trifluoroethylhydroxylamine"
        )
