"""Phase 889: two bugs found in the same fresh probe sweep as 887/888.

1. Weinreb amides (R-C(=O)-N(R')-O-R'') -- CC(=O)N(C)OC (N-methoxy-N-
   methylacetamide, an extremely common synthesis reagent for controlled
   ketone synthesis from amides) was named "N-methylacetamide", silently
   dropping the whole N-methoxy group and colliding with the real
   N-methylacetamide (CC(=O)NC). Root cause: the amide namer's N-hydroxy
   detection (for hydroxamic acids, Phase 61) only handled a bare -OH on
   the amide nitrogen; when the O instead carried a carbon substituent
   (N-O-R, not N-O-H), it fell through every check and was silently
   dropped. Fixed by also recognizing N-O-R and naming it "{alkyl}oxy"
   (reusing the existing _make_oxy_name helper), landing it in the same
   N-/N,N- prefix list as ordinary N-alkyl substituents.

2. Isourea/isothiourea (the carbamimidic acid family, O/S esters and free
   acids of H2N-C(=NH)-OH) -- CN=C(N)OC and CN=C(N)SC (two genuinely
   different compounds) were BOTH named "N'-methylmethanimidamide",
   colliding with each other and dropping the O-methyl/S-methyl group
   entirely. Root cause: _is_amidine only checked for an imine-N and an
   amine-N on the central carbon, never checking for an extra O/S
   substituent that would make it isourea/isothiourea instead of a plain
   amidine. Fixed with a guard, then implemented the correct carbamimidic-
   acid-family group properly (new group types carbamimidic_acid/
   carbamimidate_ester/carbamimidothioic_acid/carbamimidothioate_ester)
   since disabling the wrong amidine match alone still left them falling
   into a different wrong fallback ("1-aminoformamide"-style dropped
   info). Amine-N substituents get "N-", imine-N substituents get "N'-"
   (OPSIN-verified: "methyl N'-methylcarbamimidate" -> CN=C(N)OC).

   While implementing #2, caught and fixed two self-introduced bugs before
   they ever landed: the N-substituent scan and the O/S-alkyl scan both
   initially forgot to exclude the central carbamimidic carbon itself from
   "is this a real substituent" (since that carbon IS a same-symbol
   carbon neighbor of both the nitrogens AND the ester oxygen), causing
   the namer to try to name the rest of the whole molecule as a
   substituent of itself. A reminder that any O/N-neighbor scan on a
   pivot atom that then loops back to the SAME functional group's own
   central atom needs that atom excluded explicitly, not just "is it a
   plausible substituent type."
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Weinreb amide and regressions
    ("CC(=O)N(C)OC", "N-methoxy-N-methylacetamide"),
    ("CC(=O)NOC",    "N-methoxyacetamide"),
    ("CC(=O)NO",     "N-hydroxyacetamide"),
    ("CC(=O)NC",     "N-methylacetamide"),
    ("CC(=O)N(C)C",  "N,N-dimethylacetamide"),
    # carbamimidic acid family
    ("NC(=N)O",     "carbamimidic acid"),
    ("NC(=N)S",     "carbamimidothioic acid"),
    ("CN=C(N)OC",   "methyl N'-methylcarbamimidate"),
    ("COC(=N)NC",   "methyl N-methylcarbamimidate"),
    ("CN=C(N)SC",   "methyl N'-methylcarbamimidothioate"),
    ("COC(N)=N",    "methyl carbamimidate"),
    ("CSC(N)=N",    "methyl carbamimidothioate"),
    # regression: plain amidine/guanidine/imidic acid/imidate ester unchanged
    ("CC(=N)N",    "ethanimidamide"),
    ("NC(=N)N",    "guanidine"),
    ("CC(=N)O",    "ethanimidic acid"),
    ("CC(=N)OCC",  "ethyl ethanimidate"),
])
def test_phase889_weinreb_and_carbamimidic(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase889_weinreb_not_confused_with_n_methylacetamide():
    assert smiles_to_iupac("CC(=O)N(C)OC") != smiles_to_iupac("CC(=O)NC")


def test_phase889_isourea_not_confused_with_isothiourea():
    assert smiles_to_iupac("CN=C(N)OC") != smiles_to_iupac("CN=C(N)SC")
