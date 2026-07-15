"""Phase 868: phosphoramidic / phosphonamidic / phosphorodiamidic acids.

Amide (P-N) derivatives of phosphoric/phosphonic acid dropped the entire
phosphorus group: CNP(=O)(O)O was named "aminomethane", CP(=O)(N)O was
mis-named "methylphosphinic acid".

Naming (IUPAC 2013 P-67.1.4): amide-N substituents take N-/N'- locants, and
the C substituent on P takes a P- locant (phosphonamidic acid has both P and
N substitutable positions, so locants are required):

  H2N-P(=O)(OH)2 substituted on N -> {N-...}phosphoramidic acid
  (H2N)2-P(=O)(OH)                -> {N-.../N'-...}phosphorodiamidic acid
  R-P(=O)(NH2)(OH)                -> {P-R}phosphonamidic acid

The carbon-free parent H2N-P(=O)(OH)2 is out of scope (the namer requires at
least one carbon), so only N-substituted / P-substituted members are covered.
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phosphoramidic (0 C on P, 1 amide N, 2 OH)
    ("CNP(=O)(O)O",     "N-methylphosphoramidic acid"),
    ("CN(C)P(=O)(O)O",  "N,N-dimethylphosphoramidic acid"),
    ("CCNP(=O)(O)O",    "N-ethylphosphoramidic acid"),
    # phosphonamidic (1 C on P, 1 amide N, 1 OH)
    ("CP(=O)(N)O",      "P-methylphosphonamidic acid"),
    ("CCP(=O)(N)O",     "P-ethylphosphonamidic acid"),
    ("CP(=O)(NC)O",     "N,P-dimethylphosphonamidic acid"),
    # phosphorodiamidic (0 C on P, 2 amide N, 1 OH)
    ("CNP(=O)(N)O",     "N-methylphosphorodiamidic acid"),
    ("CNP(=O)(NC)O",    "N,N'-dimethylphosphorodiamidic acid"),
    # regression: oxy / thio phosphorus acids and esters unchanged
    ("CP(=O)(O)O",      "methylphosphonic acid"),
    ("CP(C)(=O)O",      "dimethylphosphinic acid"),
    ("CP(O)O",          "methylphosphonous acid"),
    ("CP(=S)(O)O",      "methylphosphonothioic acid"),
    ("CP(C)(C)=O",      "trimethylphosphane oxide"),
    ("CP(=O)(OC)OC",    "dimethyl methylphosphonate"),
])
def test_phase868_phosphoramidic_family(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
