"""Phase 891: cationic aromatic N-heterocycles (imidazolium, pyridinium)
-- the "-ium" suffix for a ring nitrogen bearing a formal +1 charge.

Found via the same fresh probe sweep as Phase 889/890. Cn1cc[n+](C)c1
(1,3-dimethylimidazolium, the core of the extremely common ionic-liquid /
NHC-precursor cation family) was named "1,3-dimethylimidazole" -- the
formal +1 charge on the ring nitrogen was silently dropped, so the cation
was reported identical to the neutral parent heterocycle. Same bug applies
to N-alkylpyridinium salts (C[n+]1ccccc1 -> wrongly "1-methylpyridine").

Root cause: the ring-composition signature function (_atom_sig) returns
plain "N" regardless of formal charge, which is actually fine for finding
the correct retained-name PARENT (imidazole/pyridine) and for numbering --
but nothing downstream ever checked for the charge afterward, so it was
lost between detection and the final name string.

Fixed with a small post-processing step, right after the base ring name
and its locant map are established: if any ring atom carries formal charge
+1, append "-{loc}-ium" (eliding the parent name's trailing "e" first,
matching the vowel-suffix elision used throughout this codebase), using
the SAME locant_map already used for substituent numbering -- so this
composes for free with existing N-alkyl substituent prefixes. Verified via
OPSIN parse-back, including the well-known BMIM+ ionic liquid cation
(1-butyl-3-methylimidazolium).
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    ("Cn1cc[n+](C)c1",       "1,3-dimethylimidazol-3-ium"),
    ("C[n+]1ccccc1",         "1-methylpyridin-1-ium"),
    ("CCCC[n+]1ccn(C)c1",    "1-butyl-3-methylimidazol-1-ium"),
    # regression: plain neutral heterocycles unchanged
    ("c1ccncc1",   "pyridine"),
    ("Cn1ccnc1",   "1-methylimidazole"),
])
def test_phase891_cationic_heterocycle_ium(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase891_imidazolium_not_confused_with_neutral_imidazole():
    assert smiles_to_iupac("Cn1cc[n+](C)c1") != smiles_to_iupac("Cn1ccnc1")


def test_phase891_pyridinium_not_confused_with_neutral_pyridine():
    result = smiles_to_iupac("C[n+]1ccccc1")
    assert "ium" in result
