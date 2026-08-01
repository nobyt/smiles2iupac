"""Phase 882: thiocarbonic acid family convention cleanup (IUPAC 2013 P-65.2.3.2).

Carbonic acid HO-C(=O)-OH has 3 replaceable chalcogen positions (1 double
bond + 2 single, acidic, O/S-H). Substituting O->S at each position gives 6
distinct structures (by total S count 0/1/1/2/2/3, with 2 structural isomers
each at the 1 and 2 levels). Before this phase, only 3 of the 6 were
handled, and one of THOSE was named wrong:

  - OC(=O)O  (0 S)                    -> "carbonic acid"            (OK)
  - OC(=S)S  (2 S: double + 1 single) -> "carbonodithioic O-acid"   (WRONG:
    real name is "carbonodithioic O,S-acid" -- confirmed against ChEBI:36958
    -style / PubChem / ChemSpider entries. The single-letter O-acid/S-acid
    tag is Phase 149's convention for CARBOXYLIC thioic acids (R-C(=O)SH vs
    R-C(=S)OH, where R disambiguates which single position needs tagging).
    Carbonic acid's family has NO carbon substituent -- both non-double
    positions are chalcogen-bearing, so IUPAC tags BOTH of them explicitly:
    O,O- / O,S- / S,S-.)
  - SC(=S)S  (3 S, all positions)     -> "trithiocarbonic acid"     (OK)

The other 3 structures fell through to the generic carboxylic-thioic-acid
namer, which assumes a carbon substituent (R-C(=X)YH) and produced
chemically wrong names implying a real C-H that doesn't exist:

  - OC(=S)O (1 S, double bond)   was "methanethioic O-acid" (implies H-C
    bonded, i.e. thioformic acid's O-tautomer -- WRONG, this carbon has 3
    heavy neighbors (O, O, =S), no room for an H)
  - OC(=O)S (1 S, one single)    was "methanethioic S-acid" (same problem)
  - SC(=O)S (2 S, both singles)  was "methanethioic S-acid" (same problem,
    additionally collided with the OC(=O)S case above -- two different
    compounds, one wrong name)

Fixed by adding all 4 missing/wrong structures to the top-level retained-
name table (_RETAINED_NAMES in __init__.py, canonical-SMILES-keyed, already
the mechanism used for the 3 pre-existing entries -- this is a closed set of
exactly 6 possible structures, not a general detector, matching how the
pre-existing carbonic-acid-family retained names were already implemented).
"""

import pytest

from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # the full 0/1/1/2/2/3-thio family, all now internally consistent
    ("OC(=O)O", "carbonic acid"),
    ("OC(=S)O", "carbonothioic O,O-acid"),    # double bond is S
    ("OC(=O)S", "carbonothioic O,S-acid"),    # one single is S
    ("SC(=O)S", "carbonodithioic S,S-acid"),  # both singles are S
    ("OC(=S)S", "carbonodithioic O,S-acid"),  # double + one single are S (was mistagged "O-acid")
    ("SC(=S)S", "trithiocarbonic acid"),      # all 3 positions S
])
def test_phase882_thiocarbonic_family(smiles, expected):
    assert smiles_to_iupac(smiles) == expected


def test_phase882_not_confused_with_thioformic_acid():
    # real thioformic acid (H-C(=X)YH, an actual C-H) must be unaffected --
    # different compound, different name, same generic thioic-acid namer.
    assert smiles_to_iupac("C(=O)S") == "methanethioic S-acid"
    assert smiles_to_iupac("C(=S)O") == "methanethioic O-acid"
    # and must differ from the carbonic-family compounds of the same total-S count
    assert smiles_to_iupac("OC(=O)S") != smiles_to_iupac("C(=O)S")
    assert smiles_to_iupac("OC(=S)O") != smiles_to_iupac("C(=S)O")


def test_phase882_esters_unaffected():
    # RO-C(=X)-OR'/SR' esters (Phase 348, a different group_type) are
    # untouched -- this phase only concerns the free (fully protonated) acid.
    assert smiles_to_iupac("COC(=S)OC") == "dimethyl carbonothioate"
    assert smiles_to_iupac("CCOC(=S)OCC") == "diethyl carbonothioate"
