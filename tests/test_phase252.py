"""Phase 252: inorganic phosphorus/sulfur halides and complete silicon halides.

IUPAC 2013 retained/accepted names for common inorganic halides:
  PCl3 → phosphorus trichloride; POCl3 → phosphoryl trichloride
  SOCl2 → thionyl chloride; SO2Cl2 → sulfuryl chloride
  SiHCl3 → trichlorosilane; SiCl4 → tetrachlorosilane
"""

import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # phosphorus trihalides
    ("ClP(Cl)Cl",       "phosphorus trichloride"),
    ("FP(F)F",          "phosphorus trifluoride"),
    ("BrP(Br)Br",       "phosphorus tribromide"),
    ("IP(I)I",          "phosphorus triiodide"),
    # phosphoryl / phosphorothioyl
    ("ClP(=O)(Cl)Cl",   "phosphoryl trichloride"),
    ("FP(=O)(F)F",      "phosphoryl trifluoride"),
    ("ClP(=S)(Cl)Cl",   "phosphorothioyl trichloride"),
    # thionyl / sulfuryl
    ("ClS(Cl)=O",       "thionyl chloride"),
    ("FS(F)=O",         "thionyl fluoride"),
    ("ClS(Cl)(=O)=O",   "sulfuryl chloride"),
    ("FS(F)(=O)=O",     "sulfuryl fluoride"),
    # silicon halides
    ("[SiH](Cl)(Cl)Cl", "trichlorosilane"),
    ("[Si](Cl)(Cl)(Cl)Cl", "tetrachlorosilane"),
    ("[SiH2](F)(F)",    "difluorosilane"),
    # regression: previously existing entries unchanged
    ("[SiH3]Cl",        "chlorosilane"),
    ("Cl[SiH2]Cl",      "dichlorosilane"),
    ("[PH3]",           "phosphane"),
])
def test_phase252_inorganic_halides(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
