"""Phase 619: non-aromatic double-bond compounds — 1H-indene, 2H/4H-chromene,
indan-1-one, indane-1,3-dione, isobenzofuran-1,3-dione, furan-2,5-dione, perimidine."""
import pytest
from smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # 1H-indene (positions 1–7 all distinct)
    ("C1=Cc2ccccc2C1",    "1H-indene"),
    ("CC1C=Cc2ccccc21",   "1-methyl-1H-indene"),
    ("CC1=Cc2ccccc2C1",   "2-methyl-1H-indene"),
    ("CC1=CCc2ccccc21",   "3-methyl-1H-indene"),
    ("Cc1cccc2c1C=CC2",   "4-methyl-1H-indene"),
    ("Cc1ccc2c(c1)C=CC2", "5-methyl-1H-indene"),
    ("Cc1ccc2c(c1)CC=C2", "6-methyl-1H-indene"),
    ("Cc1cccc2c1CC=C2",   "7-methyl-1H-indene"),
    # 2H-chromene (positions 2–8 all distinct)
    ("C1=Cc2ccccc2OC1",    "2H-chromene"),
    ("CC1C=Cc2ccccc2O1",   "2-methyl-2H-chromene"),
    ("CC1=Cc2ccccc2OC1",   "3-methyl-2H-chromene"),
    ("CC1=CCOc2ccccc21",   "4-methyl-2H-chromene"),
    ("Cc1cccc2c1C=CCO2",   "5-methyl-2H-chromene"),
    ("Cc1ccc2c(c1)C=CCO2", "6-methyl-2H-chromene"),
    ("Cc1ccc2c(c1)OCC=C2", "7-methyl-2H-chromene"),
    ("Cc1cccc2c1OCC=C2",   "8-methyl-2H-chromene"),
    # 4H-chromene (positions 2–8 all distinct)
    ("C1=COc2ccccc2C1",    "4H-chromene"),
    ("CC1=CCc2ccccc2O1",   "2-methyl-4H-chromene"),
    ("CC1=COc2ccccc2C1",   "3-methyl-4H-chromene"),
    ("CC1C=COc2ccccc21",   "4-methyl-4H-chromene"),
    ("Cc1cccc2c1CC=CO2",   "5-methyl-4H-chromene"),
    ("Cc1ccc2c(c1)CC=CO2", "6-methyl-4H-chromene"),
    ("Cc1ccc2c(c1)OC=CC2", "7-methyl-4H-chromene"),
    ("Cc1cccc2c1OC=CC2",   "8-methyl-4H-chromene"),
    # indan-1-one (positions 2–7 all distinct)
    ("O=C1CCc2ccccc21",    "indan-1-one"),
    ("CC1Cc2ccccc2C1=O",   "2-methylindan-1-one"),
    ("CC1CC(=O)c2ccccc21", "3-methylindan-1-one"),
    ("Cc1cccc2c1CCC2=O",   "4-methylindan-1-one"),
    ("Cc1ccc2c(c1)CCC2=O", "5-methylindan-1-one"),
    ("Cc1ccc2c(c1)C(=O)CC2","6-methylindan-1-one"),
    ("Cc1cccc2c1C(=O)CC2", "7-methylindan-1-one"),
    # indane-1,3-dione (2, 4=7, 5=6 by symmetry → lowest locant wins)
    ("O=C1CC(=O)c2ccccc21",  "indane-1,3-dione"),
    ("CC1C(=O)c2ccccc2C1=O", "2-methylindane-1,3-dione"),
    ("Cc1cccc2c1C(=O)CC2=O", "4-methylindane-1,3-dione"),
    ("Cc1ccc2c(c1)C(=O)CC2=O","5-methylindane-1,3-dione"),
    # isobenzofuran-1,3-dione (4=7, 5=6 by symmetry → lowest locant wins)
    ("O=C1OC(=O)c2ccccc21",   "isobenzofuran-1,3-dione"),
    ("Cc1cccc2c1C(=O)OC2=O",  "4-methylisobenzofuran-1,3-dione"),
    ("Cc1ccc2c(c1)C(=O)OC2=O","5-methylisobenzofuran-1,3-dione"),
    # furan-2,5-dione (3=4 by symmetry → both map to locant 3)
    ("O=C1C=CC(=O)O1",  "furan-2,5-dione"),
    ("CC1=CC(=O)OC1=O", "3-methylfuran-2,5-dione"),
    # perimidine (2 unique; 4=10, 5=9, 6=8 by symmetry → lowest wins)
    ("C1=Nc2cccc3cccc(c23)N1", "perimidine"),
    ("CC1=Nc2cccc3cccc(c23)N1","2-methylperimidine"),
    ("Cc1ccc2cccc3c2c1N=CN3",  "4-methylperimidine"),
    ("Cc1cc2c3c(cccc3c1)NC=N2","5-methylperimidine"),
    ("Cc1ccc2c3c(cccc13)NC=N2","6-methylperimidine"),
])
def test_phase619_nonaromatic_double_bond(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
