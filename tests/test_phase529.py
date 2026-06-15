"""Phase 529: Heteroaromatic azide naming (IUPAC 2013 substitutive nomenclature).

Rings bearing -N3 were routed through find_principal_ring/assemble_ring_name
which treated pyridine as benzene ("azidobenzene") and thiophene as
cyclopentane ("azidocyclopentane"). Now correctly named as
"{loc}-azido{hetname}" per substitutive convention.
"""
import pytest
from src.smiles2iupac import smiles_to_iupac


@pytest.mark.parametrize("smiles,expected", [
    # Heteroaromatic azides
    ("[N-]=[N+]=Nc1ccccn1",   "2-azidopyridine"),
    ("[N-]=[N+]=Nc1cccnc1",   "3-azidopyridine"),
    ("[N-]=[N+]=Nc1ccncc1",   "4-azidopyridine"),
    ("[N-]=[N+]=Nc1cccs1",    "2-azidothiophene"),
    ("[N-]=[N+]=Nc1ccco1",    "2-azidofuran"),
    ("[N-]=[N+]=Nc1ccc[nH]1", "2-azido-1H-pyrrole"),
    # Regression: benzene unchanged
    ("[N-]=[N+]=Nc1ccccc1",   "azidobenzene"),
    # Regression: aliphatic azides unchanged
    ("[N-]=[N+]=NC",          "azidomethane"),
    ("[N-]=[N+]=NCC",         "azidoethane"),
])
def test_phase529_heteroaromatic_azide(smiles, expected):
    assert smiles_to_iupac(smiles) == expected
