"""Utilities for validating IUPAC names using the PubChem PUG REST API.

This module uses only the Python standard library so no extra dependency is needed.
"""
from __future__ import annotations

from typing import Optional
import json
import urllib.request
import urllib.parse

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def _http_get_json(url: str, timeout: float = 5.0) -> Optional[dict]:
    req = urllib.request.Request(url, headers={"User-Agent": "smiles2iupac/0.1"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = resp.read()
            return json.loads(data)
    except Exception:
        return None


def _normalize_name(s: str) -> str:
    # simple normalization for comparison: lower-case and collapse whitespace
    return " ".join(s.strip().lower().split())


def is_valid_iupac(name: str, timeout: float = 5.0) -> bool:
    """Return True if PubChem recognizes the given name and reports an IUPAC name that
    matches the provided one (best-effort).

    Strategy:
    - Lookup compound CIDs by name. If none found -> False.
    - Fetch IUPACName property for the first CID and compare (normalized). If match -> True.
    - If fetching property fails or names don't match, return False.
    """
    if not name:
        return False

    q = urllib.parse.quote(name, safe="")
    cids_url = f"{BASE}/compound/name/{q}/cids/JSON"
    j = _http_get_json(cids_url, timeout=timeout)
    if not j:
        return False
    try:
        cids = j["IdentifierList"]["CID"]
    except Exception:
        return False
    if not cids:
        return False

    cid = cids[0]
    prop_url = f"{BASE}/compound/cid/{cid}/property/IUPACName/JSON"
    pj = _http_get_json(prop_url, timeout=timeout)
    if not pj:
        return False
    try:
        prop = pj["PropertyTable"]["Properties"][0]
        pub_name = prop.get("IUPACName", "")
    except Exception:
        return False

    return _normalize_name(pub_name) == _normalize_name(name)


# rate limiting: ensure at least 200ms between requests to PubChem
_last_request_time: float | None = None
_min_interval = 0.2
import time


def _rate_limited_get_json(url: str, timeout: float = 5.0) -> Optional[dict]:
    global _last_request_time
    now = time.time()
    if _last_request_time is not None:
        elapsed = now - _last_request_time
        if elapsed < _min_interval:
            time.sleep(_min_interval - elapsed)
    res = _http_get_json(url, timeout=timeout)
    _last_request_time = time.time()
    return res


def get_inchikey_for_smiles(smiles: str, timeout: float = 5.0) -> Optional[str]:
    """Return InChIKey for the given SMILES.

    Try RDKit if available; otherwise ask PubChem for InChIKey via the SMILES endpoint.
    """
    if not smiles:
        return None
    # try RDKit first
    try:
        from rdkit import Chem  # type: ignore

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise Exception("RDKit failed to parse SMILES")
        ik = Chem.MolToInchiKey(mol)
        if ik:
            return ik
    except Exception:
        pass

    # fallback to PubChem: compound/smiles/{smiles}/property/InChIKey/JSON
    q = urllib.parse.quote(smiles, safe="")
    url = f"{BASE}/compound/smiles/{q}/property/InChIKey/JSON"
    j = _rate_limited_get_json(url, timeout=timeout)
    if not j:
        return None
    try:
        ik = j["PropertyTable"]["Properties"][0].get("InChIKey")
        return ik
    except Exception:
        return None


def get_iupac_by_inchikey(inchikey: str, timeout: float = 5.0) -> Optional[str]:
    """Fetch IUPACName from PubChem using InChIKey. Returns None if not found."""
    if not inchikey:
        return None
    q = urllib.parse.quote(inchikey, safe="")
    url = f"{BASE}/compound/inchikey/{q}/property/IUPACName/JSON"
    j = _rate_limited_get_json(url, timeout=timeout)
    if not j:
        return None
    try:
        return j["PropertyTable"]["Properties"][0].get("IUPACName")
    except Exception:
        return None


def get_iupac_by_smiles_using_inchikey(smiles: str, timeout: float = 5.0) -> Optional[str]:
    """Compute InChIKey for SMILES then fetch IUPACName from PubChem. """
    ik = get_inchikey_for_smiles(smiles, timeout=timeout)
    if not ik:
        return None
    return get_iupac_by_inchikey(ik, timeout=timeout)
