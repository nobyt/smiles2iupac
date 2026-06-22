"""pytest 設定。

- RDKit が未インストールの場合は全テストをスキップ。
- PubChem の検証とキャッシュ（tests/pubchem_cache.json）を行い、差異は tests/pubchem_mismatch_report.txt に追記します。
"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Dict, Any

import pytest


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "rdkit: marks tests that require RDKit"
    )


CACHE_PATH = Path(__file__).parent / "pubchem_cache.json"
MISMATCH_PATH = Path(__file__).parent / "pubchem_mismatch_report.txt"


def _load_cache() -> Dict[str, Any]:
    if not CACHE_PATH.exists():
        return {}
    try:
        with open(CACHE_PATH, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return {}


def _save_cache(cache: Dict[str, Any]) -> None:
    tmp = tempfile.NamedTemporaryFile(delete=False, dir=str(CACHE_PATH.parent))
    try:
        with open(tmp.name, "w", encoding="utf-8") as fh:
            json.dump(cache, fh, ensure_ascii=False, indent=2)
        os.replace(tmp.name, CACHE_PATH)
    finally:
        try:
            os.unlink(tmp.name)
        except Exception:
            pass


def pytest_collection_modifyitems(config, items):
    # First: existing behavior — skip all tests if RDKit not available
    try:
        import rdkit  # noqa: F401
    except ImportError:
        skip_rdkit = pytest.mark.skip(reason="RDKit not installed. Run: uv sync")
        for item in items:
            item.add_marker(skip_rdkit)
        return

    # Then: PubChem validation & caching for parametrized tests with (smiles, expected)
    pairs = []
    for item in items:
        callspec = getattr(item, "callspec", None)
        if not callspec:
            continue
        params = callspec.params
        if "smiles" in params and "expected" in params:
            smiles = params["smiles"]
            expected = params["expected"]
            pairs.append((smiles, expected))

    if not pairs:
        return

    try:
        from smiles2iupac.pubchem import (
            get_inchikey_for_smiles,
            get_iupac_by_inchikey,
        )
    except Exception:
        # pubchem utilities not available; skip validation
        return

    cache = _load_cache()
    mismatches = []
    updated = False

    seen_inchikeys = set()
    for smiles, expected in pairs:
        try:
            ik = get_inchikey_for_smiles(smiles)
        except Exception:
            ik = None
        if not ik:
            key = f"SMILES:{smiles}"
            entry = cache.get(key, {})
            if not entry.get("not_found"):
                entry.update({"smiles": smiles, "not_found": True})
                cache[key] = entry
                updated = True
            continue

        if ik in seen_inchikeys:
            continue
        seen_inchikeys.add(ik)

        entry = cache.get(ik, {})
        if entry.get("validated"):
            continue

        try:
            pub_name = get_iupac_by_inchikey(ik)
        except Exception:
            pub_name = None

        if not pub_name:
            entry.update({"inchikey": ik, "smiles": smiles, "not_found": True})
            cache[ik] = entry
            updated = True
            continue

        norm_pub = " ".join(pub_name.strip().lower().split())
        norm_expected = " ".join(str(expected).strip().lower().split())
        entry.update({"inchikey": ik, "smiles": smiles, "pubchem_iupac": pub_name})
        if norm_pub == norm_expected:
            entry["validated"] = True
            updated = True
        else:
            entry["validated"] = False
            entry["expected_in_test"] = expected
            updated = True
            mismatches.append({"inchikey": ik, "smiles": smiles, "expected": expected, "pubchem": pub_name})
        cache[ik] = entry

    if updated:
        _save_cache(cache)

    if mismatches:
        try:
            with open(MISMATCH_PATH, "a", encoding="utf-8") as fh:
                fh.write("=== PubChem mismatch report\n")
                for m in mismatches:
                    fh.write(f"InChIKey: {m['inchikey']}\n")
                    fh.write(f"SMILES: {m['smiles']}\n")
                    fh.write(f"expected (test): {m['expected']}\n")
                    fh.write(f"pubchem IUPAC: {m['pubchem']}\n")
                    fh.write("Refer to IUPAC 2013 recommendations (https://iupac.org) to reconcile differences.\n\n")
        except Exception:
            pass
