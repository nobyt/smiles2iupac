"""pytest 設定。

- RDKit が未インストールの場合は全テストをスキップ。
- PubChem の検証とキャッシュ（tests/pubchem_cache.json）を行う。
- 差異は IUPAC 2013 規則で分類し、tests/pubchem_mismatch_report.txt に出力。

レポートの verdict 分類:
  [PUBCHEM_WRONG]    PubChem が IUPAC 2013 に従っていない (我々の名前が PIN)
  [OUR_RETAINED]     我々の名前が保留名 (PIN は PubChem 名)
  [PUBCHEM_RETAINED] PubChem が保留名を使用; 我々の名前が系統的 PIN
  [TAUTOMER]         互変異性体の命名差 (どちらも正当)
  [NEEDS_REVIEW]     自動分類不能 (手動確認が必要)
"""
from __future__ import annotations

import json
import os
import tempfile
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Any, Dict, List

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


def _classify_cache_entries(cache: Dict[str, Any]) -> bool:
    """キャッシュ内の未分類 mismatch エントリを IUPAC 2013 規則で分類する。

    verdict が未設定で pubchem_iupac と expected_in_test が揃っているエントリを対象とする。
    Returns True if any entry was updated.
    """
    try:
        from smiles2iupac.iupac_verdict import classify_mismatch
    except Exception:
        return False

    updated = False
    for entry in cache.values():
        if not isinstance(entry, dict):
            continue
        if entry.get("validated"):
            continue
        if entry.get("iupac_verdict"):
            continue  # 既に分類済み
        pub = entry.get("pubchem_iupac")
        our = entry.get("expected_in_test")
        if not pub or not our:
            continue
        result = classify_mismatch(str(our), str(pub))
        entry.update(result)
        updated = True
    return updated


def _write_mismatch_report(cache: Dict[str, Any]) -> None:
    """キャッシュ全体から verdict 付き mismatch レポートを生成・上書き。

    セクション:
      [PUBCHEM_WRONG]    — PubChem 名が IUPAC 2013 に非準拠 (我々の名前が正しい)
      [OUR_RETAINED]     — 我々の名前が保留名; PubChem 名が PIN
      [PUBCHEM_RETAINED] — PubChem が保留名を使用; 我々の名前が系統的 PIN
      [TAUTOMER]         — 互変異性体の命名差
      [NEEDS_REVIEW]     — 未分類
    """
    buckets: Dict[str, List[dict]] = defaultdict(list)
    for entry in cache.values():
        if not isinstance(entry, dict):
            continue
        if entry.get("validated"):
            continue
        if entry.get("not_found"):
            continue
        pub = entry.get("pubchem_iupac")
        our = entry.get("expected_in_test")
        if not pub or not our:
            continue
        verdict = entry.get("iupac_verdict", "needs_review")
        buckets[verdict].append(entry)

    section_order = [
        "pubchem_wrong",
        "our_retained",
        "pubchem_retained",
        "both_pin",
        "tautomer",
        "needs_review",
    ]
    section_labels = {
        "pubchem_wrong":    "[PUBCHEM_WRONG]    PubChem 名が IUPAC 2013 非準拠 — 我々の名前が正しい PIN",
        "our_retained":     "[OUR_RETAINED]     我々の名前が IUPAC 2013 保留名 — PubChem 名が優先 IUPAC 名 (PIN)",
        "pubchem_retained": "[PUBCHEM_RETAINED] PubChem が保留名を使用 — 我々の名前が系統的 PIN",
        "both_pin":         "[BOTH_PIN]         両名前とも有効な PIN",
        "tautomer":         "[TAUTOMER]         互変異性体の命名差 — どちらも正当",
        "needs_review":     "[NEEDS_REVIEW]     自動分類不能 — IUPAC 2013 による手動確認が必要",
    }

    total = sum(len(v) for v in buckets.values())
    lines: List[str] = []
    lines.append(f"PubChem × IUPAC 2013 命名比較レポート — {date.today()}")
    lines.append(f"対象化合物数: {total}")
    lines.append("")
    lines.append("凡例:")
    for v, label in section_labels.items():
        cnt = len(buckets.get(v, []))
        lines.append(f"  {label}: {cnt} 件")
    lines.append("")
    lines.append("=" * 72)

    for verdict in section_order:
        entries = buckets.get(verdict, [])
        if not entries:
            continue
        lines.append("")
        lines.append(f"--- {section_labels[verdict]} ({len(entries)} 件) ---")
        lines.append("")

        for e in entries:
            ik = e.get("inchikey", "")
            smiles = e.get("smiles", "")
            our = e.get("expected_in_test", "")
            pub = e.get("pubchem_iupac", "")
            note = e.get("iupac_note", "")
            source = e.get("iupac_source", "")
            iupac_pin = e.get("iupac_pin", "")

            lines.append(f"InChIKey: {ik}")
            lines.append(f"SMILES:   {smiles}")

            if verdict == "pubchem_wrong":
                lines.append(f"PIN (我々): {our}")
                lines.append(f"PubChem:    {pub}  ← 非準拠")
            elif verdict == "our_retained":
                lines.append(f"保留名 (我々): {our}")
                lines.append(f"PIN (PubChem): {pub}")
                if iupac_pin:
                    lines.append(f"IUPAC PIN:     {iupac_pin}")
            elif verdict == "pubchem_retained":
                lines.append(f"系統的 PIN (我々): {our}")
                lines.append(f"保留名 (PubChem): {pub}")
            elif verdict == "both_pin":
                lines.append(f"PIN (我々):    {our}")
                lines.append(f"PIN (PubChem): {pub}")
            elif verdict == "tautomer":
                lines.append(f"我々の名前:    {our}")
                lines.append(f"PubChem 名:    {pub}")
            else:  # needs_review
                lines.append(f"我々の名前: {our}")
                lines.append(f"PubChem 名: {pub}")
                lines.append("→ IUPAC 2013 規則による手動確認が必要")

            if note:
                lines.append(f"根拠: {note}")
            if source:
                lines.append(f"参照: {source}")
            lines.append("")

    try:
        with open(MISMATCH_PATH, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines))
            fh.write("\n")
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
        return

    try:
        from smiles2iupac.iupac_verdict import classify_mismatch
        _has_verdict = True
    except Exception:
        _has_verdict = False

    cache = _load_cache()
    updated = False

    # キャッシュ内の既存 mismatch エントリを一括分類 (初回のみ実行コストあり)
    if _has_verdict and _classify_cache_entries(cache):
        updated = True

    seen_inchikeys: set[str] = set()
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
            entry.pop("iupac_verdict", None)
            entry.pop("iupac_note", None)
            entry.pop("iupac_source", None)
            entry.pop("iupac_pin", None)
            updated = True
        else:
            entry["validated"] = False
            entry["expected_in_test"] = expected
            # 新規 mismatch を即座に分類
            if _has_verdict and not entry.get("iupac_verdict"):
                verdict_info = classify_mismatch(str(expected), pub_name)
                entry.update(verdict_info)
            updated = True
        cache[ik] = entry

    if updated:
        _save_cache(cache)
        _write_mismatch_report(cache)
