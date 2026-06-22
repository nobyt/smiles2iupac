"""IUPACResult: IUPAC名とその命名ステータスを保持するラッパークラス。

smiles_to_iupac() が返す文字列の代わりに使用でき、
  - 保留名 (retained) と優先IUPAC名 (PIN) の両方を保持
  - 複数の有効な PIN が存在する場合はリストとして保持
  - PubChem との照合結果 (verdict) を保持

後方互換性: str との == 比較は self.name で行う。
"""
from __future__ import annotations

from typing import List, Optional


class IUPACResult:
    """IUPAC名とその命名ステータスを保持するラッパー。

    Attributes:
        name: 本ライブラリが出力する名称（PIN または保留名）。
        iupac_pin: IUPAC 2013 優先IUPAC名 (PIN)。name と異なる場合は
                   name が保留名であることを意味する。None = 不明/未分類。
        retained: この化合物に有効な保留名のリスト。
        alternatives: その他の有効な PIN のリスト。
        verdict: PubChem 比較結果の分類。
            'pin'            — name が IUPAC 2013 PIN
            'retained'       — name が保留名 (iupac_pin に真の PIN)
            'pubchem_wrong'  — PubChem 名が IUPAC 2013 に従っていない
            'pubchem_retained' — PubChem が保留名を使用; name が系統的 PIN
            'both_pin'       — name と PubChem 名の両方が有効な PIN
            'tautomer'       — 互変異性体の命名差
            'needs_review'   — 自動分類不能
            None             — 未分類
        note: 分類の説明（日本語可）。
        iupac_source: 参照した IUPAC 2013 条項 (例: "P-16.3.4.1")。
    """

    __slots__ = ("name", "iupac_pin", "retained", "alternatives",
                 "verdict", "note", "iupac_source")

    def __init__(
        self,
        name: str,
        *,
        iupac_pin: Optional[str] = None,
        retained: Optional[List[str]] = None,
        alternatives: Optional[List[str]] = None,
        verdict: Optional[str] = None,
        note: Optional[str] = None,
        iupac_source: Optional[str] = None,
    ) -> None:
        self.name = name
        self.iupac_pin = iupac_pin
        self.retained = retained or []
        self.alternatives = alternatives or []
        self.verdict = verdict
        self.note = note
        self.iupac_source = iupac_source

    # --- string-compatible interface ---

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        parts = [repr(self.name)]
        if self.verdict:
            parts.append(f"verdict={self.verdict!r}")
        if self.iupac_pin and self.iupac_pin != self.name:
            parts.append(f"iupac_pin={self.iupac_pin!r}")
        if self.retained:
            parts.append(f"retained={self.retained!r}")
        if self.alternatives:
            parts.append(f"alternatives={self.alternatives!r}")
        return f"IUPACResult({', '.join(parts)})"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return self.name == other
        if isinstance(other, IUPACResult):
            return self.name == other.name
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self.name)

    # --- extended interface ---

    def all_names(self) -> List[str]:
        """全有効名を返す: 主名 + 代替 PIN + 保留名。重複除去・順序保持。"""
        seen: set[str] = {self.name}
        result = [self.name]
        for n in self.alternatives + self.retained:
            if n not in seen:
                seen.add(n)
                result.append(n)
        return result

    def is_retained(self) -> bool:
        """name が保留名 (retained) かどうか。"""
        return self.verdict == "retained"

    def is_pin(self) -> bool:
        """name が PIN かどうか (pubchem_wrong / pin / both_pin / pubchem_retained)。"""
        return self.verdict in ("pin", "pubchem_wrong", "both_pin", "pubchem_retained")
