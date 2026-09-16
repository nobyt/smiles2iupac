"""非環式骨格置換 ('a') 命名法 (Skeletal Replacement / 'a' Nomenclature, IUPAC 2013 P-15.4, P-55).

4個以上のヘテロ原子（O, S, N, Si, P等）を含む非環式鎖状化合物を
'a' 命名法（tetraoxadodecane 等）により Preferred IUPAC Name (PIN) として命名する。
"""

from __future__ import annotations

from typing import TYPE_CHECKING
from collections import defaultdict

if TYPE_CHECKING:
    from .molecule_analyzer import MoleculeGraph

from .constants import CHAIN_PREFIX, MULTIPLIER

_HET_PRIORITY: dict[str, int] = {
    "O": 0,
    "S": 1,
    "Se": 2,
    "Te": 3,
    "N": 4,
    "P": 5,
    "As": 6,
    "Si": 7,
    "B": 8,
}

_HET_NAME: dict[str, str] = {
    "O": "oxa",
    "S": "thia",
    "Se": "selena",
    "Te": "tellura",
    "N": "aza",
    "P": "phospha",
    "As": "arsa",
    "Si": "sila",
    "B": "bora",
}


def _name_acyclic_skeletal_replacement(graph: "MoleculeGraph", get_atom) -> str | None:
    """4個以上のヘテロ原子を含む最長非環式鎖に対して 'a' 命名法を適用する。"""
    # 環を含む場合は対象外（環は別ハンドラ）
    if any(a.in_ring for a in graph.atoms):
        return None

    # 重原子の集合
    heavy_atoms = [a.idx for a in graph.atoms if a.symbol != "H"]
    if len(heavy_atoms) < 5:
        return None

    # ヘテロ原子数を確認
    het_atoms = [idx for idx in heavy_atoms if get_atom(graph, idx).symbol in _HET_PRIORITY]
    if len(het_atoms) < 4:
        return None

    from .molecule_analyzer import get_bond_order

    # P-55 骨格置換ヘテロ原子の妥当性検証:
    # 1. すべてのヘテロ原子が単結合 (ether -O-, thioether -S-, amine -NH-/-NR-, etc.)
    # 2. 形式電荷なし (イリド、オニウム、アジド等の特殊種を除外)
    # 3. ヘテロ原子同士が直接結合していない (ペルオキシド -O-O-, 多硫化物 -S-S-, アゾ -N=N-, アジド等を除外)
    for h_idx in het_atoms:
        h_atom = get_atom(graph, h_idx)
        if h_atom.formal_charge != 0:
            return None
        for nb in graph.adjacency[h_idx]:
            if get_atom(graph, nb).symbol == "H":
                continue
            if get_bond_order(graph, h_idx, nb) != 1.0:
                return None
            if get_atom(graph, nb).symbol in _HET_PRIORITY:
                return None

    # 重原子グラフでの最長単純パスを探索
    best_path = _find_longest_acyclic_path(graph, heavy_atoms, get_atom)
    if not best_path:
        return None

    # パス上のヘテロ原子数
    path_hets = [idx for idx in best_path if get_atom(graph, idx).symbol in _HET_PRIORITY]
    if len(path_hets) < 4:
        return None

    # パス外に非水素置換基があるか確認 (純粋な直鎖 PEG / ポリアミン等を優先サポート)
    path_set = set(best_path)
    ext_heavy = [idx for idx in heavy_atoms if idx not in path_set]
    if ext_heavy:
        # 分岐置換基がある場合は現状スキップして通常命名法にフォールバック
        return None

    n = len(best_path)
    alkane_stem = CHAIN_PREFIX.get(n)
    if alkane_stem is None:
        return None
    alkane_name = f"{alkane_stem}ane"

    # 順方向と逆方向のロカント付けと優先度比較
    fwd_path = best_path
    rev_path = list(reversed(best_path))

    def _eval_path(path: list[int]) -> tuple:
        lmap = {idx: i + 1 for i, idx in enumerate(path)}
        by_sym: dict[str, list[int]] = defaultdict(list)
        for idx in path:
            sym = get_atom(graph, idx).symbol
            if sym in _HET_PRIORITY:
                by_sym[sym].append(lmap[idx])
        # 比較キー: (優先順位0のロカントリスト, 優先順位1のロカントリスト, ...)
        key = []
        for sym in sorted(_HET_PRIORITY.keys(), key=lambda s: _HET_PRIORITY[s]):
            locs = sorted(by_sym.get(sym, []))
            key.append(tuple(locs))
        return tuple(key), lmap

    key_fwd, lmap_fwd = _eval_path(fwd_path)
    key_rev, lmap_rev = _eval_path(rev_path)

    chosen_lmap = lmap_fwd if key_fwd <= key_rev else lmap_rev

    # 接頭辞の組み立て
    by_sym: dict[str, list[int]] = defaultdict(list)
    for idx in best_path:
        sym = get_atom(graph, idx).symbol
        if sym in _HET_PRIORITY:
            by_sym[sym].append(chosen_lmap[idx])

    prefix_parts: list[str] = []
    for sym in sorted(by_sym.keys(), key=lambda s: _HET_PRIORITY.get(s, 99)):
        locs = sorted(by_sym[sym])
        loc_str = ",".join(str(l) for l in locs)
        n_het = len(locs)
        mult = MULTIPLIER.get(n_het, f"{n_het}")
        het_name = _HET_NAME.get(sym, sym.lower())
        prefix_parts.append(f"{loc_str}-{mult}{het_name}")

    # 複数種のヘテロ原子接頭辞を連結するとき、後続の接頭辞の先頭ロカントを
    # ハイフンで区切る (例: "2,5,11-trioxa" + "8-aza" → "2,5,11-trioxa-8-aza")。
    # 単一種の場合はハイフンは付かない (例: "2,5,8,11-tetraoxa")。
    prefix_str = "-".join(prefix_parts)
    return f"{prefix_str}{alkane_name}"


def _find_longest_acyclic_path(graph: "MoleculeGraph", heavy_atoms: list[int], get_atom) -> list[int]:
    """非環式重原子グラフ内の最長単純パスを探索する。"""
    adj: dict[int, list[int]] = {}
    for idx in heavy_atoms:
        adj[idx] = [nb for nb in graph.adjacency[idx] if get_atom(graph, nb).symbol != "H"]

    leaves = [idx for idx in heavy_atoms if len(adj[idx]) <= 1]
    if not leaves:
        leaves = heavy_atoms

    best_path: list[int] = []

    for start in leaves:
        stack = [(start, [start], {start})]
        while stack:
            curr, path, visited = stack.pop()
            if len(path) > len(best_path):
                best_path = path
            for nb in adj.get(curr, []):
                if nb not in visited:
                    stack.append((nb, path + [nb], visited | {nb}))

    return best_path
