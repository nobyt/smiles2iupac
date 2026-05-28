"""
スピロ化合物・架橋二環式化合物（von Baeyer）の検出と命名。

IUPAC 2013 Blue Book:
  P-31.1.3.1  spiro[m.n]alkane
  P-31.1.3.2  bicyclo[l.m.n]alkane (von Baeyer)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .molecule_analyzer import MoleculeGraph

from .constants import CHAIN_PREFIX


def _ring_neighbors(graph: "MoleculeGraph") -> dict[int, list[int]]:
    """環原子インデックス → 環内隣接原子リスト を返す。"""
    from .molecule_analyzer import get_atom
    ring_set = {a.idx for a in graph.atoms if a.in_ring}
    result: dict[int, list[int]] = {}
    for idx in ring_set:
        result[idx] = [nb for nb in graph.adjacency[idx] if nb in ring_set]
    return result


def _connected_components(
    atom_set: set[int],
    ring_nbrs: dict[int, list[int]],
    exclude: set[int],
) -> list[set[int]]:
    """
    atom_set 内の原子を ring_nbrs で接続し、exclude を除いた
    連結成分のリストを返す。
    """
    remaining = atom_set - exclude
    visited: set[int] = set()
    components: list[set[int]] = []

    for start in remaining:
        if start in visited:
            continue
        component: set[int] = set()
        stack = [start]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            component.add(node)
            for nb in ring_nbrs.get(node, []):
                if nb not in visited and nb not in exclude:
                    stack.append(nb)
        components.append(component)

    return components


def _alkane_name(n: int) -> str:
    """炭素数 n のアルカン名を返す（例: 7 → heptane）。"""
    stem = CHAIN_PREFIX.get(n, f"C{n}")
    return f"{stem}ane"


def _try_spiro(graph: "MoleculeGraph") -> str | None:
    """
    スピロ化合物なら 'spiro[m.n]alkane' を返す。
    スピロ原子 = 環結合数 4 の原子がちょうど 1 つ存在する純アリ環系。
    """
    from .molecule_analyzer import get_atom

    ring_nbrs = _ring_neighbors(graph)
    ring_atoms = set(ring_nbrs.keys())

    # 純炭化水素・非芳香族のみ対応
    for idx in ring_atoms:
        a = get_atom(graph, idx)
        if a.symbol != "C" or a.is_aromatic:
            return None

    # ring 以外の重原子があれば置換基あり → スキップ
    heavy = [a for a in graph.atoms if a.symbol not in ("H",)]
    if len(heavy) != len(ring_atoms):
        return None

    deg4 = [idx for idx, nbs in ring_nbrs.items() if len(nbs) == 4]
    deg3 = [idx for idx, nbs in ring_nbrs.items() if len(nbs) == 3]

    if len(deg4) != 1 or deg3:
        return None

    spiro_atom = deg4[0]
    components = _connected_components(ring_atoms, ring_nbrs, exclude={spiro_atom})

    if len(components) != 2:
        return None

    sizes = sorted(len(c) for c in components)
    m, n = sizes[0], sizes[1]
    total = m + n + 1
    return f"spiro[{m}.{n}]{_alkane_name(total)}"


def _try_bicyclo(graph: "MoleculeGraph") -> str | None:
    """
    von Baeyer 式二環化合物なら 'bicyclo[l.m.n]alkane' を返す。
    架橋頭炭素 = 環結合数 3 の原子がちょうど 2 つ存在する純アリ環系。
    """
    from .molecule_analyzer import get_atom

    ring_nbrs = _ring_neighbors(graph)
    ring_atoms = set(ring_nbrs.keys())

    # 純炭化水素・非芳香族のみ対応
    for idx in ring_atoms:
        a = get_atom(graph, idx)
        if a.symbol != "C" or a.is_aromatic:
            return None

    # ring 以外の重原子があれば置換基あり → スキップ
    heavy = [a for a in graph.atoms if a.symbol not in ("H",)]
    if len(heavy) != len(ring_atoms):
        return None

    deg3 = [idx for idx, nbs in ring_nbrs.items() if len(nbs) == 3]
    deg4 = [idx for idx, nbs in ring_nbrs.items() if len(nbs) >= 4]

    if len(deg3) != 2 or deg4:
        return None

    b1, b2 = deg3
    components = _connected_components(ring_atoms, ring_nbrs, exclude={b1, b2})

    bridge_lengths = sorted((len(c) for c in components), reverse=True)

    # ゼロ橋: 架橋頭同士が直接結合している場合
    if b2 in ring_nbrs[b1]:
        bridge_lengths.append(0)
        bridge_lengths.sort(reverse=True)

    if len(bridge_lengths) != 3:
        return None

    l, m, n = bridge_lengths
    total = l + m + n + 2
    return f"bicyclo[{l}.{m}.{n}]{_alkane_name(total)}"


def name_polycyclic(graph: "MoleculeGraph") -> str | None:
    """
    スピロまたは架橋二環式化合物なら IUPAC 名を返す。
    対象外の場合は None を返す。
    """
    result = _try_spiro(graph)
    if result is not None:
        return result
    return _try_bicyclo(graph)
