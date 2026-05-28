"""
環状化合物（シクロアルカン・芳香族・縮合多環）の検出と命名。

IUPAC 2013 Blue Book P-31.1.3 (cycloalkanes) / P-31.1.3.4 (benzene)
                    P-31.1.2 (naphthalene)
                    P-31.1.3.4 (anthracene, phenanthrene)

対応範囲:
  - シクロアルカン (cyclopropane ~ cyclodecane)
  - 置換シクロアルカン (methylcyclohexane, 1,2-dimethylcyclopentane ...)
  - ベンゼン (benzene) / 置換ベンゼン
  - ベンゼン官能基誘導体 (phenol, benzoic acid, benzaldehyde)
  - ナフタレン (naphthalene) / 置換ナフタレン
  - アントラセン (anthracene) / フェナントレン (phenanthrene) / 置換体
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .constants import CHAIN_PREFIX, HALOGEN_PREFIX, MULTIPLIER

if TYPE_CHECKING:
    from .molecule_analyzer import MoleculeGraph


@dataclass
class RingChain:
    """環状主鎖の表現。"""
    ring_atoms: list[int]       # ロカント順の環原子インデックス
    locant_map: dict[int, int]  # atom_idx -> locant (1始まり)
    is_aromatic: bool
    ring_size: int              # 実際の環原子数 (naphthalene=10, benzene=6, ...)
    base_name: str | None = None  # 保留名オーバーライド ("naphthalene" 等)
    double_bond_locants: list[int] = None  # 環内二重結合の最小ロカント（各結合の小さい方）

    def __post_init__(self) -> None:
        if self.double_bond_locants is None:
            self.double_bond_locants = []

    @property
    def length(self) -> int:
        return self.ring_size


# ─── 環検出 ─────────────────────────────────────────────────────────

_HETERO_SYMBOLS = frozenset({"N", "O", "S", "P", "Se", "Te"})


def has_ring(graph: "MoleculeGraph") -> bool:
    """分子に環が含まれるか判定する。"""
    return any(a.in_ring for a in graph.atoms if a.symbol == "C")


def has_hetero_ring(graph: "MoleculeGraph") -> bool:
    """分子にヘテロ原子を含む環があるか判定する。"""
    return any(a.in_ring and a.symbol in _HETERO_SYMBOLS for a in graph.atoms)


def find_hetero_rings(graph: "MoleculeGraph") -> list[list[int]]:
    """
    ヘテロ原子を含む環を隣接順リストで返す。
    ring_atom_sets (RDKit SSSR) を使用。
    """
    hetero_in_ring = {a.idx for a in graph.atoms
                      if a.in_ring and a.symbol in _HETERO_SYMBOLS}
    if not hetero_in_ring:
        return []

    if graph.ring_atom_sets:
        result: list[list[int]] = []
        seen: set[frozenset[int]] = set()
        for ring_tuple in graph.ring_atom_sets:
            if not any(idx in hetero_in_ring for idx in ring_tuple):
                continue
            fs = frozenset(ring_tuple)
            if fs in seen:
                continue
            seen.add(fs)
            ordered = _order_ring(list(ring_tuple), graph)
            if ordered:
                result.append(ordered)
        return result

    return []


def has_aromatic(graph: "MoleculeGraph") -> bool:
    """分子に芳香族環が含まれるか判定する。"""
    return any(a.is_aromatic for a in graph.atoms)


def find_rings(graph: "MoleculeGraph") -> list[list[int]]:
    """
    RDKit の SSSR (最小環系) から炭素環のリストを返す。

    graph.ring_atom_sets に RDKit の AtomRings() が格納されている場合はそれを使用。
    なければ DFS フォールバック。

    Returns:
        炭素原子のみで構成される各環を隣接順リストで返す。
    """
    c_set = {a.idx for a in graph.atoms if a.symbol == "C" and a.in_ring}
    if not c_set:
        return []

    # RDKit SSSR を優先使用（縮合環系で正しい最小環を返す）
    if graph.ring_atom_sets:
        result: list[list[int]] = []
        seen: set[frozenset[int]] = set()
        for ring_tuple in graph.ring_atom_sets:
            # 炭素のみの環のみ（ヘテロ環は除外）
            if any(idx not in c_set for idx in ring_tuple):
                continue
            fs = frozenset(ring_tuple)
            if fs in seen:
                continue
            seen.add(fs)
            ordered = _order_ring(list(ring_tuple), graph)
            if ordered:
                result.append(ordered)
        return result

    # ─── DFS フォールバック（ring_atom_sets が空の場合）────────────
    c_adj: dict[int, list[int]] = {
        c: [nb for nb in graph.adjacency[c] if nb in c_set]
        for c in c_set
    }
    found_rings: list[frozenset[int]] = []
    visited: set[int] = set()
    parent: dict[int, int] = {}

    def dfs(u: int, path: list[int]) -> None:
        visited.add(u)
        path.append(u)
        for v in c_adj[u]:
            if v not in visited:
                parent[v] = u
                dfs(v, path)
            elif v != parent.get(u) and v in path:
                idx = path.index(v)
                ring = frozenset(path[idx:])
                if ring not in found_rings:
                    found_rings.append(ring)
        path.pop()

    for start in c_set:
        if start not in visited:
            dfs(start, [])

    result_dfs: list[list[int]] = []
    for ring_set in found_rings:
        ordered = _order_ring(list(ring_set), graph)
        if ordered:
            result_dfs.append(ordered)
    return result_dfs


def _order_ring(atoms: list[int], graph: MoleculeGraph) -> list[int]:
    """
    環の原子を隣接順に並べる。
    DFS で1周するリストを返す。
    """
    if not atoms:
        return []
    atom_set = set(atoms)
    start = atoms[0]
    ordered = [start]
    visited = {start}

    current = start
    for _ in range(len(atoms) - 1):
        found_next = False
        for nb in graph.adjacency[current]:
            if nb in atom_set and nb not in visited:
                ordered.append(nb)
                visited.add(nb)
                current = nb
                found_next = True
                break
        if not found_next:
            break

    return ordered if len(ordered) == len(atoms) else atoms


# ─── 主鎖選択 ────────────────────────────────────────────────────────

def _find_principal_group_anchor(
    ring: list[int],
    ring_set: set[int],
    principal_grp_atoms: list[int],
    graph: "MoleculeGraph",
) -> int | None:
    """
    主官能基のアンカーとなる環原子を返す。

    - 環内型 (cyclohexanol, phenol): principal_grp_atoms 中に環原子がある
    - 環外型 (benzoic acid, benzaldehyde): 環原子の隣に principal_grp_atom がある
    """
    if not principal_grp_atoms:
        return None

    pgrp_set = set(principal_grp_atoms)

    # 環内型: 環原子自身が主官能基の一部
    for c_idx in ring:
        if c_idx in pgrp_set:
            return c_idx

    # 環外型: 環原子に隣接する主官能基原子（環外）
    for c_idx in ring:
        for nb in graph.adjacency[c_idx]:
            if nb in pgrp_set and nb not in ring_set:
                return c_idx

    return None


# ─── 4環以上の縮合多環芳香族 (PAH) フィンガープリント照合 ────────────

def _ring_graph_degrees(rings: list[list[int]]) -> tuple[int, ...]:
    """SSSR の環グラフ次数列（ソート済み）を返す。"""
    n = len(rings)
    ring_sets = [set(r) for r in rings]
    deg = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if len(ring_sets[i] & ring_sets[j]) >= 2:
                deg[i] += 1
                deg[j] += 1
    return tuple(sorted(deg))


# フィンガープリント: (n_rings, n_aromatic_C, ring_graph_degrees) → 保留名
_PAH_RETAINED: dict[tuple[int, int, tuple[int, ...]], str] = {
    (4, 16, (2, 2, 3, 3)): "pyrene",
    (4, 18, (1, 1, 1, 3)): "triphenylene",
    (4, 18, (1, 1, 2, 2)): "chrysene",
    (5, 20, (2, 2, 2, 2, 4)): "perylene",
}


def _try_pah_retained(
    arom_rings: list[list[int]],
    graph: "MoleculeGraph",
) -> "RingChain | None":
    """
    4環以上の縮合芳香族保留名を照合する。

    フィンガープリント (n_rings, n_aromatic_C, ring_graph_degrees) で同定。
    一致しない場合は None を返す。
    """
    n_rings = len(arom_rings)
    if n_rings < 4:
        return None

    from .molecule_analyzer import get_atom
    n_c = sum(
        1 for a in graph.atoms
        if a.symbol == "C" and a.is_aromatic and a.in_ring
    )
    deg = _ring_graph_degrees(arom_rings)
    key = (n_rings, n_c, deg)
    name = _PAH_RETAINED.get(key)
    if name is None:
        return None

    # 全芳香族 C をロカントマップ (ロカント値は名前参照用 dummy)
    all_arom_c = [a.idx for a in graph.atoms
                  if a.symbol == "C" and a.is_aromatic and a.in_ring]
    locant_map = {idx: i + 1 for i, idx in enumerate(all_arom_c)}

    return RingChain(
        ring_atoms=all_arom_c,
        locant_map=locant_map,
        is_aromatic=True,
        ring_size=n_c,
        base_name=name,
    )


# ─── 三環縮合芳香族 (anthracene / phenanthrene) 検出と命名 ───────────

def _find_triple_fused_aromatic(
    rings: list[list[int]],
    graph: "MoleculeGraph",
) -> "tuple[list[int], list[int], list[int], set[int], set[int]] | None":
    """
    3 つの 6 員芳香族環が 14 原子で縮合した三環系を検出する。

    Returns:
        (middle_ring, outer1, outer2, shared1, shared2) または None
        - middle_ring: 中央環の原子リスト
        - outer1, outer2: 外側 2 環の原子リスト
        - shared1: middle_ring と outer1 の共有原子集合 (橋頭 2 個)
        - shared2: middle_ring と outer2 の共有原子集合 (橋頭 2 個)
    """
    from itertools import combinations
    from .molecule_analyzer import get_atom

    def is_aromatic(ring: list[int]) -> bool:
        return all(get_atom(graph, idx).is_aromatic for idx in ring)

    arom_six = [r for r in rings if len(r) == 6 and is_aromatic(r)]
    if len(arom_six) < 3:
        return None

    for r0, r1, r2 in combinations(arom_six, 3):
        union = set(r0) | set(r1) | set(r2)
        if len(union) != 14:
            continue

        s01 = set(r0) & set(r1)
        s12 = set(r1) & set(r2)
        s02 = set(r0) & set(r2)

        # 中央環を特定: 2 つのペアに属する環が中央
        if len(s01) == 2 and len(s12) == 2 and len(s02) == 0:
            return list(r1), list(r0), list(r2), s01, s12
        if len(s01) == 2 and len(s02) == 2 and len(s12) == 0:
            return list(r0), list(r1), list(r2), s01, s02
        if len(s12) == 2 and len(s02) == 2 and len(s01) == 0:
            return list(r2), list(r0), list(r1), s12, s02

    return None


def _is_linear_tricyclic(
    middle_ring: list[int],
    shared1: "set[int]",
    shared2: "set[int]",
    graph: "MoleculeGraph",
) -> bool:
    """
    三環縮合が直線型 (anthracene) か角型 (phenanthrene) かを判定する。

    中央環内で 2 つの共有結合の位置関係を調べる。
    - para (距離 ≈ 3) → True  → anthracene (直線)
    - meta (距離 ≈ 2) → False → phenanthrene (角)
    """
    ordered = _order_ring(middle_ring, graph)
    n = len(ordered)  # 6

    sh1_pos = [ordered.index(a) for a in shared1 if a in ordered]
    sh2_pos = [ordered.index(a) for a in shared2 if a in ordered]

    if len(sh1_pos) < 2 or len(sh2_pos) < 2:
        return True  # デフォルト: anthracene

    mid1 = (sh1_pos[0] + sh1_pos[1]) / 2.0
    mid2 = (sh2_pos[0] + sh2_pos[1]) / 2.0

    dist = abs(mid1 - mid2)
    if dist > n / 2:
        dist = n - dist

    # anthracene: 3, phenanthrene: 2
    return dist >= 2.5


def _build_peripheral_order(
    peripheral: list[int],
    bridgeheads: "set[int]",
    all_set: "set[int]",
    graph: "MoleculeGraph",
) -> list[int]:
    """
    橋頭原子を経由しながら 10 個の周辺原子を環状順序に並べる。
    """
    periph_set = set(peripheral)
    bh_set = bridgeheads

    # 各周辺原子の「到達可能な隣接周辺原子」を構築（橋頭経由含む）
    periph_adj: dict[int, list[int]] = {a: [] for a in peripheral}
    for a in peripheral:
        for nb in graph.adjacency[a]:
            if nb in periph_set:
                if nb not in periph_adj[a]:
                    periph_adj[a].append(nb)
            elif nb in bh_set:
                for nb2 in graph.adjacency[nb]:
                    if nb2 in periph_set and nb2 != a and nb2 not in periph_adj[a]:
                        periph_adj[a].append(nb2)

    # DFS で 1 周する順序を探索
    best: list[int] = []
    for start in peripheral:
        order: list[int] = [start]
        visited: set[int] = {start}
        current = start
        for _ in range(9):
            found = False
            for nb in periph_adj[current]:
                if nb not in visited:
                    order.append(nb)
                    visited.add(nb)
                    current = nb
                    found = True
                    break
            if not found:
                break
        if len(order) == 10:
            return order

    return peripheral  # フォールバック


def _assign_tricyclic_locants(
    middle_ring: list[int],
    outer1: list[int],
    outer2: list[int],
    shared1: "set[int]",
    shared2: "set[int]",
    graph: "MoleculeGraph",
    pgrp_atoms: list[int],
    base_name: str,
) -> "RingChain":
    """
    アントラセン / フェナントレンの IUPAC ロカントを周辺 10 原子に割り当てる。

    anthracene ルール:
      - meso 位 (中央環の非橋頭原子 2 個) → locant 9, 10
      - 外環周辺原子 8 個 → locant 1–8
      - 置換基のロカント集合が最小になる開始・方向を選択

    phenanthrene ルール:
      - 周辺 10 原子に 1–10 を最小ロカント集合で割り当て
    """
    from .molecule_analyzer import get_atom

    bridgeheads = shared1 | shared2
    all_set = set(middle_ring) | set(outer1) | set(outer2)
    peripheral = [a for a in all_set if a not in bridgeheads]  # 10 個

    # meso 原子: 中央環の中で橋頭でない原子
    meso_set = {a for a in middle_ring if a not in bridgeheads}  # 2 個

    # 周辺原子の環状順序
    periph_order = _build_peripheral_order(peripheral, bridgeheads, all_set, graph)

    # 置換基を持つ周辺原子を収集
    pgrp_set = set(pgrp_atoms)
    sub_atoms: set[int] = set()
    for a in peripheral:
        for nb in graph.adjacency[a]:
            nb_atom = get_atom(graph, nb)
            if nb_atom.symbol == "H":
                continue
            if nb in all_set:
                continue
            sub_atoms.add(a)

    # 主官能基アンカーの特定
    anchor: "int | None" = None
    if pgrp_atoms:
        for a in peripheral:
            if a in pgrp_set:
                anchor = a
                break
        if anchor is None:
            for a in peripheral:
                for nb in graph.adjacency[a]:
                    if nb in pgrp_set:
                        anchor = a
                        break
                if anchor is not None:
                    break

    n = len(periph_order)
    best_map: "dict[int, int] | None" = None
    best_locs: list[int] = [11] * 100

    # アルファ原子: 橋頭に隣接する周辺原子 (IUPAC 位置 1 の候補)
    alpha_atoms = {
        a for a in peripheral
        if any(nb in bridgeheads for nb in graph.adjacency[a])
    }

    for start in range(n):
        for direction in (1, -1):
            arrangement = [periph_order[(start + direction * i) % n] for i in range(n)]

            # アンカー制約: 主官能基がある場合はそれを C1 に固定
            if anchor is not None and arrangement[0] != anchor:
                continue

            if base_name == "anthracene":
                # IUPAC 規則: アンカーがない場合、位置 1 はアルファ位から始める
                if anchor is None and arrangement[0] not in alpha_atoms:
                    continue

                # meso 原子が arrangement を 4+4 に分割するか検証
                meso_pos = [i for i, a in enumerate(arrangement) if a in meso_set]
                if len(meso_pos) != 2:
                    continue
                p1, p2 = sorted(meso_pos)
                seg1 = p2 - p1 - 1
                seg2 = (n - p2) + p1 - 1
                if {seg1, seg2} != {4, 4}:
                    continue

                # ロカント割り当て: 外環 → 1–8, meso → 9–10
                loc: dict[int, int] = {}
                outer_count = 0
                meso_count = 9
                for a in arrangement:
                    if a in meso_set:
                        loc[a] = meso_count
                        meso_count += 1
                    else:
                        outer_count += 1
                        loc[a] = outer_count
            else:
                # phenanthrene: 周辺に 1–10 を連番
                # アルファ位から始める (IUPAC 慣習)
                if anchor is None and arrangement[0] not in alpha_atoms:
                    continue
                loc = {a: i + 1 for i, a in enumerate(arrangement)}

            locs = sorted(loc[a] for a in sub_atoms if a in loc)
            if locs < best_locs:
                best_locs = locs
                best_map = dict(loc)

    if best_map is None:
        # フォールバック: アルファ位があればそこから始める
        start_atom = next((a for a in periph_order if a in alpha_atoms), periph_order[0])
        offset = periph_order.index(start_atom)
        fallback_order = [periph_order[(offset + i) % n] for i in range(n)]
        best_map = {a: i + 1 for i, a in enumerate(fallback_order)}

    return RingChain(
        ring_atoms=periph_order,          # 10 個の周辺原子のみ
        locant_map=best_map,
        is_aromatic=True,
        ring_size=14,                     # 実際の環炭素数
        base_name=base_name,
    )


# ─── ナフタレン検出と命名 ────────────────────────────────────────────

def _find_fused_pair(rings: list[list[int]]) -> tuple[list[int], list[int]] | None:
    """
    2つの6員環が1結合（2原子）を共有するペアを返す。
    ナフタレン検出に使用。
    """
    six_rings = [r for r in rings if len(r) == 6]
    for i in range(len(six_rings)):
        for j in range(i + 1, len(six_rings)):
            shared = set(six_rings[i]) & set(six_rings[j])
            if len(shared) == 2:
                return six_rings[i], six_rings[j]
    return None


def _try_biphenyl(
    arom_six: list[list[int]],
    graph: "MoleculeGraph",
    pgrp_atoms: list[int],
) -> "RingChain | None":
    """
    2 つの 6 員芳香族環が非環内の単結合で連結 → biphenyl 系。

    片方の環を主環として locant を振り、もう片方を "phenyl" 置換基として扱う。
    主官能基なし (alkane) の場合のみ呼ぶ。
    """
    from .molecule_analyzer import get_atom as _ga

    _ring_sets = graph.ring_atom_sets or []

    # 2つの benzene 環が非縮合結合で連結する edge を探す
    for i in range(len(arom_six)):
        for j in range(i + 1, len(arom_six)):
            r1_set = set(arom_six[i])
            r2_set = set(arom_six[j])
            if r1_set & r2_set:
                continue  # 共有原子があれば縮合 → biphenyl ではない
            # r1 → r2 の inter-ring 単結合を探す
            bridge: list[tuple[int, int]] = []
            for a in arom_six[i]:
                for b in graph.adjacency[a]:
                    if b in r2_set:
                        bridge.append((a, b))
            if len(bridge) != 1:
                continue  # 単一結合のみ対象
            a_bridge, b_bridge = bridge[0]

            # 置換基が多い方の環を主環に選ぶ (置換ビフェニルで正しい locant を得るため)
            def _ext_sub_count(ring_atoms: list[int], skip_nb: int) -> int:
                r_set = set(ring_atoms)
                return sum(
                    1
                    for a in ring_atoms
                    for nb in graph.adjacency[a]
                    if nb not in r_set and nb != skip_nb and _ga(graph, nb).symbol != "H"
                )
            subs_i = _ext_sub_count(arom_six[i], b_bridge)
            subs_j = _ext_sub_count(arom_six[j], a_bridge)
            ring1 = arom_six[i] if subs_i >= subs_j else arom_six[j]
            return _assign_ring_locants(graph, ring1, True, "alkane", pgrp_atoms)

    return None


def _assign_naphthalene_locants(
    ring1: list[int],
    ring2: list[int],
    graph: "MoleculeGraph",
    principal_grp_atoms: list[int],
) -> RingChain:
    """
    ナフタレンの周辺原子（8個）にロカント 1〜8 を割り当てる。

    IUPAC ナフタレン番号付け:
      橋頭原子を 4a/8a として、周辺を 1→2→3→4→[4a]→5→6→7→8→[8a] と周回。
      置換基ロカント集合が最小になる開始点・方向を選ぶ。
      主官能基アンカーがある場合はそこを C1 に固定。
    """
    shared = set(ring1) & set(ring2)          # bridgehead atoms (4a, 8a)
    all_atoms = set(ring1) | set(ring2)       # 10 atoms
    peripheral = [a for a in all_atoms if a not in shared]  # 8 atoms

    # 周辺原子の隣接順序付け（bridgehead を経由して環を周回）
    # 開始点から環を辿り、bridgehead は通過して番号なしとする
    def traverse_peripheral(start: int) -> list[int] | None:
        """start から橋頭を経由して8個の周辺原子を順序付けする。"""
        order = [start]
        visited = {start}
        current = start
        for _ in range(9):  # 最大9ステップ（10原子周回）
            found = False
            for nb in graph.adjacency[current]:
                if nb not in all_atoms or nb in visited:
                    continue
                order.append(nb)
                visited.add(nb)
                current = nb
                found = True
                break
            if not found:
                break
        # 10原子全部辿れたか
        if len(order) == 10 and set(order) == all_atoms:
            return order
        return None

    # 全周辺原子から開始を試みる
    ring_order: list[int] | None = None
    for start in peripheral:
        o = traverse_peripheral(start)
        if o is not None:
            ring_order = o
            break

    if ring_order is None:
        # フォールバック: ring1 + ring2 合体
        ring_order = ring1 + [a for a in ring2 if a not in set(ring1)]

    # 周辺原子だけを抽出（bridgehead をスキップした位置にロカントを振る）
    # ring_order は [p1, p2, p3, p4, BH1, p5, p6, p7, p8, BH2] のような順
    # bridgehead は 4a / 8a に相当し、ロカントなし
    def peripheral_sequence(order: list[int]) -> list[int]:
        """橋頭を除いた周辺原子の順序リスト（8個）を返す。"""
        return [a for a in order if a not in shared]

    pgrp_set = set(principal_grp_atoms)
    peripheral_list = peripheral_sequence(ring_order)  # 8個の周辺原子

    # ─── Alpha 原子（橋頭原子に隣接する周辺原子）を特定 ──────────────
    # IUPAC ナフタレン番号付けでは位置 1 は常に alpha 位 (橋頭隣接)。
    # この制約により 1-chloro (alpha) と 2-chloro (beta) が区別される。
    alpha_atoms: set[int] = {
        a for a in peripheral_list
        if any(nb in shared for nb in graph.adjacency[a])
    }
    alpha_indices = [i for i, a in enumerate(peripheral_list) if a in alpha_atoms]

    # ─── 主官能基アンカーの特定 ──────────────────────────────────────
    anchor: int | None = None
    if principal_grp_atoms:
        for c_idx in peripheral_list:
            if c_idx in pgrp_set:
                anchor = c_idx
                break
        if anchor is None:
            for c_idx in peripheral_list:
                for nb in graph.adjacency[c_idx]:
                    if nb in pgrp_set and nb not in all_atoms:
                        anchor = c_idx
                        break
                if anchor is not None:
                    break

    # ─── 置換基を持つ周辺原子 ─────────────────────────────────────────
    sub_atoms: set[int] = set()
    for a in peripheral_list:
        for nb in graph.adjacency[a]:
            if nb not in all_atoms:
                from .molecule_analyzer import get_atom
                if get_atom(graph, nb).symbol != "H":
                    sub_atoms.add(a)

    # ─── 開始インデックスの決定 ───────────────────────────────────────
    # 主官能基アンカーがある場合: アンカーが alpha ならそこのみ試す。
    #   beta アンカーの場合: アンカーに隣接する alpha から始める。
    # なければ: 全 alpha 位を試す。
    if anchor is not None and anchor in alpha_atoms:
        # アンカーが alpha: そこを C1 に固定
        start_positions = [peripheral_list.index(anchor)]
    elif anchor is not None:
        # アンカーが beta: 隣接 alpha の中で最小ロカントになる方向を選ぶ
        start_positions = alpha_indices
    else:
        start_positions = alpha_indices if alpha_indices else list(range(8))

    best_arrangement: list[int] | None = None
    best_locs: list[int] = [9] * 100

    for si in start_positions:
        for direction in (1, -1):
            arrangement = [peripheral_list[(si + direction * i) % 8] for i in range(8)]
            locs = sorted(i + 1 for i, a in enumerate(arrangement) if a in sub_atoms)
            if locs < best_locs:
                best_locs = locs
                best_arrangement = arrangement

    if best_arrangement is None:
        best_arrangement = peripheral_list

    locant_map = {atom: i + 1 for i, atom in enumerate(best_arrangement)}
    return RingChain(
        ring_atoms=best_arrangement,
        locant_map=locant_map,
        is_aromatic=True,
        ring_size=10,
        base_name="naphthalene",
    )


def find_principal_ring(
    graph: "MoleculeGraph",
    rings: list[list[int]],
    principal_grp_type: str,
    principal_grp_atoms: list[int] | None = None,
) -> RingChain:
    """
    主となる環を選択してロカントを割り当てる。

    選択規則（IUPAC P-44 環への適用）:
      1. ナフタレン（6+6 縮合芳香族） → _assign_naphthalene_locants
      2. 最大単環を選択
      3. 置換基数が最多 / ロカント集合が最小
    """
    from .molecule_analyzer import get_atom  # 遅延インポート

    pgrp_atoms = principal_grp_atoms or []

    if not rings:
        raise ValueError("No rings found in molecule.")

    # 芳香族フラグ
    def ring_is_aromatic(ring: list[int]) -> bool:
        return all(get_atom(graph, idx).is_aromatic for idx in ring)

    # ─── 4環以上 PAH 保留名 ───────────────────────────────────────
    arom_six = [r for r in rings if len(r) == 6 and ring_is_aromatic(r)]
    if len(arom_six) >= 4:
        pah = _try_pah_retained(arom_six, graph)
        if pah is not None:
            return pah

    # ─── 三環縮合芳香族判定 (anthracene / phenanthrene) ──────────
    # ナフタレン判定より先に行う（三環はナフタレンのペアも含むため）
    if len(arom_six) >= 3:
        triple = _find_triple_fused_aromatic(arom_six, graph)
        if triple is not None:
            middle, outer1, outer2, shared1, shared2 = triple
            is_linear = _is_linear_tricyclic(middle, shared1, shared2, graph)
            base = "anthracene" if is_linear else "phenanthrene"
            return _assign_tricyclic_locants(
                middle, outer1, outer2, shared1, shared2,
                graph, pgrp_atoms, base,
            )

    # ─── ナフタレン判定 ───────────────────────────────────────────
    # 2つの6員芳香族環が1結合を共有 → ナフタレン系
    fused_pair = _find_fused_pair(arom_six) if len(arom_six) >= 2 else None
    if fused_pair is not None:
        return _assign_naphthalene_locants(fused_pair[0], fused_pair[1],
                                           graph, pgrp_atoms)

    # ─── ビフェニル判定 (Phase 37) ───────────────────────────────
    # 2つの6員芳香族環が縮合せず単結合で連結 → biphenyl 系
    if len(arom_six) >= 2 and principal_grp_type == "alkane":
        bp = _try_biphenyl(arom_six, graph, pgrp_atoms)
        if bp is not None:
            return bp

    # ─── 単環 ────────────────────────────────────────────────────
    if len(rings) == 1:
        ring = rings[0]
        is_arom = ring_is_aromatic(ring)
        return _assign_ring_locants(graph, ring, is_arom, principal_grp_type, pgrp_atoms)

    # 複数環: 最大環を選択
    best = max(rings, key=lambda r: len(r))
    is_arom = ring_is_aromatic(best)
    return _assign_ring_locants(graph, best, is_arom, principal_grp_type, pgrp_atoms)


def _alpha_sub_key(
    arrangement: list[int],
    ring_set: set[int],
    pgrp_atoms_set: set[int],
    graph: "MoleculeGraph",
) -> list[tuple[int, str]]:
    """
    Phase 156: アルファベット順タイブレーク用の置換基キーを生成する。

    各環原子の置換基を (ロカント, 置換基名) のリストとして返す。
    ロカント集合が等しい場合に、辞書順で小さい方の配置を選ぶ。
    """
    from .substituent import name_substituent
    from .molecule_analyzer import get_atom as _ga

    excl = ring_set | pgrp_atoms_set
    result: list[tuple[int, str]] = []
    for i, atom_idx in enumerate(arrangement):
        locant = i + 1
        for nb in graph.adjacency[atom_idx]:
            if _ga(graph, nb).symbol == "H" or nb in excl:
                continue
            sname = name_substituent(graph, nb, excl)
            if sname:
                result.append((locant, sname))
    return sorted(result)


def _assign_ring_locants(
    graph: "MoleculeGraph",
    ring: list[int],
    is_aromatic: bool,
    principal_grp_type: str,
    principal_grp_atoms: list[int] | None = None,
) -> RingChain:
    """
    環原子にロカントを割り当てる。

    規則:
      1. principal characteristic group のアンカー炭素を C1 にする
      2. 環内二重結合のロカント集合が最小になる方向を優先
      3. その後、置換基のロカント集合が最小になる方向を選ぶ
    """
    from .molecule_analyzer import get_atom, get_bond_order  # 遅延インポート

    pgrp_atoms = principal_grp_atoms or []
    n = len(ring)
    ring_set = set(ring)

    # 環内二重結合の検出: (a, b) のペアで a < b
    db_pairs: list[tuple[int, int]] = []
    db_atoms: set[int] = set()
    for c_idx in ring:
        for nb in graph.adjacency[c_idx]:
            if nb in ring_set and nb > c_idx:
                if get_bond_order(graph, c_idx, nb) == 2.0:
                    db_pairs.append((c_idx, nb))
                    db_atoms.add(c_idx)
                    db_atoms.add(nb)

    # 置換基を持つ炭素のインデックスを収集
    substituent_atoms: set[int] = set()
    for c_idx in ring:
        for nb in graph.adjacency[c_idx]:
            nb_atom = get_atom(graph, nb)
            if nb not in ring_set and nb_atom.symbol not in ("H",):
                substituent_atoms.add(c_idx)

    # 主官能基アンカー: C1 に固定する環原子
    anchor = _find_principal_group_anchor(ring, ring_set, pgrp_atoms, graph)

    has_content = bool(substituent_atoms or db_pairs or anchor)
    if not has_content:
        # 置換基なし・二重結合なし・アンカーなし: そのままのリスト順
        locant_map = {atom: i + 1 for i, atom in enumerate(ring)}
        return RingChain(ring_atoms=ring, locant_map=locant_map,
                        is_aromatic=is_aromatic, ring_size=n)

    # アンカーがある場合はそこだけを開始点とする。なければ全開始点を試す。
    if anchor is not None:
        start_indices = [ring.index(anchor)]
    else:
        start_indices = list(range(n))

    best_arrangement: list[int] | None = None
    best_key: tuple[list[int], list[int]] = ([n + 1] * 100, [n + 1] * 100)
    best_alpha_key: list[tuple[int, str]] = []

    pgrp_atoms_set = set(pgrp_atoms)

    for start_idx in start_indices:
        for direction in (1, -1):
            arrangement = _arrange_ring(ring, start_idx, direction)
            atom_to_loc = {atom: i + 1 for i, atom in enumerate(arrangement)}

            # 二重結合ロカント: ラップアラウンド結合を正しく処理
            db_locs = sorted(
                _db_locant_for_ring(atom_to_loc[a], atom_to_loc[b], n)
                for a, b in db_pairs
            )
            sub_locs = _ring_substituent_locants(arrangement, substituent_atoms)

            # 比較キー: (二重結合ロカント優先, 置換基ロカント)
            key = (db_locs, sub_locs)
            if _locant_key_less(key, best_key):
                best_key = key
                best_arrangement = arrangement
                best_alpha_key = []  # リセット（次の同順位で再計算）
            elif key == best_key and best_arrangement is not None:
                # Phase 156: 同一ロカント集合のとき置換基名のアルファベット順で決定
                if not best_alpha_key:
                    best_alpha_key = _alpha_sub_key(
                        best_arrangement, ring_set, pgrp_atoms_set, graph)
                cur_alpha = _alpha_sub_key(
                    arrangement, ring_set, pgrp_atoms_set, graph)
                if cur_alpha < best_alpha_key:
                    best_alpha_key = cur_alpha
                    best_arrangement = arrangement

    if best_arrangement is None:
        best_arrangement = ring

    locant_map = {atom: i + 1 for i, atom in enumerate(best_arrangement)}

    # 最終的な二重結合ロカントを計算
    final_db_locs: list[int] = []
    if db_pairs:
        final_db_locs = sorted(
            _db_locant_for_ring(locant_map[a], locant_map[b], n)
            for a, b in db_pairs
        )

    return RingChain(
        ring_atoms=best_arrangement,
        locant_map=locant_map,
        is_aromatic=is_aromatic,
        ring_size=n,
        double_bond_locants=final_db_locs,
    )


def _locant_key_less(
    key: tuple[list[int], list[int]],
    other: tuple[list[int], list[int]],
) -> bool:
    """(db_locs, sub_locs) の辞書比較で key < other かどうか返す。"""
    a_db, a_sub = key
    b_db, b_sub = other
    if a_db != b_db:
        return a_db < b_db
    return a_sub < b_sub


def _db_locant_for_ring(p: int, q: int, n: int) -> int:
    """
    環内二重結合のロカントを返す。

    通常: low = min(p, q)
    ラップアラウンド結合 (高番号-低番号 = n-1) の場合: locant = high
    例: n=6, (1,6) → wrap-around → 6  (not 1)
        n=6, (1,2) → 1
        n=6, (5,6) → 5
    """
    low, high = min(p, q), max(p, q)
    if high - low == n - 1:
        return high
    return low


def _arrange_ring(ring: list[int], start_idx: int, direction: int) -> list[int]:
    """
    環原子リストを start_idx から始めて direction (+1/-1) 順に並べる。
    """
    n = len(ring)
    return [ring[(start_idx + direction * i) % n] for i in range(n)]


def _ring_substituent_locants(arrangement: list[int], sub_atoms: set[int]) -> list[int]:
    """
    並べた環原子リストで、置換基を持つ原子のロカントを返す（昇順）。
    """
    locs = []
    for i, atom in enumerate(arrangement):
        if atom in sub_atoms:
            locs.append(i + 1)
    return sorted(locs)


# ─── 置換基収集（環用）──────────────────────────────────────────────

def collect_ring_substituents(
    graph: MoleculeGraph,
    ring_chain: RingChain,
    principal_grp_atoms: list[int],
) -> list[tuple[int, str]]:
    """
    環炭素に付く置換基を収集する。

    Returns:
        [(locant, substituent_name), ...] ロカント昇順
    """
    from .substituent import name_substituent
    from .molecule_analyzer import get_atom  # 遅延インポート

    # ring_set: 命名原子 + 縮合環系の橋頭原子（in_ring フラグで拡張）
    # Phase 37: 非縮合結合 (biphenyl の C-C 単結合) は ring_set に含めない
    ring_set = set(ring_chain.ring_atoms)
    _ring_atom_sets = graph.ring_atom_sets or []
    for c_idx in ring_chain.ring_atoms:
        for nb in graph.adjacency[c_idx]:
            if not get_atom(graph, nb).in_ring:
                continue
            # nb が c_idx と同じ環に属する（縮合結合）場合のみ追加
            if any(c_idx in rs and nb in rs for rs in _ring_atom_sets):
                ring_set.add(nb)

    principal_set = set(principal_grp_atoms)
    result: list[tuple[int, str]] = []

    for c_idx in ring_chain.ring_atoms:
        locant = ring_chain.locant_map[c_idx]
        for nb_idx in graph.adjacency[c_idx]:
            nb_atom = get_atom(graph, nb_idx)
            if nb_atom.symbol == "H":
                continue
            if nb_idx in ring_set:
                continue
            if nb_idx in principal_set:
                continue
            sub_name = name_substituent(graph, nb_idx, ring_set | principal_set)
            if sub_name:
                result.append((locant, sub_name))

    result.sort(key=lambda x: (x[0], x[1]))
    return result


# ─── 環の名前組み立て ────────────────────────────────────────────────

def assemble_ring_name(
    ring_chain: RingChain,
    substituents: list[tuple[int, str]],
    principal_grp_type: str,
    suffix_locant: int | None,
    stereo_descriptors: list[str],
) -> str:
    """
    環状化合物の IUPAC 名を組み立てる。

    Args:
        ring_chain: 環の情報
        substituents: [(locant, name)] 置換基リスト
        principal_grp_type: 'alkane', 'alcohol', 'benzene', ...
        suffix_locant: alcohol/ketone の suffix ロカント
        stereo_descriptors: 立体記述子リスト

    Returns:
        IUPAC 名文字列
    """
    from .name_assembler import _build_prefix
    from .constants import FUNCTIONAL_GROUPS

    n = ring_chain.ring_size
    spec = FUNCTIONAL_GROUPS.get(principal_grp_type)

    # ─── ベース名 ────────────────────────────────────────────────
    # base_name オーバーライド (naphthalene 等)
    if ring_chain.base_name is not None:
        # Phase 31: 縮合多環基 (naphthalene 等) に官能基 suffix を付加
        raw_base = ring_chain.base_name
        if spec is not None and principal_grp_type != "alkane":
            suffix = spec.suffix
            if spec.needs_locant and suffix_locant is not None:
                loc = suffix_locant
                # IUPAC 母音省略: suffix が母音始まりなら base の末尾 'e' を省略
                if suffix[0].lower() in "aeiou" and raw_base.endswith("e"):
                    base = f"{raw_base[:-1]}-{loc}-{suffix}"
                else:
                    base = f"{raw_base}-{loc}-{suffix}"
            elif suffix_locant is not None and spec.anchor_c1 and not spec.needs_locant:
                # 環外官能基 (carboxylic acid, aldehyde, amide) on polycyclic
                _exo_sfx = {
                    "carboxylic_acid": "carboxylic acid",
                    "aldehyde": "carbaldehyde",
                    "amide": "carboxamide",
                }
                _sfx = _exo_sfx.get(principal_grp_type)
                base = f"{raw_base}-{suffix_locant}-{_sfx}" if _sfx else raw_base
            else:
                base = raw_base
        else:
            base = raw_base
    elif ring_chain.is_aromatic and n == 6:
        # ベンゼン環の官能基: IUPAC 保留名を使用 (spec.benzene_name)
        if spec is not None and spec.benzene_name is not None:
            base = spec.benzene_name
        else:
            base = "benzene"
    else:
        stem = CHAIN_PREFIX.get(n)
        if stem is None:
            raise ValueError(f"Unsupported ring size: {n}")

        db_locs = ring_chain.double_bond_locants
        if db_locs:
            db_loc_str = ",".join(str(l) for l in sorted(db_locs))
            num_db = len(db_locs)
            # 官能基あり (alcohol, ketone 等) かつ非アルカン系
            if (spec is not None
                    and spec.cyclic_template is not None
                    and spec.suffix not in ("ane", "ene", "yne")):
                loc = suffix_locant if suffix_locant is not None else 1
                suffix_str = spec.suffix
                if num_db == 1:
                    # cyclohex-2-en-1-ol
                    base = f"cyclo{stem}-{db_loc_str}-en-{loc}-{suffix_str}"
                else:
                    # cyclohexa-1,3-dien-1-ol (複数二重結合+官能基)
                    base = f"cyclo{stem}a-{db_loc_str}-dien-{loc}-{suffix_str}"
            else:
                # 純シクロアルケン / シクロアルカジエン
                if num_db == 1:
                    base = f"cyclo{stem}-{db_loc_str}-ene"
                elif num_db == 2:
                    base = f"cyclo{stem}a-{db_loc_str}-diene"
                else:
                    _mult = {3: "triene", 4: "tetraene"}.get(num_db, f"{num_db}ene")
                    base = f"cyclo{stem}a-{db_loc_str}-{_mult}"
        elif spec is not None and spec.cyclic_template is not None:
            loc = suffix_locant if suffix_locant is not None else 1
            base = spec.cyclic_template.format(stem=stem, loc=loc)
        else:
            base = f"cyclo{stem}ane"

    # ─── 接頭辞 ──────────────────────────────────────────────────
    prefix_part = _build_prefix(substituents)

    # ─── Phase 37/106: biphenyl 保留名 ────────────────────────────
    # benzene + phenyl 置換基 → "biphenyl" / 置換ビフェニル
    if base == "benzene" and not ring_chain.double_bond_locants:
        phenyl_subs = [(loc, nm) for loc, nm in substituents if nm == "phenyl"]
        other_subs = [(loc, nm) for loc, nm in substituents if nm != "phenyl"]
        if len(phenyl_subs) == 1:
            if not other_subs:
                return "biphenyl"
            # 置換ビフェニル: phenyl 接続点を C1 として他置換基の locant を振り直す
            phenyl_loc = phenyl_subs[0][0]
            def _bip_renumber(direction: int) -> list[tuple[int, str]]:
                result = []
                for loc, nm in other_subs:
                    delta = (loc - phenyl_loc) % 6 if direction == 1 else (phenyl_loc - loc) % 6
                    result.append((delta + 1, nm))
                return result
            subs_fwd = _bip_renumber(1)
            subs_rev = _bip_renumber(-1)
            locs_fwd = sorted(t[0] for t in subs_fwd)
            locs_rev = sorted(t[0] for t in subs_rev)
            best_subs = subs_fwd if locs_fwd <= locs_rev else subs_rev
            pfx = _build_prefix(best_subs)
            return f"{pfx}biphenyl"

    # ─── 一置換環: ロカント 1 省略 ───────────────────────────────────
    # シクロアルカン・ベンゼン系のみ適用。
    # 環内二重結合があるとき / ナフタレン等縮合環系は常にロカント表示。
    if (substituents
            and ring_chain.base_name is None
            and not ring_chain.double_bond_locants
            and len(set(loc for loc, _ in substituents)) == 1):
        prefix_part = _build_prefix_no_locant_if_single(substituents)

    # ─── 立体記述子 ──────────────────────────────────────────────
    # 複数不斉中心: (1R)-(2S) ではなく (1R,2S) にまとめる
    if stereo_descriptors:
        combined = ",".join(d.strip("()") for d in stereo_descriptors)
        stereo_part = f"({combined})"
    else:
        stereo_part = ""

    # ─── 結合 ────────────────────────────────────────────────────
    if prefix_part:
        result = f"{prefix_part}{base}"
    else:
        result = base

    if stereo_part:
        result = f"{stereo_part}-{result}"

    return result


def _build_prefix_no_locant_if_single(substituents: list[tuple[int, str]]) -> str:
    """
    単一置換ベンゼンの接頭辞: ロカント 1 を省略する。
    例: [(1, 'chloro')] → 'chloro' (ロカントなし)
    複数置換の場合はロカントあり。
    """
    from .name_assembler import _build_prefix

    by_name: dict[str, list[int]] = {}
    for loc, name in substituents:
        by_name.setdefault(name, []).append(loc)

    # 全置換基がロカント 1 だけ → ロカント省略
    all_at_1 = all(all(l == 1 for l in locs) for locs in by_name.values())
    if all_at_1:
        # ロカントなし
        import re as _re
        from .constants import MULTIPLIER

        def _cpd(nm: str) -> bool:
            if _re.search(r"[0-9]", nm):
                return True
            for pfx in ("di", "tri", "tetra", "penta", "fluoro", "chloro", "bromo", "iodo"):
                if nm.startswith(pfx) and len(nm) > len(pfx):
                    return True
            return False

        parts = []
        for name in sorted(by_name.keys()):
            n = len(by_name[name])
            mult = MULTIPLIER.get(n, "")
            entry = f"{mult}{name}"
            parts.append(f"({entry})" if _cpd(entry) else entry)
        return "".join(parts)  # e.g., "chloro", "dimethyl", "(trifluoromethyl)"

    return _build_prefix(substituents)


# ─── 環 vs 側鎖の主鎖選択 ────────────────────────────────────────────

def ring_or_chain_is_principal(
    graph: MoleculeGraph,
    ring: list[int],
    chain_carbons: list[int],
) -> str:
    """
    環と側鎖のどちらを主鎖にするか決定する。

    IUPAC 規則: 炭素数が多い方が主鎖。同じ場合は環を優先。

    Returns:
        'ring' または 'chain'
    """
    ring_size = len(ring)
    chain_size = len(chain_carbons)
    return "ring" if ring_size >= chain_size else "chain"
