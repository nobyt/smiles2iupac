"""
ヘテロ環保留名の検出と命名。

IUPAC 2013 Blue Book P-31.1.3.4 / P-31.1.3.5 / P-31.1.3.6

対応範囲:
  - 芳香族ヘテロ環: pyridine, furan, thiophene, 1H-pyrrole, 1H-imidazole,
                    pyrimidine, pyrazine, pyridazine
  - 飽和ヘテロ環: piperidine, morpholine
  - 置換ヘテロ環: 2-methylpyridine, 3-chloropyridine 等
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .molecule_analyzer import MoleculeGraph


# ─── 優先度付きヘテロ原子記号 ─────────────────────────────────────────
# NH (芳香族 NH) は別扱い: 他の N より優先し "NH" として記録する

_PRIORITY: dict[str, int] = {
    "NH": 0,   # 最高優先 (1H-pyrrole, 1H-imidazole)
    "O":  1,
    "S":  2,
    "N":  3,
    "P":  4,
    "Se": 5,
    "Te": 6,
}


def _atom_sig(graph: "MoleculeGraph", idx: int) -> str:
    """原子インデックスをシグネチャ文字列に変換する。"""
    from .molecule_analyzer import get_atom
    a = get_atom(graph, idx)
    if a.symbol == "N" and a.is_aromatic:
        # build_molecule_graph は AddHs() 後にグラフを構築するため
        # num_hs ではなく隣接 H 原子の有無で NH を判定する
        has_h = any(get_atom(graph, nb).symbol == "H" for nb in graph.adjacency[idx])
        if has_h:
            return "NH"
    if a.symbol == "C":
        return "C"
    return a.symbol


def _canonical_sig(ordered: list[int], graph: "MoleculeGraph") -> tuple[str, ...]:
    """
    環の原子リスト（始端=ヘテロ原子）から正準シグネチャを生成する。

    2つの方向（時計回り・反時計回り）の両方を試し、
    辞書順で小さい方を選択する。
    """
    n = len(ordered)
    fwd = tuple(_atom_sig(graph, ordered[i % n]) for i in range(n))
    rev = tuple(_atom_sig(graph, ordered[(-i) % n]) for i in range(n))
    return min(fwd, rev)


def _hetero_priority(sig: str) -> int:
    return _PRIORITY.get(sig, 99)


def _find_best_start(ring: list[int], graph: "MoleculeGraph") -> list[int]:
    """
    ヘテロ原子の中から最高優先度のものを始端にした環の回転を返す。
    複数ある場合はロカント集合が最小になる回転を選択する。
    同点の場合は外環 C=O ロカント集合が最小の回転を選択する（IUPAC P-14.5）。
    """
    from .molecule_analyzer import get_atom, get_bond_order

    n = len(ring)
    ring_set = set(ring)
    # ヘテロ原子の位置を特定
    hetero_positions = [
        i for i, idx in enumerate(ring)
        if _atom_sig(graph, idx) != "C"
    ]
    if not hetero_positions:
        return ring

    # 最高優先ヘテロ原子の種類を特定
    best_priority = min(_hetero_priority(_atom_sig(graph, ring[i])) for i in hetero_positions)
    best_positions = [
        i for i in hetero_positions
        if _hetero_priority(_atom_sig(graph, ring[i])) == best_priority
    ]

    # 候補の回転を生成: 各最高優先ヘテロ原子を位置1に置く (順/逆の両方向を試す)
    candidates: list[tuple[tuple[str, ...], list[int]]] = []
    for start in best_positions:
        fwd = ring[start:] + ring[:start]
        rev = [ring[(start - i) % n] for i in range(n)]
        for rotated in (fwd, rev):
            sig = _canonical_sig(rotated, graph)
            candidates.append((sig, rotated))

    # ロカント集合が最小の回転を選択
    candidates.sort(key=lambda x: x[0])
    best_sig = candidates[0][0]
    best_candidates = [rot for sig, rot in candidates if sig == best_sig]

    def _hetero_locant_key(rot: list[int]) -> list[int]:
        return sorted(i + 1 for i, a in enumerate(rot) if _atom_sig(graph, a) != "C")

    def _exo_co_locant_key(rot: list[int]) -> list[int]:
        locs = []
        for i, idx in enumerate(rot):
            if get_atom(graph, idx).symbol != "C":
                continue
            has_exo_co = any(
                get_atom(graph, nb).symbol == "O"
                and get_bond_order(graph, idx, nb) == 2.0
                for nb in graph.adjacency[idx]
                if nb not in ring_set
            )
            if has_exo_co:
                locs.append(i + 1)
        return sorted(locs)

    def _composite_key(rot: list[int]) -> tuple:
        return (_hetero_locant_key(rot), _exo_co_locant_key(rot))

    best_candidates.sort(key=_composite_key)
    return best_candidates[0]


# ─── 保留名テーブル ────────────────────────────────────────────────────
# key: (is_aromatic, canonical_signature_tuple)
# value: (retained_name, nh_indicator)
#   nh_indicator: True → "1H-" プレフィクスを付ける
#   注: canonical_sig は始端=最高優先ヘテロ原子、辞書順最小方向

# ─── Hantzsch-Widman 名テーブル ───────────────────────────────────────
# key: (ring_size, heteroatom_sig)  — 単一ヘテロ原子・非芳香族のみ
_HW_NAMES: dict[tuple[int, str], str] = {
    (3, "O"): "oxirane",
    (3, "S"): "thiirane",
    (3, "N"): "aziridine",
    (3, "P"): "phosphirane",
    (4, "O"): "oxetane",
    (4, "S"): "thietane",
    (4, "N"): "azetidine",
    (4, "P"): "phosphetane",
    (5, "O"): "oxolane",
    (5, "S"): "thiolane",
    (5, "N"): "pyrrolidine",
    (5, "P"): "phospholane",
    (6, "O"): "oxane",
    (6, "S"): "thiane",
    (6, "P"): "phosphinane",
    # (6, "N") → piperidine は保留名テーブルで処理済み
    (7, "N"): "azepane",
    (7, "O"): "oxepane",
    (7, "S"): "thiepane",
    (7, "P"): "phosphepane",
    (8, "N"): "azocane",
    (8, "O"): "oxocane",
    (8, "S"): "thiocane",
    (8, "P"): "phosphocane",
}


def _match_hantzsch_widman(ring: list[int], graph: "MoleculeGraph") -> str | None:
    """
    単一ヘテロ原子の飽和環を Hantzsch-Widman 名で返す。
    一致しない場合は None。
    """
    if _is_aromatic_ring(ring, graph):
        return None
    n = len(ring)
    if n not in (3, 4, 5, 6, 7, 8):
        return None
    # ヘテロ原子を収集
    hetero_sigs: list[str] = []
    for idx in ring:
        sig = _atom_sig(graph, idx)
        if sig != "C":
            hetero_sigs.append(sig)
    # 単一ヘテロ原子のみ対応
    if len(hetero_sigs) != 1:
        return None
    return _HW_NAMES.get((n, hetero_sigs[0]))


# ─── Phase 272: 部分不飽和単環ヘテロ環 ─────────────────────────────────
# 非芳香族で環内二重結合を持つ 5/6員環 → "X,Y-dihydro[parent]" 形式

_PARTIAL_UNSAT_PARENT: dict[tuple[int, str], str] = {
    # 5-membered
    (5, "O"): "furan",
    (5, "S"): "thiophene",
    (5, "N"): "pyrrole",
    # 6-membered
    (6, "N"): "pyridine",
    (6, "O"): "pyran",
    (6, "S"): "thiopyran",
}


def _match_partial_unsat(ring: list[int], graph: "MoleculeGraph") -> str | None:
    """
    部分不飽和単一ヘテロ原子環の dihydro 名を返す。
    例: 2,3-dihydrofuran, 1,2-dihydropyridine
    """
    from .molecule_analyzer import get_atom, get_bond_order

    if _is_aromatic_ring(ring, graph):
        return None
    n = len(ring)
    if n not in (5, 6):
        return None

    # 環内二重結合を収集
    ring_set = set(ring)
    seen: set[frozenset] = set()
    db_pairs: list[tuple[int, int]] = []
    for idx in ring:
        for nb in graph.adjacency[idx]:
            if nb in ring_set:
                key = frozenset((idx, nb))
                if key not in seen:
                    seen.add(key)
                    if get_bond_order(graph, idx, nb) == 2.0:
                        db_pairs.append((idx, nb))

    if not db_pairs:
        return None  # 完全飽和 → HW 名で処理

    # 単一ヘテロ原子のみ対応
    hetero_atoms = [idx for idx in ring if _atom_sig(graph, idx) != "C"]
    if len(hetero_atoms) != 1:
        return None
    hetero_idx = hetero_atoms[0]
    hetero_sym = get_atom(graph, hetero_idx).symbol
    has_nh = any(get_atom(graph, nb).symbol == "H" for nb in graph.adjacency[hetero_idx])

    parent = _PARTIAL_UNSAT_PARENT.get((n, hetero_sym))
    if parent is None:
        return None

    # sp3 原子 = 環内二重結合に関与しない炭素 (ヘテロ原子除く)
    db_atoms = {a for pair in db_pairs for a in pair}
    sp3_atoms = [idx for idx in ring if idx not in db_atoms and idx != hetero_idx]

    # 5員環: 1 db → sp3 が 2個 / 6員環: 1 db → sp3 が 3個, 2 db → sp3 が 1個(+N-H)
    # 6員環 + N-H の場合: N がsp3になる → sp3_atoms にも含める可能性
    sp3_with_hetero = sp3_atoms[:]
    if has_nh and hetero_idx not in db_atoms:
        sp3_with_hetero.append(hetero_idx)

    if not sp3_with_hetero:
        return None

    # ヘテロ原子を 1 位として両方向で sp3 ロカントを求め最小を選ぶ
    hi = ring.index(hetero_idx)
    best_locs: list[int] | None = None
    for direction in (1, -1):
        ordered = [ring[(hi + direction * k) % n] for k in range(n)]
        loc_map = {ordered[i]: i + 1 for i in range(n)}
        locs = sorted(loc_map[a] for a in sp3_with_hetero)
        if best_locs is None or locs < best_locs:
            best_locs = locs

    if best_locs is None:
        return None

    # N-H を持つ pyrrole/pyridine 親の indicated hydrogen
    indicated_h = ""
    if hetero_sym == "N" and has_nh and 1 not in best_locs:
        indicated_h = "1H-"

    mult = {2: "di", 3: "tri", 4: "tetra", 5: "penta"}.get(len(best_locs), "")
    locs_str = ",".join(str(l) for l in best_locs)
    return f"{locs_str}-{mult}hydro{indicated_h}{parent}"


# ─── 保留名テーブル ────────────────────────────────────────────────────
# key: (is_aromatic, canonical_signature_tuple)
# value: (retained_name, nh_indicator)
#   nh_indicator: True → "1H-" プレフィクスを付ける
#   注: canonical_sig は始端=最高優先ヘテロ原子、辞書順最小方向

_RETAINED_NAMES: dict[tuple[bool, tuple[str, ...]], tuple[str, bool]] = {
    # 6員芳香族ヘテロ環
    # pyridine: N-C-C-C-C-C (canonical min of fwd/rev)
    (True,  ("N", "C", "C", "C", "C", "C")): ("pyridine",  False),
    # pyrazine: N at 1,4 → min("N","C","C","N","C","C") vs ("N","C","C","N","C","C") same
    (True,  ("N", "C", "C", "N", "C", "C")): ("pyrazine",  False),
    # pyrimidine: N at 1,3 → best_start=[N,C,N,C,C,C], rev=[N,C,C,C,N,C], min=("N","C","C","C","N","C")
    (True,  ("N", "C", "C", "C", "N", "C")): ("pyrimidine", False),
    # pyridazine: N at 1,2 → best_start=[N,N,C,C,C,C], rev=[N,C,C,C,C,N], min=("N","C","C","C","C","N")
    (True,  ("N", "C", "C", "C", "C", "N")): ("pyridazine", False),
    # 5員芳香族ヘテロ環
    (True,  ("O", "C", "C", "C", "C")): ("furan",      False),
    (True,  ("S", "C", "C", "C", "C")): ("thiophene",  False),
    # pyrrole: 1 NH → ("NH","C","C","C","C")
    (True,  ("NH", "C", "C", "C", "C")): ("pyrrole",   True),   # 1H-pyrrole
    # N-substituted pyrrole: N has no H (replaced by alkyl/acyl)
    (True,  ("N",  "C", "C", "C", "C")): ("pyrrole",   False),  # 1-substituted pyrrole
    # imidazole (1H): NH at 1, N at 3 → best_start=[NH,C,N,C,C], rev=[NH,C,C,N,C], min=("NH","C","C","N","C")
    (True,  ("NH", "C", "C", "N", "C")): ("imidazole", True),   # 1H-imidazole
    # N-substituted imidazole
    (True,  ("N",  "C", "C", "N", "C")): ("imidazole", False),  # 1-substituted imidazole
    # 飽和6員ヘテロ環
    (False, ("N", "C", "C", "C", "C", "C")): ("piperidine", False),
    (False, ("O", "C", "C", "N", "C", "C")): ("morpholine",  False),
    (False, ("N", "C", "C", "O", "C", "C")): ("morpholine",  False),
    # ピペラジン: N at 1,4 (飽和)
    (False, ("N", "C", "C", "N", "C", "C")): ("piperazine",  False),
    # NH ピペラジン (N-H あり)
    (False, ("NH", "C", "C", "N",  "C", "C")): ("piperazine",  False),
    (False, ("NH", "C", "C", "NH", "C", "C")): ("piperazine",  False),
    # Phase 154: 5員 O,N 飽和環 (oxazolidine = 1,3-oxazolidine; isoxazolidine = 1,2-oxazolidine)
    (False, ("O", "C", "C", "N", "C")):           ("oxazolidine",   False),
    (False, ("O", "C", "C", "C", "N")):           ("isoxazolidine", False),
    # Phase 154: 5員 S,N 飽和環 (thiazolidine = 1,3-; isothiazolidine = 1,2-)
    (False, ("S", "C", "C", "N", "C")):           ("thiazolidine",   False),
    (False, ("S", "C", "C", "C", "N")):           ("isothiazolidine", False),
    # Phase 154: 6員 O,N 飽和環 (1,3-oxazinane, 1,4-oxazinane 等)
    (False, ("O", "C", "C", "C", "N", "C")):      ("1,3-oxazinane",   False),
    (False, ("O", "C", "N", "C", "C", "C")):      ("1,3-oxazinane",   False),  # alt sig
    # Phase 154: 6員 S,N 飽和環
    (False, ("S", "C", "C", "N", "C", "C")):      ("1,4-thiazinane",  False),
    (False, ("S", "C", "C", "C", "N", "C")):      ("1,3-thiazinane",  False),
    # Phase 152: imidazolidine (飽和5員 N,N 環)
    (False, ("N", "C", "C", "N", "C")):           ("imidazolidine",   False),
    # Phase 153: テトラゾール (5員 aromatic 4-N 環)
    (True,  ("NH", "C", "N", "N", "N")):          ("tetrazole",       True),
    (True,  ("N",  "C", "N", "N", "N")):          ("tetrazole",       False),  # N-substituted
    # Phase 153: 1,3,5-triazinane (6員 tri-N 環; kekulization 失敗で aromatic 扱い)
    (True,  ("NH", "C", "NH", "C", "NH", "C")):   ("1,3,5-triazinane", False),
    (False, ("N",  "C", "N",  "C", "N",  "C")):   ("1,3,5-triazinane", False),
    # Phase 153: hexahydropyrimidine 型 (非芳香族6員 N1,N3 環)
    (False, ("N", "C", "C", "C", "N", "C")):      ("hexahydropyrimidine", False),
    # Phase 269: 5員二ヘテロ芳香族
    (True,  ("NH", "C", "C", "C", "N")): ("pyrazole",    True),   # 1H-pyrazole
    (True,  ("N",  "C", "C", "C", "N")): ("pyrazole",    False),  # 1-substituted
    (True,  ("O",  "C", "C", "C", "N")): ("isoxazole",   False),
    (True,  ("O",  "C", "C", "N", "C")): ("oxazole",     False),
    (True,  ("S",  "C", "C", "N", "C")): ("thiazole",    False),
    (True,  ("S",  "C", "C", "C", "N")): ("isothiazole", False),
    # Phase 269: 6員 triazine/tetrazine
    (True,  ("N",  "C", "C", "C", "N", "N")): ("1,2,3-triazine",    False),
    (True,  ("N",  "C", "C", "N", "C", "N")): ("1,2,4-triazine",    False),
    (True,  ("N",  "C", "N", "C", "N", "C")): ("1,3,5-triazine",    False),
    (True,  ("N",  "C", "N", "N", "C", "N")): ("1,2,4,5-tetrazine", False),
    # Phase 151: 二ヘテロ原子飽和環 (1,3-dioxolane, 1,4-dioxane 等)
    (False, ("O", "C", "C", "O", "C")):           ("1,3-dioxolane",   False),
    (False, ("O", "C", "C", "O", "C", "C")):      ("1,4-dioxane",     False),
    (False, ("O", "C", "C", "C", "O", "C")):      ("1,3-dioxane",     False),
    (False, ("S", "C", "C", "S", "C")):           ("1,3-dithiolane",  False),
    (False, ("S", "C", "C", "C", "S")):           ("1,2-dithiolane",  False),
    (False, ("S", "C", "C", "S", "C", "C")):      ("1,4-dithiane",    False),
    (False, ("S", "C", "C", "C", "C", "S")):      ("1,2-dithiane",    False),
    (False, ("O", "C", "C", "S", "C")):           ("1,3-oxathiolane", False),
    (False, ("O", "C", "C", "S", "C", "C")):      ("1,4-oxathiane",   False),
    # Phase 232: 1,3,5-trioxane (6員 tri-O 環)
    (False, ("O", "C", "O", "C", "O", "C")): ("1,3,5-trioxane", False),
    # Phase 196: 4員環 二ヘテロ原子 (1,3-dioxetane 等)
    (False, ("O", "C", "O", "C")):                ("1,3-dioxetane",   False),
    (False, ("O", "C", "N", "C")):                ("1,3-oxazetidine", False),
    (False, ("N", "C", "O", "C")):                ("1,3-oxazetidine", False),
    (False, ("O", "C", "C", "N")):                ("1,2-oxazetidine", False),
    (False, ("S", "C", "S", "C")):                ("1,3-dithietane",  False),
    (False, ("O", "C", "S", "C")):                ("1,3-oxathietane", False),
}


def _is_aromatic_ring(ring: list[int], graph: "MoleculeGraph") -> bool:
    from .molecule_analyzer import get_atom
    return all(get_atom(graph, idx).is_aromatic for idx in ring)


def _match_retained(ring: list[int], graph: "MoleculeGraph") -> tuple[str, bool, list[int]] | None:
    """
    環が保留名に一致するか確認する。

    Returns:
        (name, is_nh, best_rotation) または None
        - is_nh: True → "1H-" プレフィクス
        - best_rotation: ロカント1が始端のリスト
    """
    is_arom = _is_aromatic_ring(ring, graph)
    best_rotation = _find_best_start(ring, graph)
    sig = _canonical_sig(best_rotation, graph)

    entry = _RETAINED_NAMES.get((is_arom, sig))
    if entry is not None:
        return entry[0], entry[1], best_rotation

    return None


# ─── ロカントマップ構築 ────────────────────────────────────────────────

def _build_locant_map(rotation: list[int]) -> dict[int, int]:
    """始端=1 のロカントマップを返す。"""
    return {idx: i + 1 for i, idx in enumerate(rotation)}


# ─── 置換基収集（ヘテロ環用） ──────────────────────────────────────────

def _collect_hetero_substituents(
    graph: "MoleculeGraph",
    ring_atoms: list[int],
    locant_map: dict[int, int],
    excluded_atoms: "set[int] | None" = None,
) -> list[tuple[int, str]]:
    """
    ヘテロ環の置換基を収集する。
    excluded_atoms: Phase 22 外環官能基などを除外するためのインデックス集合。
    Returns: [(locant, substituent_name), ...]  ソート済み
    """
    from .molecule_analyzer import get_atom, get_bond_order
    from .substituent import name_substituent

    ring_set = set(ring_atoms)
    result: list[tuple[int, str]] = []

    for ring_idx in ring_atoms:
        locant = locant_map[ring_idx]
        ring_atom = get_atom(graph, ring_idx)

        for nb_idx in graph.adjacency[ring_idx]:
            if nb_idx in ring_set:
                continue
            if excluded_atoms and nb_idx in excluded_atoms:
                continue
            nb = get_atom(graph, nb_idx)
            # H は無視
            if nb.symbol == "H":
                continue
            # ヘテロ原子上のHは無視（NH）
            if ring_atom.symbol in ("N", "O", "S") and nb.symbol == "H":
                continue
            sub_name = name_substituent(graph, nb_idx, ring_set)
            result.append((locant, sub_name))

    result.sort()
    return result


# ─── 名前組み立て ─────────────────────────────────────────────────────

def _format_substituents(
    base: str,
    substituents: list[tuple[int, str]],
) -> str:
    """置換基プレフィクスを付与した名前を組み立てる。"""
    if not substituents:
        return base

    import re
    from collections import defaultdict
    from .constants import MULTIPLIER

    # 同一置換基をまとめる: name -> [locants]
    grouped: dict[str, list[int]] = defaultdict(list)
    for locant, name in substituents:
        grouped[name].append(locant)

    def _cpd_parens(nm: str) -> bool:
        """複合置換基名 → True (括弧が必要)"""
        if re.search(r"[0-9]", nm):
            return True
        for pfx in ("di", "tri", "tetra", "penta", "fluoro", "chloro", "bromo", "iodo"):
            if nm.startswith(pfx) and len(nm) > len(pfx):
                return True
        return False

    parts: list[str] = []
    for name in sorted(grouped):
        locs = sorted(grouped[name])
        parens = _cpd_parens(name)
        n_locs = len(locs)
        if n_locs == 1:
            parts.append(f"{locs[0]}-({name})" if parens else f"{locs[0]}-{name}")
        else:
            loc_str = ",".join(str(l) for l in locs)
            mult = MULTIPLIER.get(n_locs, f"{n_locs}")
            parts.append(f"{loc_str}-{mult}({name})" if parens else f"{loc_str}-{mult}{name}")

    # アルファベット順（先頭の数字・ハイフンを除いて比較）
    def alpha_key(s: str) -> str:
        return re.sub(r"^[\d,\-]+", "", s)

    sorted_parts = sorted(parts, key=alpha_key)
    # 置換基プレフィクスを base 名に連結 ("1H-" prefix の場合はハイフンを挟む)
    # 例: "3-chloro" + "2-methyl" + "pyridine" → "3-chloro-2-methylpyridine"
    # 例: "4-methyl" + "1H-imidazole" → "4-methyl-1H-imidazole"
    sep = "-" if base[:1].isdigit() else ""
    if len(sorted_parts) == 1:
        return f"{sorted_parts[0]}{sep}{base}"
    prefix = "-".join(sorted_parts[:-1]) + "-" + sorted_parts[-1]
    return f"{prefix}{sep}{base}"


def _apply_hetero_suffixes(
    full_base: str,
    substituents: list[tuple[int, str]],
) -> str:
    """
    ヘテロ芳香族の principal group を -amine/-ol/-thiol サフィックス形式に変換して返す。
    いずれも含まない場合は _format_substituents の結果を返す。
    """
    if not substituents:
        return full_base

    from .constants import MULTIPLIER

    for sub_nm, suffix, elide_e in (
        ("amino", "amine", True),
        ("hydroxy", "ol", True),
        ("sulfanyl", "thiol", False),
    ):
        entries = [(loc, nm) for loc, nm in substituents if nm == sub_nm]
        if not entries:
            continue
        other = [(loc, nm) for loc, nm in substituents if nm != sub_nm]
        locs = sorted(loc for loc, _ in entries)
        loc_str = ",".join(str(l) for l in locs)
        mult = MULTIPLIER.get(len(locs), "") if len(locs) > 1 else ""
        stem = full_base[:-1] if (elide_e and full_base.endswith("e")) else full_base
        base_with_suffix = f"{stem}-{loc_str}-{mult}{suffix}"
        if not other:
            return base_with_suffix
        return _format_substituents(base_with_suffix, other)

    return _format_substituents(full_base, substituents)


# ─── 縮合ヘテロ芳香族保留名テーブル (Phase 17) ───────────────────────────
# key: RDKit canonical SMILES (AddHs 前)
# value: IUPAC 保留名
_FUSED_HETERO_RETAINED: dict[str, str] = {
    "c1ccc2ncccc2c1":   "quinoline",
    "c1ccc2cnccc2c1":   "isoquinoline",
    "c1ccc2[nH]ccc2c1": "1H-indole",
    "c1ccc2[nH]cnc2c1": "1H-benzimidazole",
    "c1ccc2nccnc2c1":   "quinoxaline",        # benzo[g]pyrazine (Phase 130)
    "c1ccc2occc2c1":    "benzofuran",
    "c1ccc2sccc2c1":    "benzo[b]thiophene",
    "c1ncc2[nH]cnc2n1": "9H-purine",
    "c1ccc2nc3ccccc3cc2c1": "acridine",
    # Phase 130: 追加縮合ヘテロ芳香族
    "c1ccc2ocnc2c1":    "1,3-benzoxazole",
    "c1ccc2scnc2c1":    "1,3-benzothiazole",
    "c1ccc2ncncc2c1":   "quinazoline",
    "c1ccc2nnccc2c1":   "cinnoline",
    "c1ccc2cnncc2c1":   "phthalazine",
    "c1ccc2nc3ccccc3nc2c1": "phenazine",
    "c1cnc2nccnc2n1":   "pteridine",
    # Phase 131: 三環式縮合ヘテロ芳香族
    "c1ccc2c(c1)[nH]c1ccccc12": "9H-carbazole",
    "c1ccc2c(c1)sc1ccccc12":    "dibenzothiophene",
    "c1ccc2c(c1)oc1ccccc12":    "dibenzofuran",
    # Phase 132: ナフチリジン全6異性体 (canonical SMILES = RDKit canonical)
    "c1cnc2cccnc2c1":   "1,5-naphthyridine",
    "c1cnc2ccncc2c1":   "1,6-naphthyridine",
    "c1cnc2cnccc2c1":   "1,7-naphthyridine",
    "c1cnc2ncccc2c1":   "1,8-naphthyridine",
    "c1cc2cnccc2cn1":   "2,6-naphthyridine",
    "c1cc2ccncc2cn1":   "2,7-naphthyridine",
    # Phase 132: 縮合ヘテロ芳香族追加 (IUPAC 2013 P-31.1.3)
    "c1ccn2cccc2c1":    "indolizine",
    "c1ccc2[nH]ncc2c1": "1H-indazole",
    "c1ccc2[nH]nnc2c1": "1H-benzotriazole",
    "c1ccc2nocc2c1":    "1,2-benzisoxazole",
    "c1ccc2[nH]occ2c1": "1,2-benzisoxazole",
    "c1ccc2oncc2c1":    "2,1,3-benzoxadiazole",
    "c1ccc2onnc2c1":    "2,1,3-benzoxadiazole",
    # Phase 133: 部分飽和縮合環 保留名 (IUPAC 2013 P-31.1.2, P-31.1.6)
    "c1ccc2c(c1)CCC2":  "indane",
    "c1ccc2c(c1)CCCC2": "1,2,3,4-tetrahydronaphthalene",
    "c1ccc2c(c1)CCN2":  "indoline",
    "c1ccc2c(c1)CCO2":  "1,3-dihydro-2-benzofuran",
    "c1ccc2c(c1)CCS2":  "1,3-dihydro-2-benzothiophene",
    "c1ccc2c(c1)CCCO2": "chromane",
    "c1ccc2c(c1)CCCN2": "1,2,3,4-tetrahydroquinoline",
    "c1ccc2c(c1)CCCS2": "thiochroman",
    "O=C1CCc2ccccc21":  "indan-1-one",
    "O=C1CC(=O)c2ccccc21": "indane-1,3-dione",
    "O=C1Cc2ccccc2N1":  "indolin-2-one",
    # Phase 264: 縮合無水物・マレイン酸無水物
    "O=C1OC(=O)c2ccccc21": "isobenzofuran-1,3-dione",
    "O=C1C=CC(=O)O1":      "furan-2,5-dione",
    # Phase 267: 1,3-ベンゾジオキソール (メチレンジオキシベンゼン骨格)
    "c1ccc2c(c1)OCO2":  "1,3-benzodioxole",
    # Phase 134: 追加縮合環保留名 (IUPAC 2013 P-31.1.3, fluorene/xanthene 系)
    "c1ccc2c(c1)Cc1ccccc1-2":  "fluorene",
    "c1ccc2c(c1)Cc1ccccc1O2":  "xanthene",
    "c1ccc2c(c1)Cc1ccccc1S2":  "thioxanthene",
    "O=c1ccc2ccccc2o1":        "coumarin",
    "O=c1occc2ccccc12":        "isocoumarin",
    "C1=COc2ccccc2C1":         "2H-chromene",
    "C1=Cc2ccccc2OC1":         "4H-chromene",
    "c1ccc2c(c1)Nc1ccccc1O2":  "phenoxazine",
    "c1ccc2c(c1)Nc1ccccc1S2":  "phenothiazine",
    "O=c1c2ccccc2oc2ccccc12":  "xanthen-9-one",
    "C1=Nc2cccc3cccc1c23":     "perimidine",
    # Phase 142: 追加ヘテロ芳香族 (セレノフェン、縮合二環式)
    "c1cc[se]c1":       "selenophene",
    # Phase 255: テルロフェン
    "c1cc[te]c1":       "tellurophene",
    "c1ccn2ccnc2c1":    "imidazo[1,2-a]pyridine",
    "c1cnc2ccnn2c1":    "imidazo[1,2-b]pyridazine",
    "c1cnc2ccsc2c1":    "thieno[2,3-b]pyridine",
    "c1cnc2cn[nH]c2c1": "pyrazolo[1,5-a]pyrimidine",
    "c1cnc2[nH]ccc2c1": "1H-pyrrolo[2,3-b]pyridine",
    "c1cnc2cc[nH]c2c1": "1H-pyrrolo[3,2-b]pyridine",
    "c1cc2ccsc2s1":     "thieno[3,2-b]thiophene",
    "c1cnc2nc[nH]c2c1": "3H-imidazo[4,5-b]pyridine",
    "c1ncc2nc[nH]c2n1": "6H-purine",
    # Phase 141: 単環式ヘテロ芳香族 (Hantzsch-Widman が未対応) IUPAC 2013 P-31.1.3
    # 5-membered, two different heteroatoms
    "c1cnoc1":   "isoxazole",
    "c1cocn1":   "oxazole",
    "c1cscn1":   "thiazole",
    "c1cnsc1":   "isothiazole",
    # 5-membered oxadiazoles
    "c1conn1":   "1,2,3-oxadiazole",
    "c1ncon1":   "1,2,4-oxadiazole",
    "c1cnon1":   "1,2,5-oxadiazole",
    "c1nnco1":   "1,3,4-oxadiazole",
    # 5-membered thiadiazoles
    "c1csnn1":   "1,2,3-thiadiazole",
    "c1ncsn1":   "1,2,4-thiadiazole",
    "c1cnsn1":   "1,2,5-thiadiazole",
    "c1nncs1":   "1,3,4-thiadiazole",
    # 5-membered with adjacent N (pyrazole, triazoles, tetrazole)
    "c1cn[nH]c1":  "1H-pyrazole",
    "c1cn[nH]n1":  "1H-1,2,3-triazole",
    "c1c[nH]nn1":  "2H-1,2,3-triazole",
    "c1nc[nH]n1":  "1H-1,2,4-triazole",
    "c1nn[nH]n1":  "1H-tetrazole",
    "c1nnn[nH]1":  "1H-tetrazole",  # alternate canonical form
    # 6-membered triazines and tetrazine
    "c1cnnnc1":    "1,2,3-triazine",
    "c1cnncn1":    "1,2,4-triazine",
    "c1ncncn1":    "1,3,5-triazine",
    "c1nncnn1":    "1,2,4,5-tetrazine",
    # Phase 140: 二ヘテロ原子飽和環 (Hantzsch-Widman, IUPAC 2013 P-31.1.3.4)
    "C1COCCO1":  "1,4-dioxane",
    "C1COCO1":   "1,3-dioxolane",
    "C1COCOC1":  "1,3-dioxane",
    "C1CSCS1":   "1,3-dithiolane",
    "C1CSSC1":   "1,2-dithiolane",
    "C1CSCCS1":  "1,4-dithiane",
    "C1CCSSC1":  "1,2-dithiane",
    "C1CSCCO1":  "1,4-oxathiane",
    "C1CSCO1":   "1,3-oxathiolane",
    # Phase 260: 1,2-ジオキサン・1,2-ジオキソラン (隣接 O-O を持つ環)
    "C1CCOOC1":  "1,2-dioxane",
    "C1COOC1":   "1,2-dioxolane",
    # Phase 138: 多環芳香族炭化水素 保留名 (IUPAC 2013 P-31.1.2)
    "c1ccc2cccc-2cc1":                            "azulene",
    "C1=Cc2ccccc2C1":                             "1H-indene",
    "C1=Cc2cccc3cccc1c23":                        "acenaphthylene",
    "c1cc2c3c(cccc3c1)CC2":                       "acenaphthene",
    "c1cc2cccc3ccc4ccccc4c3c-2c1":                "fluoranthene",
    "c1ccc2cc3c(ccc4ccccc43)cc2c1":               "chrysene",
    "c1ccc2c(c1)c1ccccc1c1ccccc21":               "triphenylene",
    "c1cc2ccc3cccc4ccc(c1)c2c34":                 "pyrene",
    "C1=Cc2cccc3ccc4ccc5cccc1c5c4c23":            "perylene",
    "c1ccc2c3c4c(cccc-3cc2c1)ccc1ccccc14":        "benzo[a]pyrene",
    "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61":     "coronene",
}

# atom index in canonical base SMILES → IUPAC locant (None = ring junction, skip)
_FUSED_LOCANT_MAP: dict[str, dict[int, int | None]] = {
    "c1ccc2ncccc2c1":   {0: 6, 1: 7, 2: 8, 3: None, 4: 1, 5: 2, 6: 3, 7: 4, 8: None, 9: 5},
    "c1ccc2cnccc2c1":   {0: 6, 1: 7, 2: 8, 3: None, 4: 1, 5: 2, 6: 3, 7: 4, 8: None, 9: 5},
    "c1ccc2[nH]ccc2c1": {0: 5, 1: 6, 2: 7, 3: None, 4: 1, 5: 2, 6: 3, 7: None, 8: 4},
    "c1ccc2[nH]cnc2c1": {0: 5, 1: 6, 2: 7, 3: None, 4: 1, 5: 2, 6: 3, 7: None, 8: 4},
    "c1ccc2occc2c1":    {0: 5, 1: 6, 2: 7, 3: None, 4: 1, 5: 2, 6: 3, 7: None, 8: 4},
    "c1ccc2sccc2c1":    {0: 5, 1: 6, 2: 7, 3: None, 4: 1, 5: 2, 6: 3, 7: None, 8: 4},
    # Phase 158: 9H-carbazole (N at atom 6 = IUPAC position 9)
    "c1ccc2c(c1)[nH]c1ccccc12": {
        0: 2, 1: 3, 2: 4, 3: None, 4: None, 5: 1,
        6: 9, 7: None, 8: 8, 9: 7, 10: 6, 11: 5, 12: None,
    },
}


def _try_fused_hetero_retained(graph: "MoleculeGraph") -> str | None:
    """
    縮合ヘテロ芳香族保留名テーブルに一致する場合、名前を返す。
    置換体にも対応: 環コア SMILES でマッチし、置換基を命名する。
    """
    if graph.rdkit_mol is None:
        return None
    from rdkit.Chem import MolToSmiles, MolFromSmiles, MolFragmentToSmiles
    canon = MolToSmiles(graph.rdkit_mol)
    # 完全一致チェック (置換なし)
    name = _FUSED_HETERO_RETAINED.get(canon)
    if name is not None:
        return name

    # 置換体チェック: 環原子のみの core SMILES を取得して照合
    all_ring_atoms = sorted({a for ring in graph.ring_atom_sets for a in ring})
    if not all_ring_atoms:
        return None
    core_smi_raw = MolFragmentToSmiles(graph.rdkit_mol, all_ring_atoms, canonical=True)
    # MolFragmentToSmiles は必ずしも MolToSmiles と同じ正準形にならないので正規化する
    _tmp_mol = MolFromSmiles(core_smi_raw)
    if _tmp_mol is not None:
        core_smi = MolToSmiles(_tmp_mol)
    else:
        # Phase 158: N-置換縮合環 (1-methylindole 等): 各 'n' を順に '[nH]' に置換して再試行
        import re as _re
        _bare_n_positions = [m.start() for m in _re.finditer(r"(?<!\[)n(?!\])", core_smi_raw)]
        _found_core = None
        for _pos in _bare_n_positions:
            _alt_raw = (core_smi_raw[:_pos] + "[nH]"
                        + core_smi_raw[_pos + 1:])
            _tmp_alt = MolFromSmiles(_alt_raw)
            if _tmp_alt is not None:
                _alt_canon = MolToSmiles(_tmp_alt)
                if _alt_canon in _FUSED_HETERO_RETAINED:
                    _found_core = _alt_canon
                    break
        core_smi = _found_core if _found_core is not None else core_smi_raw

    base_name = _FUSED_HETERO_RETAINED.get(core_smi)
    locant_map_def = _FUSED_LOCANT_MAP.get(core_smi)
    if base_name is None or locant_map_def is None:
        return None

    # substructure match: base_mol_atom_idx → graph_rdkit_atom_idx
    base_mol = MolFromSmiles(core_smi)
    if base_mol is None:
        return None
    match = graph.rdkit_mol.GetSubstructMatch(base_mol)
    if not match:
        return None
    # match[base_idx] = rdkit_atom_idx in graph.rdkit_mol

    # rdkit atom idx → locant
    rdkit_to_locant: dict[int, int] = {}
    for base_idx, rdkit_idx in enumerate(match):
        loc = locant_map_def.get(base_idx)
        if loc is not None:
            rdkit_to_locant[rdkit_idx] = loc

    # graph atom idx → rdkit atom idx (by atom order; same since no H reordering)
    # Non-ring neighbors of ring atoms are substituents
    from .molecule_analyzer import get_atom
    from .substituent import name_substituent

    ring_set = set(all_ring_atoms)
    substituents: list[tuple[int, str]] = []

    # Phase 148: 外環主官能基（カルボン酸・アルデヒド・アミド・ニトリル）の検出
    _EXOCYCLIC_SUFFIX: dict[str, str] = {
        "carboxylic_acid": "carboxylic acid",
        "aldehyde": "carbaldehyde",
        "amide": "carboxamide",
        "nitrile": "carbonitrile",
    }
    from .functional_group import detect_groups, principal_group as _pgrp_fn
    fg_groups = detect_groups(graph)
    pgrp_ex = _pgrp_fn(fg_groups)
    exo_loc: int | None = None
    exo_suffix: str | None = None
    exo_anchor_c: int | None = None

    if pgrp_ex is not None and pgrp_ex.group_type in _EXOCYCLIC_SUFFIX:
        anchor_c = pgrp_ex.atom_indices[0]
        for nb_idx in graph.adjacency[anchor_c]:
            if nb_idx in ring_set:
                loc_candidate = rdkit_to_locant.get(nb_idx)
                if loc_candidate is not None:
                    exo_loc = loc_candidate
                    exo_suffix = _EXOCYCLIC_SUFFIX[pgrp_ex.group_type]
                    exo_anchor_c = anchor_c
                    break

    for base_idx, rdkit_idx in enumerate(match):
        loc = locant_map_def.get(base_idx)
        if loc is None:
            continue  # junction atom
        # Find non-ring neighbors of this ring atom in graph
        for nb_idx in graph.adjacency[rdkit_idx]:
            if nb_idx in ring_set:
                continue
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "H":
                continue
            if nb_idx == exo_anchor_c:
                continue  # handled as principal group suffix
            sub_name = name_substituent(graph, nb_idx, ring_set)
            substituents.append((loc, sub_name))

    if exo_suffix is not None and exo_loc is not None:
        result = f"{base_name}-{exo_loc}-{exo_suffix}"
        if substituents:
            result = _apply_hetero_suffixes(result, substituents)
        return result

    if not substituents:
        return base_name
    return _apply_hetero_suffixes(base_name, substituents)


def _is_fused_hetero_aromatic(graph: "MoleculeGraph") -> bool:
    """
    複数の芳香族環が縮合（2原子以上共有）しており、かつヘテロ原子を含む系かどうかを判定する。
    ビフェニル型（単結合で連結、非縮合）は False を返す。
    テーブル未収録の縮合ヘテロ芳香族に対して誤命名を防ぐ用途。
    """
    aromatic_rings = [r for r in graph.ring_atom_sets if _is_aromatic_ring(list(r), graph)]
    if len(aromatic_rings) < 2:
        return False
    # 縮合チェック: 任意の2環が 2原子以上を共有しているか
    has_fusion = False
    for i in range(len(aromatic_rings)):
        for j in range(i + 1, len(aromatic_rings)):
            if len(set(aromatic_rings[i]) & set(aromatic_rings[j])) >= 2:
                has_fusion = True
                break
        if has_fusion:
            break
    if not has_fusion:
        return False
    # 少なくとも1つのヘテロ原子が芳香族環に含まれるか
    from .molecule_analyzer import get_atom
    ring_atoms = {idx for ring in aromatic_rings for idx in ring}
    return any(
        get_atom(graph, idx).symbol not in ("C", "H") and get_atom(graph, idx).is_aromatic
        for idx in ring_atoms
    )


def name_heterocycle(graph: "MoleculeGraph") -> str | None:
    """
    分子がヘテロ環保留名に一致する場合、その IUPAC 名を返す。
    一致しない場合は None を返す。
    """
    from .ring_handler import find_hetero_rings

    # Phase 17: 縮合ヘテロ芳香族保留名を先行チェック
    fused_name = _try_fused_hetero_retained(graph)
    if fused_name is not None:
        return fused_name

    # 未収録の縮合ヘテロ芳香族は誤命名を防ぐため None を返す
    if _is_fused_hetero_aromatic(graph):
        return None

    hetero_rings = find_hetero_rings(graph)
    if not hetero_rings:
        return None

    # 単環のみ対応
    if len(hetero_rings) != 1:
        return None

    ring = hetero_rings[0]

    # 1. Phase 272: 部分不飽和環を最初にチェック (環内二重結合あり)
    # 飽和環や芳香族環は db_pairs が空 or is_aromatic で None を返すので安全
    pu_name = _match_partial_unsat(ring, graph)
    if pu_name is not None:
        rotation = _find_best_start(ring, graph)
        base_name, is_nh = pu_name, False
    else:
        # 2. 保留名テーブルを試みる
        match = _match_retained(ring, graph)
        if match is not None:
            base_name, is_nh, rotation = match
        else:
            # 3. Hantzsch-Widman 名を試みる
            hw_name = _match_hantzsch_widman(ring, graph)
            if hw_name is None:
                return None
            rotation = _find_best_start(ring, graph)
            base_name, is_nh = hw_name, False

    base_name, is_nh, rotation = (base_name, is_nh, rotation)
    locant_map = _build_locant_map(rotation)

    # 完全ベース名
    nh_prefix = "1H-" if is_nh else ""
    full_base = f"{nh_prefix}{base_name}"

    # Phase 22: 外環主官能基（カルボン酸・アルデヒド・アミド・ニトリル）の検出
    _EXOCYCLIC_SUFFIX: dict[str, str] = {
        "carboxylic_acid": "carboxylic acid",
        "aldehyde": "carbaldehyde",
        "amide": "carboxamide",
        "nitrile": "carbonitrile",
    }
    from .functional_group import detect_groups, principal_group as _principal_group
    fg_groups = detect_groups(graph)
    pgrp_ex = _principal_group(fg_groups)
    ring_set_ex = set(ring)

    if pgrp_ex is not None and pgrp_ex.group_type in _EXOCYCLIC_SUFFIX:
        anchor_c = pgrp_ex.atom_indices[0]
        # 主官能基アンカーがこのヘテロ芳香環に隣接していない場合 → _name_cyclic に委譲
        if not any(nb in ring_set_ex for nb in graph.adjacency[anchor_c]):
            return None
        for nb_idx in graph.adjacency[anchor_c]:
            if nb_idx in ring_set_ex:
                L = locant_map[nb_idx]
                n = len(rotation)
                L_rev = n + 2 - L
                if L_rev < L:
                    rev_rotation = [rotation[0]] + list(reversed(rotation[1:]))
                    locant_map = _build_locant_map(rev_rotation)
                    L = L_rev
                fg_suffix = _EXOCYCLIC_SUFFIX[pgrp_ex.group_type]
                other_subs = _collect_hetero_substituents(
                    graph, ring, locant_map, excluded_atoms={anchor_c}
                )
                result = f"{full_base}-{L}-{fg_suffix}"
                if other_subs:
                    result = _format_substituents(result, other_subs)
                return result

    # Phase 25/30: ラクタム (N 含有) / ラクトン (O 含有) 検出
    # 環内 C=O (exocyclic =O) を持つ環 → "{base}-{loc}-one"
    ring_set_lactam = set(ring)
    from .molecule_analyzer import get_bond_order as _get_bo_lact, get_atom as _get_atom_lact

    def _has_exo_dbl_o(c_idx: int) -> int | None:
        """環外に =O を持つ C のその O インデックスを返す。なければ None。"""
        for nb in graph.adjacency[c_idx]:
            if nb in ring_set_lactam:
                continue
            if (_get_atom_lact(graph, nb).symbol == "O"
                    and _get_bo_lact(graph, c_idx, nb) == 2.0):
                return nb
        return None

    # Phase 152/153: 複数 exo C=O グループを持つ環 → "-dione"/"-trione"
    # Phase 58/127 より前に実行して 3+ C=O (trione) も正しく処理する
    exo_co_ring_cs_e = []
    exo_co_oxygens_e = []
    for ring_idx_e in ring:
        if _get_atom_lact(graph, ring_idx_e).symbol != "C":
            continue
        exo_o_e = _has_exo_dbl_o(ring_idx_e)
        if exo_o_e is not None:
            exo_co_ring_cs_e.append(ring_idx_e)
            exo_co_oxygens_e.append(exo_o_e)
    if len(exo_co_ring_cs_e) >= 2:
        locs_e = [locant_map[c] for c in exo_co_ring_cs_e]
        rev_rotation_e = [rotation[0]] + list(reversed(rotation[1:]))
        locant_map_rev_e = _build_locant_map(rev_rotation_e)
        locs_rev_e = [locant_map_rev_e[c] for c in exo_co_ring_cs_e]
        if sorted(locs_rev_e) < sorted(locs_e):
            locant_map = locant_map_rev_e
            locs_e = locs_rev_e
        locs_sorted_e = sorted(locs_e)
        loc_str_e = ",".join(str(l) for l in locs_sorted_e)
        suffix_e = "dione" if len(locs_sorted_e) == 2 else "trione"
        _eb = full_base  # -dione/-trione starts with consonant: no elision
        dione_name_e = f"{_eb}-{loc_str_e}-{suffix_e}"
        excl_e = set(exo_co_oxygens_e)
        other_subs_e = _collect_hetero_substituents(
            graph, ring, locant_map, excluded_atoms=excl_e
        )
        if not other_subs_e:
            return dione_name_e
        return _format_substituents(dione_name_e, other_subs_e)

    # Phase 58: 環状イミド検出 (N が両側の C(=O) に挟まれている)
    # O=C1NC(=O)CC1 → pyrrolidine-2,5-dione (exo C=O 1 個の場合のフォールバック)
    for n_idx_im in ring:
        if _get_atom_lact(graph, n_idx_im).symbol != "N":
            continue
        n_ring_c = [nb for nb in graph.adjacency[n_idx_im]
                    if nb in ring_set_lactam and _get_atom_lact(graph, nb).symbol == "C"]
        if len(n_ring_c) == 2:
            o1 = _has_exo_dbl_o(n_ring_c[0])
            o2 = _has_exo_dbl_o(n_ring_c[1])
            if o1 is not None and o2 is not None:
                loc1 = locant_map[n_ring_c[0]]
                loc2 = locant_map[n_ring_c[1]]
                # ロカントが小さい方から順に並べる
                locs = sorted([loc1, loc2])
                # より小さいロカント集合になる向きを選ぶ (ロカント最小化)
                n_ring = len(ring)
                locs_rev = sorted([n_ring + 2 - loc1, n_ring + 2 - loc2])
                if locs_rev < locs:
                    rev_rotation = [rotation[0]] + list(reversed(rotation[1:]))
                    locant_map = _build_locant_map(rev_rotation)
                    locs = locs_rev
                loc_str = ",".join(str(l) for l in locs)
                _imide_base = full_base  # -dione starts with consonant: no elision
                imide_name = f"{_imide_base}-{loc_str}-dione"
                excl_o = {o1, o2}
                other_subs = _collect_hetero_substituents(
                    graph, ring, locant_map, excluded_atoms=excl_o
                )
                if not other_subs:
                    return imide_name
                return _format_substituents(imide_name, other_subs)

    # Phase 127: 環状無水物 (oxa-dione): 環内 O の両隣 C が共に exo C=O → "{base}-{loc1},{loc2}-dione"
    for o_idx_ca in ring:
        if _get_atom_lact(graph, o_idx_ca).symbol != "O":
            continue
        o_ring_c = [nb for nb in graph.adjacency[o_idx_ca]
                    if nb in ring_set_lactam and _get_atom_lact(graph, nb).symbol == "C"]
        if len(o_ring_c) == 2:
            exo1 = _has_exo_dbl_o(o_ring_c[0])
            exo2 = _has_exo_dbl_o(o_ring_c[1])
            if exo1 is not None and exo2 is not None:
                loc1 = locant_map[o_ring_c[0]]
                loc2 = locant_map[o_ring_c[1]]
                locs = sorted([loc1, loc2])
                loc_str = ",".join(str(l) for l in locs)
                _anhy_base = full_base  # -dione starts with consonant: no elision
                dione_name = f"{_anhy_base}-{loc_str}-dione"
                excl_o = {exo1, exo2}
                other_subs = _collect_hetero_substituents(
                    graph, ring, locant_map, excluded_atoms=excl_o
                )
                if not other_subs:
                    return dione_name
                return _format_substituents(dione_name, other_subs)

    # Phase 152: 複数 exo C=O グループを持つ環 → "-dione"/"-trione" (例: piperazine-2,5-dione)
    exo_co_ring_cs = []
    exo_co_oxygens = []
    for ring_idx_m in ring:
        if _get_atom_lact(graph, ring_idx_m).symbol != "C":
            continue
        exo_o = _has_exo_dbl_o(ring_idx_m)
        if exo_o is not None:
            exo_co_ring_cs.append(ring_idx_m)
            exo_co_oxygens.append(exo_o)
    if len(exo_co_ring_cs) >= 2:
        # 複数 C=O → dione/trione suffix
        locs_m = [locant_map[c] for c in exo_co_ring_cs]
        n_ring_m = len(ring)
        # ロカント集合を最小化: すべての C=O ロカントの和が最小の向きを選ぶ
        rev_rotation_m = [rotation[0]] + list(reversed(rotation[1:]))
        locant_map_rev_m = _build_locant_map(rev_rotation_m)
        locs_rev_m = [locant_map_rev_m[c] for c in exo_co_ring_cs]
        if sorted(locs_rev_m) < sorted(locs_m):
            locant_map = locant_map_rev_m
            locs_m = locs_rev_m
        locs_sorted = sorted(locs_m)
        loc_str_m = ",".join(str(l) for l in locs_sorted)
        suffix_m = "dione" if len(locs_sorted) == 2 else "trione"
        _dm = full_base  # -dione/-trione starts with consonant: no elision
        dione_name_m = f"{_dm}-{loc_str_m}-{suffix_m}"
        excl_m = set(exo_co_oxygens)
        other_subs_m = _collect_hetero_substituents(
            graph, ring, locant_map, excluded_atoms=excl_m
        )
        if not other_subs_m:
            return dione_name_m
        return _format_substituents(dione_name_m, other_subs_m)

    has_ring_hetero = any(_get_atom_lact(graph, idx).symbol not in ("C", "H") for idx in ring)
    if has_ring_hetero:
        for ring_idx in ring:
            ring_atom = _get_atom_lact(graph, ring_idx)
            if ring_atom.symbol != "C":
                continue
            for nb_idx in graph.adjacency[ring_idx]:
                if nb_idx in ring_set_lactam:
                    continue
                nb = _get_atom_lact(graph, nb_idx)
                if nb.symbol == "O" and _get_bo_lact(graph, ring_idx, nb_idx) == 2.0:
                    loc = locant_map[ring_idx]
                    n_ring = len(ring)
                    loc_rev = n_ring + 2 - loc
                    if loc_rev < loc:
                        rev_rotation = [rotation[0]] + list(reversed(rotation[1:]))
                        locant_map = _build_locant_map(rev_rotation)
                        loc = locant_map[ring_idx]
                    lactam_base = (full_base[:-1] if full_base.endswith("e") else full_base)
                    lactam_name = f"{lactam_base}-{loc}-one"
                    other_subs = _collect_hetero_substituents(
                        graph, ring, locant_map, excluded_atoms={nb_idx}
                    )
                    if not other_subs:
                        return lactam_name
                    return _format_substituents(lactam_name, other_subs)

    # 置換基収集: 全有効回転候補を試して最小ロカント集合 → 同点時アルファベット順
    n_ring = len(rotation)
    fwd_hetero = sorted(
        i + 1 for i, a in enumerate(rotation) if _atom_sig(graph, a) != "C"
    )

    # 各ヘテロ原子を起点として逆方向回転を生成し、ヘテロ原子ロカント集合が同じものを収集
    rot_candidates = [rotation]
    for i in range(n_ring):
        if _atom_sig(graph, rotation[i]) == "C":
            continue
        alt = [rotation[(i - j) % n_ring] for j in range(n_ring)]
        if alt == rotation:
            continue
        alt_hetero = sorted(
            k + 1 for k, a in enumerate(alt) if _atom_sig(graph, a) != "C"
        )
        # Only swap to an alternative that starts with the same element (e.g. N↔N ok,
        # N↔O not allowed: in isoxazole O must remain at locant 1, not N).
        if alt_hetero == fwd_hetero and (
            _get_atom_lact(graph, alt[0]).symbol == _get_atom_lact(graph, rotation[0]).symbol
        ):
            rot_candidates.append(alt)

    # 各候補の置換基ロカント集合を収集
    candidate_subs = [
        _collect_hetero_substituents(graph, ring, _build_locant_map(rot))
        for rot in rot_candidates
    ]

    # 最小ロカント集合を選択し、同点時はアルファベット順で tie-break
    def _subs_key(subs: list[tuple[int, str]]) -> tuple[list[int], list[str]]:
        ordered = sorted(subs, key=lambda x: x[0])
        return (sorted(loc for loc, _ in subs), [nm for _, nm in ordered])

    candidate_subs.sort(key=_subs_key)
    substituents = candidate_subs[0]

    if not substituents:
        return full_base

    return _apply_hetero_suffixes(full_base, substituents)
