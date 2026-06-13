"""
置換基の命名。

主鎖から分岐する置換基を IUPAC 規則に従って命名する。
炭素分岐は再帰的に命名（複合置換基対応）。
ハロゲンは HALOGEN_PREFIX から直接取得。
"""

from __future__ import annotations

from .constants import CHAIN_PREFIX, HALOGEN_PREFIX, MULTIPLIER
from .molecule_analyzer import MoleculeGraph, get_atom


def collect_substituents(
    graph: MoleculeGraph,
    chain_atom_indices: list[int],
    locant_map: dict[int, int],
    principal_group_atom_indices: list[int],
) -> list[tuple[int, str]]:
    """
    主鎖上の各炭素に付く置換基を列挙する。

    Args:
        graph: MoleculeGraph
        chain_atom_indices: 主鎖の炭素インデックスリスト (ロカント順)
        locant_map: atom_idx -> locant
        principal_group_atom_indices: 主官能基の原子インデックス (除外対象)

    Returns:
        [(locant, substituent_name), ...] のリスト（ロカント昇順）
    """
    from .molecule_analyzer import get_bond_order as _gbo_cs

    chain_set = set(chain_atom_indices)
    principal_set = set(principal_group_atom_indices)
    excluded = chain_set | principal_set
    result: list[tuple[int, str]] = []

    for c_idx in chain_atom_indices:
        locant = locant_map[c_idx]
        for nb_idx in graph.adjacency[c_idx]:
            nb = get_atom(graph, nb_idx)
            # H・主鎖炭素・主官能基の原子はスキップ
            if nb.symbol == "H":
                continue
            if nb_idx in chain_set:
                continue
            if nb_idx in principal_set:
                continue
            # Phase 280: 主鎖炭素への exo C=C 二重結合 → alkylidene 接頭辞
            if nb.symbol == "C" and _gbo_cs(graph, c_idx, nb_idx) == 2.0:
                sub_name = _name_chain_alkylidene(graph, nb_idx, excluded)
            else:
                sub_name = name_substituent(graph, nb_idx, excluded)
            if sub_name:
                result.append((locant, sub_name))

    result.sort(key=lambda x: (x[0], x[1]))
    return result


def _name_chain_alkylidene(
    graph: MoleculeGraph,
    root_idx: int,
    excluded: set[int],
) -> str:
    """exo C=C に対するアルキリデン接頭辞 (=CH2 → 'methylidene' 等)。"""
    carbons = _collect_substituent_carbons(graph, root_idx, excluded)
    n = len(carbons)
    if n == 0:
        return ""
    if n == 1:
        return "methylidene"
    stem = CHAIN_PREFIX.get(n, f"C{n}")
    return f"{stem}ylidene"


def name_substituent(
    graph: MoleculeGraph,
    root_idx: int,
    excluded: set[int],
) -> str:
    """
    root_idx から始まる置換基の IUPAC 名を返す。

    Args:
        graph: MoleculeGraph
        root_idx: 置換基の付け根原子インデックス
        excluded: 除外する原子インデックスのセット（主鎖・主官能基）

    Returns:
        置換基名 ('methyl', 'chloro', '1-methylethyl' など)
    """
    atom = get_atom(graph, root_idx)

    # ハロゲン
    if atom.symbol in HALOGEN_PREFIX:
        return HALOGEN_PREFIX[atom.symbol]

    # 酸素置換基: -OH → hydroxy, =O → oxo, -O-R → (R)oxy
    if atom.symbol == "O":
        from .molecule_analyzer import get_bond_order as _gbo_o
        neighbors = graph.adjacency[root_idx]
        # =O (exocyclic ketone/aldehyde substituent) → oxo
        for nb in neighbors:
            if _gbo_o(graph, root_idx, nb) == 2.0:
                return "oxo"
        has_h = any(get_atom(graph, nb).symbol == "H" for nb in neighbors)
        if has_h:
            return "hydroxy"
        # エーテル酸素: 主鎖/環外側の炭素を特定して (R)oxy を構築
        c_other = [
            nb for nb in neighbors
            if nb not in excluded and get_atom(graph, nb).symbol == "C"
        ]
        if c_other:
            c_nb = c_other[0]
            # アシルオキシ: C 隣接がカルボニル C (C=O) → (acyl)oxy
            is_acyloxy = any(
                get_atom(graph, oxo).symbol == "O"
                and _gbo_o(graph, c_nb, oxo) == 2.0
                for oxo in graph.adjacency[c_nb]
                if oxo != root_idx
            )
            acyl_or_alkyl = _name_carbon_substituent(graph, c_nb, excluded | {root_idx})
            if is_acyloxy:
                return f"({acyl_or_alkyl}oxy)"
            return _make_oxy_name(acyl_or_alkyl)
        return "oxy"

    # 硫黄置換基: -SH → sulfanyl, -S-R → (R)sulfanyl,
    #   -S(=O)-R → (R)sulfinyl, -S(=O)₂-R → (R)sulfonyl
    if atom.symbol == "S":
        from .molecule_analyzer import get_bond_order as _gbo_s
        neighbors = graph.adjacency[root_idx]
        has_h = any(get_atom(graph, nb).symbol == "H" for nb in neighbors)
        if has_h:
            return "sulfanyl"
        n_oxo = sum(
            1 for nb in neighbors
            if get_atom(graph, nb).symbol == "O"
            and _gbo_s(graph, root_idx, nb) == 2.0
        )
        c_other = [
            nb for nb in neighbors
            if nb not in excluded and get_atom(graph, nb).symbol == "C"
        ]
        # S-C≡N → thiocyanato (Phase 382)
        if c_other and n_oxo == 0:
            _c_nb = c_other[0]
            for _n_nb in graph.adjacency[_c_nb]:
                if _n_nb in excluded or _n_nb == root_idx:
                    continue
                if get_atom(graph, _n_nb).symbol == "N" and _gbo_s(graph, _c_nb, _n_nb) == 3.0:
                    _n_heavy = [x for x in graph.adjacency[_n_nb]
                                if x != _c_nb and get_atom(graph, x).symbol != "H"]
                    if not _n_heavy:
                        return "thiocyanato"
        if c_other:
            alkyl = _name_carbon_substituent(graph, c_other[0], excluded | {root_idx})
            if n_oxo == 2:
                sfx = "sulfonyl"
            elif n_oxo == 1:
                sfx = "sulfinyl"
            else:
                return _make_sulfanyl_name(alkyl)
            if any(c.isdigit() for c in alkyl) or "-" in alkyl:
                return f"({alkyl}){sfx}"
            return f"{alkyl}{sfx}"
        if n_oxo == 2:
            return "sulfonyl"
        if n_oxo == 1:
            return "sulfinyl"
        return "sulfanyl"

    # セレン/テルル置換基: -SeH → selanyl, -TeH → tellanyl, -Se-R → (R)selanyl
    if atom.symbol in ("Se", "Te"):
        neighbors = graph.adjacency[root_idx]
        has_h = any(get_atom(graph, nb).symbol == "H" for nb in neighbors)
        prefix = "selan" if atom.symbol == "Se" else "tellan"
        if has_h:
            return f"{prefix}yl"
        c_other = [
            nb for nb in neighbors
            if nb not in excluded and get_atom(graph, nb).symbol == "C"
        ]
        if c_other:
            alkyl = _name_carbon_substituent(graph, c_other[0], excluded | {root_idx})
            return _make_sulfanyl_name(alkyl).replace("sulfanyl", f"{prefix}yl")
        return f"{prefix}yl"

    # 窒素置換基: amino / nitro など
    if atom.symbol == "N":
        return _name_nitrogen_substituent(graph, root_idx, excluded)

    # 炭素置換基: DFS で最長炭素鎖を探索
    if atom.symbol == "C":
        # ヘテロ芳香族環の付け根 → heteroaryl-yl (ベンゼン環は _name_carbon_substituent で処理)
        if atom.is_aromatic and atom.in_ring:
            ring_has_hetero = any(
                get_atom(graph, a).symbol != "C"
                for rt in graph.ring_atom_sets if root_idx in rt
                for a in rt
            )
            if ring_has_hetero:
                ring_neighbors = [
                    nb for nb in graph.adjacency[root_idx]
                    if nb not in excluded and get_atom(graph, nb).is_aromatic
                ]
                if len(ring_neighbors) >= 2:
                    return _name_aryl_substituent(graph, root_idx, excluded)

        # シアノ基: C≡N (nitrile が置換基として付く場合)
        from .molecule_analyzer import get_bond_order as _gbo
        for nb_idx in graph.adjacency[root_idx]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "N" and nb_idx not in excluded:
                if _gbo(graph, root_idx, nb_idx) == 3.0:
                    # N が他の重原子と結合していないことを確認
                    n_heavy = [
                        n for n in graph.adjacency[nb_idx]
                        if n != root_idx and get_atom(graph, n).symbol != "H"
                    ]
                    if not n_heavy:
                        return "cyano"
        return _name_carbon_substituent(graph, root_idx, excluded)

    return f"({atom.symbol})"  # 未対応元素のフォールバック


def _name_aryl_substituent(graph: "MoleculeGraph", root_idx: int, excluded: set[int]) -> str:
    """
    芳香族環の付け根 root_idx を置換基として命名する。
    benzene → "phenyl"、heteroaryl → "{base}-{locant}-yl"
    縮合ヘテロ芳香族 (1H-indole 等) も対応 (Phase 159)。
    """
    from .ring_handler import _order_ring
    from .heterocycle_handler import _match_retained

    # graph.ring_atom_sets からこの原子を含む環を検索
    ring_tuple: tuple[int, ...] | None = None
    for rt in graph.ring_atom_sets:
        if root_idx in rt:
            ring_tuple = rt
            break

    if ring_tuple is None:
        return "phenyl"

    ring = _order_ring(list(ring_tuple), graph)

    # 純粋ベンゼン環かチェック (単一環の場合)
    if all(get_atom(graph, a).symbol == "C" for a in ring) and len(ring) == 6:
        # 縮合環の場合は phenyl ではなく naphthyl 等の可能性があるが簡略化
        return "phenyl"

    # Phase 159: 縮合ヘテロ芳香族システム (indole, benzimidazole 等)
    # root_idx を含む全リング系を収集
    all_ring_atoms: set[int] = set()
    for rt in (graph.ring_atom_sets or []):
        if root_idx in rt:
            all_ring_atoms.update(rt)
    # 共有原子を持つ他のリングも追加 (縮合環を全収集)
    changed_r = True
    while changed_r:
        changed_r = False
        for rt in (graph.ring_atom_sets or []):
            rt_set = set(rt)
            if rt_set & all_ring_atoms and not rt_set <= all_ring_atoms:
                all_ring_atoms.update(rt_set)
                changed_r = True

    if len(all_ring_atoms) > len(ring):
        # 縮合環系 → _FUSED_HETERO_RETAINED / _FUSED_LOCANT_MAP で命名
        from .heterocycle_handler import _FUSED_HETERO_RETAINED, _FUSED_LOCANT_MAP
        if graph.rdkit_mol is not None:
            from rdkit.Chem import MolToSmiles, MolFromSmiles, MolFragmentToSmiles
            import re as _re_a
            _core_raw = MolFragmentToSmiles(
                graph.rdkit_mol, sorted(all_ring_atoms), canonical=True)
            _tmp_a = MolFromSmiles(_core_raw)
            if _tmp_a is not None:
                _core_smi = MolToSmiles(_tmp_a)
            else:
                _bare_n = [m.start() for m in _re_a.finditer(r"(?<!\[)n(?!\])", _core_raw)]
                _core_smi = _core_raw
                for _pos in _bare_n:
                    _alt = _core_raw[:_pos] + "[nH]" + _core_raw[_pos + 1:]
                    _tmp2 = MolFromSmiles(_alt)
                    if _tmp2 is not None:
                        _alt_c = MolToSmiles(_tmp2)
                        if _alt_c in _FUSED_HETERO_RETAINED:
                            _core_smi = _alt_c
                            break
            _base_name_f = _FUSED_HETERO_RETAINED.get(_core_smi)
            _loc_map_f = _FUSED_LOCANT_MAP.get(_core_smi)
            if _base_name_f is not None and _loc_map_f is not None:
                # substructure match to find locant of root_idx
                _base_mol = MolFromSmiles(_core_smi)
                if _base_mol is not None:
                    _match_f = graph.rdkit_mol.GetSubstructMatch(_base_mol)
                    if _match_f:
                        _rdkit_to_loc: dict[int, int] = {}
                        for _bi, _ri in enumerate(_match_f):
                            _loc = _loc_map_f.get(_bi)
                            if _loc is not None:
                                _rdkit_to_loc[_ri] = _loc
                        _root_loc = _rdkit_to_loc.get(root_idx)
                        if _root_loc is not None:
                            _stem_f = (
                                _base_name_f[:-1] if _base_name_f.endswith("e")
                                else _base_name_f
                            )
                            return f"{_stem_f}-{_root_loc}-yl"

    # 単純ヘテロ芳香族: retained name を取得
    match = _match_retained(ring, graph)
    if match is None:
        return "phenyl"  # fallback
    base_name, _is_nh, rotation = match

    # "1H-" プレフィクスを付与 (is_nh フラグ)
    if _is_nh:
        base_name = f"1H-{base_name}"

    # attachment locant
    if root_idx in rotation:
        locant = rotation.index(root_idx) + 1
    else:
        locant = 1  # fallback

    # "pyridine" → "pyridin-4-yl", "furan" → "furan-2-yl", "thiophene" → "thiophen-2-yl"
    stem = base_name[:-1] if base_name.endswith("e") else base_name
    return f"{stem}-{locant}-yl"


def _count_acyl_chain(graph: "MoleculeGraph", carbonyl_c: int, excluded: set[int]) -> int:
    """
    カルボニル C から始まる acyl 側の炭素鎖の長さを返す (carbonyl_c 自身を含む)。
    O (二重結合) および excluded 方向は無視する。
    """
    from .molecule_analyzer import get_atom as _ga2, get_bond_order as _gbo2

    visited = set(excluded) | {carbonyl_c}
    stack = [carbonyl_c]
    count = 0
    while stack:
        idx = stack.pop()
        if get_atom(graph, idx).symbol != "C":
            continue
        count += 1
        for nb in graph.adjacency[idx]:
            if nb in visited:
                continue
            nb_a = _ga2(graph, nb)
            if nb_a.symbol == "H":
                continue
            # O (=O or -OH) はスキップ
            if nb_a.symbol == "O":
                continue
            # 次の C → 鎖を延ばす
            if nb_a.symbol == "C" and not nb_a.in_ring:
                visited.add(nb)
                stack.append(nb)
    return count


def _name_nitrogen_substituent(
    graph: "MoleculeGraph", n_idx: int, excluded: "set[int] | None" = None
) -> str:
    """
    N 原子から始まる置換基を命名する。

    対応パターン:
      –NH₂              → amino  (第一級アミン)
      –NHR              → {R}amino  (二級アミノ)  (Phase 520)
      –NR₁R₂            → {R1}{R2}amino  (三級アミノ)  (Phase 520)
      –NO₂ / –N⁺(=O)[O⁻] → nitro
      –N=O              → nitroso  (Phase 52)
      –N=[N+]=[N-]      → azido    (Phase 53)
      –N=C=O            → isocyanato (Phase 54)
    """
    from .molecule_analyzer import get_bond_order as _gbo
    _excl = excluded or set()
    neighbors = graph.adjacency[n_idx]
    o_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "O"]
    h_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "H"]
    n_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "N"]
    c_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "C"]

    # ニトロ基: N に O が 2 つ結合 (–NO₂ または –N⁺(=O)[O⁻])
    if len(o_neighbors) == 2:
        return "nitro"

    # ニトロソ基: N に O が 1 つ、二重結合 (–N=O) (Phase 52)
    if len(o_neighbors) == 1 and len(h_neighbors) == 0:
        o_idx = o_neighbors[0]
        if _gbo(graph, n_idx, o_idx) == 2.0:
            return "nitroso"

    # アジド基: N=N+=N- または N=[N+]=[N-] (Phase 53)
    if len(n_neighbors) >= 1 and len(o_neighbors) == 0 and len(h_neighbors) == 0:
        for nn_idx in n_neighbors:
            if _gbo(graph, n_idx, nn_idx) == 2.0:
                nn = get_atom(graph, nn_idx)
                if nn.symbol == "N":
                    # N1=N2, N2 の隣に N3 があるか
                    nn2_neighbors = [nb for nb in graph.adjacency[nn_idx]
                                     if nb != n_idx and get_atom(graph, nb).symbol == "N"]
                    if nn2_neighbors:
                        return "azido"

    # イソシアナト基: N=C=O (Phase 54) / イソチオシアナト基: N=C=S (Phase 162)
    s_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "S"]
    if len(c_neighbors) >= 1 and len(o_neighbors) == 0 and len(h_neighbors) == 0:
        for cn_idx in c_neighbors:
            if _gbo(graph, n_idx, cn_idx) == 2.0:
                for cnb in graph.adjacency[cn_idx]:
                    if cnb == n_idx:
                        continue
                    if get_atom(graph, cnb).symbol == "O" and _gbo(graph, cn_idx, cnb) == 2.0:
                        return "isocyanato"
                    if get_atom(graph, cnb).symbol == "S" and _gbo(graph, cn_idx, cnb) == 2.0:
                        return "isothiocyanato"

    # アミノ基: N に H が 2 つ以上 (–NH₂)
    if len(h_neighbors) >= 2:
        return "amino"

    # Phase 157: アシルアミノ (acetamido 等): N-H 1 つ + C=O 隣接
    if len(h_neighbors) == 1 and c_neighbors:
        from .molecule_analyzer import get_bond_order as _gbo_n
        for _cn in c_neighbors:
            _cn_atom = get_atom(graph, _cn)
            # carbonyl C かチェック
            _has_co = any(
                get_atom(graph, _nb).symbol == "O"
                and _gbo_n(graph, _cn, _nb) == 2.0
                for _nb in graph.adjacency[_cn]
            )
            if not _has_co:
                continue
            # acyl 側の鎖を収集 (N 方向を除外)
            # carbonyl C からの非-O, 非-N 方向の炭素鎖
            _acyl_chain = _count_acyl_chain(graph, _cn, {n_idx})
            # 芳香環に付いている場合 → benzamido 等
            _ring_nb = [
                _nb for _nb in graph.adjacency[_cn]
                if _nb not in {n_idx}
                and get_atom(graph, _nb).symbol == "C"
                and get_atom(graph, _nb).is_aromatic
                and get_atom(graph, _nb).in_ring
            ]
            if _ring_nb:
                _ring_set_a: set[int] = set()
                _q_a = [_ring_nb[0]]
                while _q_a:
                    _a_a = _q_a.pop()
                    if _a_a in _ring_set_a or not get_atom(graph, _a_a).in_ring:
                        continue
                    _ring_set_a.add(_a_a)
                    _q_a.extend(graph.adjacency[_a_a])
                if (len(_ring_set_a) == 6
                        and all(get_atom(graph, _a_a).symbol == "C" for _a_a in _ring_set_a)):
                    return "benzamido"
                return "aroylaminoo"  # ヘテロ芳香環 (fallback)
            from .constants import CHAIN_PREFIX as _CP_n
            if _acyl_chain == 1:
                return "formamido"
            elif _acyl_chain == 2:
                return "acetamido"
            else:
                _stem_n = _CP_n.get(_acyl_chain, f"C{_acyl_chain}")
                return f"{_stem_n}anamido"

    # シアノ基の末端 N: C≡N の N 末端 → cyano (Phase 383)
    if (len(c_neighbors) == 1 and len(o_neighbors) == 0
            and len(h_neighbors) == 0 and len(n_neighbors) == 0):
        if _gbo(graph, n_idx, c_neighbors[0]) == 3.0:
            return "cyano"

    # Phase 520: 二級/三級アミノ置換基: N-H + R → {R}amino; N(R)(R') → {R,R'}amino
    _r_cs = [c for c in c_neighbors if c not in _excl]
    if _r_cs:
        from .constants import MULTIPLIER as _MULT_N
        _r_names = sorted(_name_carbon_substituent(graph, c, {n_idx}) for c in _r_cs)
        if len(set(_r_names)) == 1:
            _mult_n = _MULT_N.get(len(_r_names), "")
            _rn = _r_names[0]
            _rn_fmt = f"({_rn})" if any(ch.isdigit() for ch in _rn) or "-" in _rn else _rn
            return f"{_mult_n}{_rn_fmt}amino"
        _parts_n = [
            f"({n})" if any(ch.isdigit() for ch in n) or "-" in n else n
            for n in _r_names
        ]
        return "(" + ",".join(_parts_n) + ")amino"

    # N-H のみ (primary amino fallback)
    if len(h_neighbors) >= 1:
        return "amino"

    return "(N)"  # 未対応


def _name_carbon_substituent(
    graph: MoleculeGraph,
    root_idx: int,
    excluded: set[int],
) -> str:
    """
    炭素置換基の命名。
    単純鎖: methyl, ethyl, propyl, butyl ... (付け根が末端炭素の直鎖)
    分岐: 1-methylethyl など (複合置換基、付け根が非末端 or 内部分岐あり)
    """
    from .molecule_analyzer import get_bond_order as _gbo

    # カルボキシ基 (root が COOH or COO-): carboxy prefix  (Phase 148)
    _o_double = [nb for nb in graph.adjacency[root_idx]
                 if get_atom(graph, nb).symbol == "O" and _gbo(graph, root_idx, nb) == 2.0]
    _o_single_oh = [nb for nb in graph.adjacency[root_idx]
                    if get_atom(graph, nb).symbol == "O"
                    and _gbo(graph, root_idx, nb) == 1.0
                    and (get_atom(graph, nb).num_hs >= 1
                         or any(get_atom(graph, hh).symbol == "H"
                                for hh in graph.adjacency[nb]))]
    _o_single_neg = [nb for nb in graph.adjacency[root_idx]
                     if get_atom(graph, nb).symbol == "O"
                     and _gbo(graph, root_idx, nb) == 1.0
                     and get_atom(graph, nb).formal_charge == -1]
    if _o_double and (_o_single_oh or _o_single_neg):
        return "carboxy"

    # アシル基 (root が C=O): ethanoyl, propanoyl 等 (Phase 49)
    for nb_idx in graph.adjacency[root_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "O" and _gbo(graph, root_idx, nb_idx) == 2.0:
            # root が カルボニル C → acyl 置換基
            # 芳香族環に直接結合したカルボニル → benzoyl (Phase 265)
            for aryl_nb in graph.adjacency[root_idx]:
                if aryl_nb in excluded or aryl_nb == nb_idx:
                    continue
                aryl_atom = get_atom(graph, aryl_nb)
                if aryl_atom.symbol == "C" and aryl_atom.is_aromatic and aryl_atom.in_ring:
                    ring_atoms = next(
                        (set(rt) for rt in (graph.ring_atom_sets or []) if aryl_nb in rt),
                        set(),
                    )
                    if (len(ring_atoms) == 6
                            and all(get_atom(graph, a).symbol == "C" for a in ring_atoms)):
                        return "benzoyl"
            acyl_carbons = _collect_substituent_carbons(graph, root_idx, excluded)
            n_acyl = len(acyl_carbons)
            if n_acyl == 1:
                return "formyl"
            if n_acyl == 2:
                return "acetyl"
            stem = CHAIN_PREFIX.get(n_acyl, f"C{n_acyl}")
            return f"{stem}anoyl"

    # 置換基の炭素を DFS で列挙
    sub_carbons = _collect_substituent_carbons(graph, root_idx, excluded)
    n = len(sub_carbons)

    if n == 0:
        return "methyl"  # フォールバック

    carbon_set = set(sub_carbons)

    # アリール基 (root が芳香環 C): phenyl or heteroaryl-yl (Phase 50)
    root_atom = get_atom(graph, root_idx)
    if root_atom.in_ring and root_atom.is_aromatic:
        # ヘテロ芳香族環チェック
        ring_has_hetero = any(
            get_atom(graph, a).symbol != "C"
            for rt in graph.ring_atom_sets if root_idx in rt
            for a in rt
        )
        if ring_has_hetero:
            return _name_aryl_substituent(graph, root_idx, excluded)

        # root から到達できる芳香族炭素環 C を収集
        aryl_cs: set[int] = set()
        q = [root_idx]
        while q:
            a = q.pop()
            if a in aryl_cs or a in excluded:
                continue
            atom_a = get_atom(graph, a)
            if not (atom_a.in_ring and atom_a.is_aromatic and atom_a.symbol == "C"):
                continue
            aryl_cs.add(a)
            q.extend(graph.adjacency[a])
        if len(aryl_cs) != 6:
            return "phenyl"  # フォールバック
        # 環外に非水素置換基があるか確認
        has_external = any(
            get_atom(graph, nb).symbol != "H"
            for ring_c in aryl_cs
            for nb in graph.adjacency[ring_c]
            if nb not in excluded and nb not in aryl_cs
        )
        if not has_external:
            return "phenyl"
        # 置換フェニル: 環外置換基を検出してロカントを付ける
        # 環の走査順序を決定 (root_idx = locant 1)
        def _order_benzene(start: int, ring_set: set[int]) -> list[int]:
            order = [start]
            prev = -1
            curr = start
            while len(order) < 6:
                nexts = [nb for nb in graph.adjacency[curr]
                         if nb in ring_set and nb != prev and nb not in order]
                if not nexts:
                    break
                prev, curr = curr, nexts[0]
                order.append(curr)
            return order
        order_fwd = _order_benzene(root_idx, aryl_cs)
        order_rev = [root_idx] + list(reversed(order_fwd[1:]))
        # 環外置換基を全種類 (C/halogen/O/N 等) 収集し最小ロカントの方向を選択
        best_subs: list[tuple[int, int]] | None = None  # [(locant, sub_atom), ...]
        best_locs2: list[int] | None = None
        for order in (order_fwd, order_rev):
            loc_map2 = {idx: i + 1 for i, idx in enumerate(order)}
            ring_subs2: list[tuple[int, int]] = []
            for ring_c in aryl_cs:
                loc = loc_map2[ring_c]
                for nb in graph.adjacency[ring_c]:
                    if nb in excluded or nb in aryl_cs:
                        continue
                    if get_atom(graph, nb).symbol != "H":
                        ring_subs2.append((loc, nb))
            cur_locs = sorted(t[0] for t in ring_subs2)
            if best_locs2 is None or cur_locs < best_locs2:
                best_locs2 = cur_locs
                best_subs = ring_subs2
        if not best_subs:
            return "phenyl"
        # 置換基名の組み立て (name_substituent で全種類に対応)
        from collections import defaultdict
        grouped2: dict[str, list[int]] = defaultdict(list)
        for loc, sub_atom in best_subs:
            sub_nm = name_substituent(graph, sub_atom, excluded | aryl_cs)
            grouped2[sub_nm].append(loc)
        prefix_parts2: list[str] = []
        for nm in sorted(grouped2):
            from .name_assembler import _needs_bis_tris as _nbt2
            locs = sorted(grouped2[nm])
            loc_str = ",".join(str(l) for l in locs)
            n_sub = len(locs)
            if n_sub > 1 and _nbt2(nm):
                from .constants import MULTIPLIER as _MULT2
                from .name_assembler import _BIS_MULTIPLIER as _BMULT2
                mult = _BMULT2.get(n_sub, "")
                prefix_parts2.append(f"{loc_str}-{mult}({nm})")
            elif _nbt2(nm):
                prefix_parts2.append(f"{loc_str}-({nm})")
            elif n_sub > 1:
                from .constants import MULTIPLIER as _MULT2
                mult = _MULT2.get(n_sub, "")
                prefix_parts2.append(f"{loc_str}-{mult}{nm}")
            else:
                prefix_parts2.append(f"{loc_str}-{nm}")
        return "-".join(prefix_parts2) + "phenyl"

    # シクロアルキル / ヘテロシクロ: root が環内かつ置換基全体が単一環 (Phase 49)
    if root_atom.in_ring and not root_atom.is_aromatic:
        # root を含む環の全原子を収集 (ヘテロ原子も含む)
        all_ring_set: set[int] | None = None
        for rt in (graph.ring_atom_sets or []):
            if root_idx in rt:
                all_ring_set = set(rt)
                break
        # 環内 C のみも収集
        ring_cs: set[int] = set()
        q = [root_idx]
        while q:
            a = q.pop()
            if a in ring_cs or a in excluded:
                continue
            if not get_atom(graph, a).in_ring:
                continue
            if get_atom(graph, a).symbol != "C":
                continue
            ring_cs.add(a)
            q.extend(graph.adjacency[a])
        ring_size = len(ring_cs)
        # 環外 C が存在しない (sub_carbons == ring_cs) → 純シクロアルキル or ヘテロ環
        if carbon_set == ring_cs and all_ring_set is not None:
            # ヘテロ原子が含まれる場合は ヘテロ環 substituent 命名 (Phase 247)
            hetero_in_ring = [a for a in all_ring_set if get_atom(graph, a).symbol != "C"]
            if hetero_in_ring:
                from .ring_handler import _order_ring
                from .heterocycle_handler import _match_hantzsch_widman, _find_best_start
                ring_ordered = _order_ring(list(all_ring_set), graph)
                hw_name = _match_hantzsch_widman(ring_ordered, graph)
                if hw_name is not None:
                    # ロカントを決定: ヘテロ原子が locant 1 になるよう並べる
                    best_start = _find_best_start(ring_ordered, graph)
                    try:
                        locant = best_start.index(root_idx) + 1
                    except ValueError:
                        locant = ring_ordered.index(root_idx) + 1
                    stem_hw = hw_name[:-1] if hw_name.endswith("e") else hw_name
                    return f"{stem_hw}-{locant}-yl"
                # 保留名テーブルを試みる (morpholine, piperazine 等)
                from .heterocycle_handler import _RETAINED_NAMES as _RN_H
                _try_sigs = []
                for _start_idx in range(len(ring_ordered)):
                    _rotated = ring_ordered[_start_idx:] + ring_ordered[:_start_idx]
                    _sig = tuple(get_atom(graph, a).symbol for a in _rotated)
                    _try_sigs.append((False, _sig, _rotated))
                for _arom_k, _sig_k, _rot_k in _try_sigs:
                    _h_match = _RN_H.get((_arom_k, _sig_k))
                    if _h_match is not None:
                        _base, _has_nh = _h_match[0], _h_match[1]
                        _loc = _rot_k.index(root_idx) + 1
                        _stem = _base[:-1] if _base.endswith("e") else _base
                        return f"{_stem}-{_loc}-yl"
            else:
                stem = CHAIN_PREFIX.get(ring_size, f"C{ring_size}")
                return f"cyclo{stem}yl"

    # 付け根 (root) が末端かどうか確認
    # 末端 = 置換基内での炭素隣接が 1 以下
    root_c_in_sub = [
        nb for nb in graph.adjacency[root_idx]
        if get_atom(graph, nb).symbol == "C" and nb in carbon_set
    ]
    root_is_terminal = len(root_c_in_sub) <= 1

    # 単純 n-アルキル: 付け根が末端 かつ 全体が直鎖状
    if root_is_terminal and _is_linear_chain(graph, sub_carbons, excluded):
        prefix = CHAIN_PREFIX.get(n, f"C{n}")
        chain_path = _find_longest_path(graph, root_idx, excluded, carbon_set)
        # 鎖末端が COOH の場合: carboxymethyl / 2-carboxyethyl 等 (Phase 161)
        from .molecule_analyzer import get_bond_order as _gbo_co
        last_c = chain_path[-1]
        _last_o_dbl = [nb for nb in graph.adjacency[last_c]
                       if nb not in excluded and nb not in carbon_set
                       and get_atom(graph, nb).symbol == "O"
                       and _gbo_co(graph, last_c, nb) == 2.0]
        _last_o_oh = [nb for nb in graph.adjacency[last_c]
                      if nb not in excluded and nb not in carbon_set
                      and get_atom(graph, nb).symbol == "O"
                      and _gbo_co(graph, last_c, nb) == 1.0
                      and (get_atom(graph, nb).num_hs >= 1
                           or any(get_atom(graph, hh).symbol == "H"
                                  for hh in graph.adjacency[nb])
                           or get_atom(graph, nb).formal_charge == -1)]
        if _last_o_dbl and _last_o_oh and n >= 2:
            alkyl_n = n - 1  # COOH C を除いたアルキル鎖の長さ
            alkyl_pfx = CHAIN_PREFIX.get(alkyl_n, f"C{alkyl_n}")
            if alkyl_n == 1:
                return "carboxymethyl"
            else:
                return f"{alkyl_n}-carboxy{alkyl_pfx}yl"

        # 鎖端に芳香環が付いている場合は分岐置換基パスへ (Phase 160)
        has_aryl_nb = any(
            get_atom(graph, nb).symbol == "C"
            and get_atom(graph, nb).in_ring
            and get_atom(graph, nb).is_aromatic
            for c_idx in chain_path
            for nb in graph.adjacency[c_idx]
            if nb not in excluded and nb not in carbon_set
            and get_atom(graph, nb).symbol != "H"
        )
        if has_aryl_nb:
            return _name_branched_substituent(graph, root_idx, excluded, sub_carbons)
        # アルコキシ置換: 鎖内 C に ether O(-C) が付いている場合 → _name_branched_substituent へ
        has_oxy_sub = False
        for _c_idx in chain_path:
            for _nb in graph.adjacency[_c_idx]:
                if _nb in carbon_set or _nb in excluded:
                    continue
                if get_atom(graph, _nb).symbol != "O":
                    continue
                _o_c_nbs = [
                    nb2 for nb2 in graph.adjacency[_nb]
                    if nb2 != _c_idx and nb2 not in excluded
                    and get_atom(graph, nb2).symbol == "C"
                ]
                if _o_c_nbs:
                    has_oxy_sub = True
                    break
            if has_oxy_sub:
                break
        if has_oxy_sub:
            return _name_branched_substituent(graph, root_idx, excluded, sub_carbons)
        # チオエーテル (-S-C) が鎖内にある場合: _name_branched_substituent へ
        has_thioether_sub = False
        for _c_idx in chain_path:
            for _nb in graph.adjacency[_c_idx]:
                if _nb in carbon_set or _nb in excluded:
                    continue
                if get_atom(graph, _nb).symbol != "S":
                    continue
                _s_c_nbs = [nb2 for nb2 in graph.adjacency[_nb]
                            if nb2 != _c_idx and nb2 not in excluded
                            and get_atom(graph, nb2).symbol == "C"]
                if _s_c_nbs:
                    has_thioether_sub = True
                    break
            if has_thioether_sub:
                break
        if has_thioether_sub:
            return _name_branched_substituent(graph, root_idx, excluded, sub_carbons)
        # 多重結合 (二重・三重) が鎖内にあるか確認
        for i in range(len(chain_path) - 1):
            bo = _gbo(graph, chain_path[i], chain_path[i + 1])
            if bo == 3.0:
                # C≡C: n=2 → "ethynyl"; n>=3 → "prop-1-yn-1-yl" 等
                loc = i + 1
                if n == 2:
                    return f"{prefix}ynyl"
                return f"{prefix}-{loc}-yn-1-yl"
            if bo == 2.0:
                # C=C: n=2 → "ethenyl"; n>=3 → "prop-2-en-1-yl" 等
                loc = i + 1
                stereo_pfx_sub = ""
                try:
                    from .stereochemistry import assign_stereochemistry
                    from .chain_finder import PrincipalChain
                    _pc_sub = PrincipalChain(
                        atom_indices=chain_path,
                        locant_map={c: j + 1 for j, c in enumerate(chain_path)})
                    _stereo_sub = assign_stereochemistry(graph, _pc_sub)
                    if _stereo_sub:
                        stereo_pfx_sub = "(" + ",".join(d.strip("()") for d in _stereo_sub) + ")-"
                except Exception:
                    pass
                if n == 2:
                    return f"{stereo_pfx_sub}{prefix}enyl"
                return f"{stereo_pfx_sub}{prefix}-{loc}-en-1-yl"
        # ハロゲン置換アルキル: CF3 → "trifluoromethyl" 等 (IUPAC P-61.7)
        _HAL_PREFIX = {"F": "fluoro", "Cl": "chloro", "Br": "bromo", "I": "iodo"}
        hal_subs: list[tuple[int, str]] = []
        for pos, c_idx in enumerate(chain_path, 1):
            for nb_idx in graph.adjacency[c_idx]:
                if nb_idx in excluded or nb_idx in carbon_set:
                    continue
                nb_sym = get_atom(graph, nb_idx).symbol
                if nb_sym in _HAL_PREFIX:
                    hal_subs.append((pos, _HAL_PREFIX[nb_sym]))
        if hal_subs:
            from collections import defaultdict
            hal_by_name: dict[str, list[int]] = defaultdict(list)
            for pos, nm in hal_subs:
                hal_by_name[nm].append(pos)
            hal_parts: list[str] = []
            for nm in sorted(hal_by_name.keys()):
                locs = sorted(hal_by_name[nm])
                mult = MULTIPLIER.get(len(locs), "")
                if n == 1:
                    hal_parts.append(f"{mult}{nm}")
                else:
                    loc_str = ",".join(str(l) for l in locs)
                    hal_parts.append(f"{loc_str}-{mult}{nm}")
            return f"{'-'.join(hal_parts)}{prefix}yl"

        # ヒドロキシ置換アルキル: HO-CH₂- → hydroxymethyl, HO-CH₂CH₂- → 2-hydroxyethyl (Phase 274)
        hy_subs: list[tuple[int, str]] = []
        for pos, c_idx in enumerate(chain_path, 1):
            for nb_idx in graph.adjacency[c_idx]:
                if nb_idx in excluded or nb_idx in carbon_set:
                    continue
                nb_atom_hy = get_atom(graph, nb_idx)
                if nb_atom_hy.symbol != "O":
                    continue
                # 自由 OH: H を持ち、他 C に結合していない (エーテルでない)
                hy_has_h = any(get_atom(graph, hh).symbol == "H"
                               for hh in graph.adjacency[nb_idx])
                hy_c_nbs = [nb2 for nb2 in graph.adjacency[nb_idx]
                            if nb2 != c_idx and nb2 not in excluded
                            and get_atom(graph, nb2).symbol == "C"]
                # C=O でない (シングル結合) かつ OH かつ エーテルでない
                from .molecule_analyzer import get_bond_order as _gbo_hy
                if (_gbo_hy(graph, c_idx, nb_idx) == 1.0
                        and hy_has_h and not hy_c_nbs):
                    hy_subs.append((pos, "hydroxy"))
        if hy_subs:
            from collections import defaultdict as _dd_hy
            hy_by_name: dict[str, list[int]] = _dd_hy(list)
            for pos, nm in hy_subs:
                hy_by_name[nm].append(pos)
            hy_parts: list[str] = []
            for nm in sorted(hy_by_name.keys()):
                locs = sorted(hy_by_name[nm])
                mult = MULTIPLIER.get(len(locs), "")
                if n == 1:
                    hy_parts.append(f"{mult}{nm}")
                else:
                    loc_str = ",".join(str(l) for l in locs)
                    hy_parts.append(f"{loc_str}-{mult}{nm}")
            return f"{'-'.join(hy_parts)}{prefix}yl"

        # アミノ置換アルキル: H₂N-CH₂- → aminomethyl, H₂N-CH₂CH₂- → 2-aminoethyl 等 (Phase 216)
        # 鎖上の N 置換基を探す
        amino_subs: list[tuple[int, str]] = []
        for pos, c_idx in enumerate(chain_path, 1):
            for nb_idx in graph.adjacency[c_idx]:
                if nb_idx in excluded or nb_idx in carbon_set:
                    continue
                if get_atom(graph, nb_idx).symbol != "N":
                    continue
                n_atom = get_atom(graph, nb_idx)
                # N が環内・非芳香族の場合はスキップ
                if n_atom.in_ring:
                    continue
                # N に付いている H 数と C 数を確認
                n_h = sum(1 for hh in graph.adjacency[nb_idx] if get_atom(graph, hh).symbol == "H")
                n_c = [cc for cc in graph.adjacency[nb_idx]
                       if cc != c_idx and get_atom(graph, cc).symbol == "C"]
                # 三重結合 N (ニトリル/イソニトリル) はアミノ扱いしない
                if _gbo(graph, c_idx, nb_idx) == 3.0:
                    continue
                if not n_c and n_h >= 0:  # primary amine (NH2 or NH)
                    amino_subs.append((pos, "amino"))
        if amino_subs:
            from collections import defaultdict
            amino_by_name: dict[str, list[int]] = defaultdict(list)
            for pos, nm in amino_subs:
                amino_by_name[nm].append(pos)
            amino_parts: list[str] = []
            for nm in sorted(amino_by_name.keys()):
                locs = sorted(amino_by_name[nm])
                mult = MULTIPLIER.get(len(locs), "")
                if n == 1:
                    amino_parts.append(f"{mult}{nm}")
                else:
                    loc_str = ",".join(str(l) for l in locs)
                    amino_parts.append(f"{loc_str}-{mult}{nm}")
            return f"{'-'.join(amino_parts)}{prefix}yl"

        # スルファニル置換アルキル: HS-CH2- → sulfanylmethyl 等 (Phase 514)
        sulfanyl_subs: list[tuple[int, str]] = []
        for pos, c_idx in enumerate(chain_path, 1):
            for nb_idx in graph.adjacency[c_idx]:
                if nb_idx in excluded or nb_idx in carbon_set:
                    continue
                nb_atom_s = get_atom(graph, nb_idx)
                if nb_atom_s.symbol != "S":
                    continue
                s_has_h = any(get_atom(graph, hh).symbol == "H"
                              for hh in graph.adjacency[nb_idx])
                s_c_nbs = [nb2 for nb2 in graph.adjacency[nb_idx]
                           if nb2 != c_idx and nb2 not in excluded
                           and get_atom(graph, nb2).symbol == "C"]
                if s_has_h and not s_c_nbs:
                    sulfanyl_subs.append((pos, "sulfanyl"))
        if sulfanyl_subs:
            from collections import defaultdict as _dd_s
            s_by_name: dict[str, list[int]] = _dd_s(list)
            for pos, nm in sulfanyl_subs:
                s_by_name[nm].append(pos)
            s_parts: list[str] = []
            for nm in sorted(s_by_name.keys()):
                locs = sorted(s_by_name[nm])
                mult = MULTIPLIER.get(len(locs), "")
                if n == 1:
                    s_parts.append(f"{mult}{nm}")
                else:
                    loc_str = ",".join(str(l) for l in locs)
                    s_parts.append(f"{loc_str}-{mult}{nm}")
            return f"{'-'.join(s_parts)}{prefix}yl"

        return f"{prefix}yl"

    # 分岐置換基 (isopropyl 等): 再帰的に命名
    return _name_branched_substituent(graph, root_idx, excluded, sub_carbons)


def _collect_substituent_carbons(
    graph: MoleculeGraph,
    root_idx: int,
    excluded: set[int],
) -> list[int]:
    """root_idx から到達できる炭素原子を DFS で全収集。

    芳香族環への侵入は、起点が非芳香族の場合に停止する。
    これにより鎖から始まる芳香環への誤った連続カウントを防ぐ。
    """
    visited: set[int] = set(excluded)
    result: list[int] = []
    root_atom = get_atom(graph, root_idx)
    root_is_aromatic = root_atom.in_ring and root_atom.is_aromatic

    def dfs(idx: int, from_aromatic: bool) -> None:
        if idx in visited:
            return
        atom = get_atom(graph, idx)
        if atom.symbol != "C":
            return
        this_is_aromatic = atom.in_ring and atom.is_aromatic
        # 非芳香族領域から芳香族環に入ったら、そこで停止 (環原子は含めない)
        if this_is_aromatic and not from_aromatic:
            return
        visited.add(idx)
        result.append(idx)
        for nb in graph.adjacency[idx]:
            dfs(nb, this_is_aromatic)

    dfs(root_idx, root_is_aromatic)
    return result


def _is_linear_chain(
    graph: MoleculeGraph,
    carbons: list[int],
    excluded: set[int],
) -> bool:
    """炭素リストが直鎖（分岐なし）かを判定。"""
    carbon_set = set(carbons)
    for c_idx in carbons:
        c_neighbors = [
            nb for nb in graph.adjacency[c_idx]
            if get_atom(graph, nb).symbol == "C"
            and nb in carbon_set
        ]
        # 末端は1つ、中間は2つの炭素隣接 → これ以上は分岐
        if len(c_neighbors) > 2:
            return False
    return True


def _name_branched_substituent(
    graph: MoleculeGraph,
    root_idx: int,
    excluded: set[int],
    sub_carbons: list[int],
) -> str:
    """
    分岐置換基の IUPAC 命名 (IUPAC 2013 P-31.1.3.4)。
    root が内部炭素の場合は propan-2-yl 形式を使用する。
    例: -CH(CH3)2 → 'propan-2-yl', -CH2CH(CH3)2 → '2-methylpropyl'
    """
    carbon_set = set(sub_carbons)

    # root を通る最長鎖を探索 (through-path)
    main_path = _find_through_path(graph, root_idx, excluded, carbon_set)

    n = len(main_path)
    if n == 0:
        return "methyl"

    root_pos = main_path.index(root_idx) + 1  # 1-indexed position of root

    prefix = CHAIN_PREFIX.get(n, f"C{n}")
    main_set = set(main_path)
    main_locant = {ai: loc for loc, ai in enumerate(main_path, start=1)}

    # 分岐置換基を収集
    sub_result: list[tuple[int, str]] = []
    new_excluded = excluded | main_set

    for c_idx in main_path:
        loc = main_locant[c_idx]
        for nb in graph.adjacency[c_idx]:
            nb_atom = get_atom(graph, nb)
            if nb in new_excluded or nb_atom.symbol == "H":
                continue
            sub_name = name_substituent(graph, nb, new_excluded)
            sub_result.append((loc, sub_name))

    by_name: dict[str, list[int]] = {}
    for loc, sname in sub_result:
        by_name.setdefault(sname, []).append(loc)

    parts: list[str] = []
    for sname in sorted(by_name.keys()):
        locs = sorted(by_name[sname])
        mult = MULTIPLIER.get(len(locs), "")
        # n==1 の場合はロカントを省略 (唯一の位置のため)
        if n == 1:
            parts.append(f"{mult}{sname}")
        else:
            loc_str = ",".join(str(l) for l in locs)
            parts.append(f"{loc_str}-{mult}{sname}")

    substituent_prefix = "-".join(parts)

    # 主鎖内の二重・三重結合を検出
    from .molecule_analyzer import get_bond_order as _gbo_b
    unsaturation = ""
    for i in range(len(main_path) - 1):
        bo = _gbo_b(graph, main_path[i], main_path[i + 1])
        if bo == 3.0:
            loc_u = i + 1
            unsaturation = f"-{loc_u}-yn" if n >= 3 else "yn"
            break
        if bo == 2.0:
            loc_u = i + 1
            unsaturation = f"-{loc_u}-en" if n >= 3 else "en"
            break

    # ロカントなし不飽和 (n<3): "en"/"yn" → suffix は直接 "yl"
    unsat_has_locant = any(c.isdigit() for c in unsaturation)

    if root_pos == 1:
        # root が鎖の末端: 従来の {prefix}{sub}yl 形式
        if not substituent_prefix:
            return f"{prefix}{unsaturation}yl"
        if unsaturation:
            if unsat_has_locant:
                return f"{substituent_prefix}{prefix}{unsaturation}-1-yl"
            return f"{substituent_prefix}{prefix}{unsaturation}yl"
        return f"{substituent_prefix}{prefix}yl"
    else:
        # root が内部炭素: propan-2-yl 形式 (IUPAC 2013 PIN)
        if not substituent_prefix:
            if unsaturation:
                return f"{prefix}an{unsaturation}-{root_pos}-yl"
            return f"{prefix}an-{root_pos}-yl"
        if unsaturation:
            return f"{substituent_prefix}{prefix}an{unsaturation}-{root_pos}-yl"
        return f"{substituent_prefix}{prefix}an-{root_pos}-yl"


def _find_through_path(
    graph: MoleculeGraph,
    root_idx: int,
    excluded: set[int],
    carbon_set: set[int],
) -> list[int]:
    """
    root を通る最長鎖を返す。
    root が内部炭素の場合は両方向の鎖を結合する。
    """
    # root からの各方向への最長経路を収集
    branches: list[list[int]] = []
    for nb in graph.adjacency[root_idx]:
        nb_atom = get_atom(graph, nb)
        if nb_atom.symbol != "C" or nb in excluded or nb not in carbon_set:
            continue
        tail = _find_longest_path(graph, nb, excluded | {root_idx}, carbon_set)
        branches.append(tail)

    if not branches:
        return [root_idx]

    # 最長と2番目の枝を選ぶ
    branches.sort(key=len, reverse=True)
    longest = branches[0]

    if len(branches) >= 2:
        second = branches[1]
        # through-path: second_reversed + [root] + longest
        return list(reversed(second)) + [root_idx] + longest
    else:
        # root は末端
        return [root_idx] + longest


def _find_longest_path(
    graph: MoleculeGraph,
    start: int,
    excluded: set[int],
    carbon_set: set[int],
) -> list[int]:
    """start から carbon_set 内を DFS して最長炭素経路を返す。"""
    best: list[int] = []

    def dfs(cur: int, visited: set[int], path: list[int]) -> None:
        nonlocal best
        path.append(cur)
        visited.add(cur)
        extended = False
        for nb in graph.adjacency[cur]:
            nb_atom = get_atom(graph, nb)
            if nb_atom.symbol == "C" and nb not in visited and nb not in excluded and nb in carbon_set:
                extended = True
                dfs(nb, visited, path)
        if not extended:
            if len(path) > len(best):
                best = list(path)
        path.pop()
        visited.remove(cur)

    dfs(start, set(excluded), [])
    return best


# ─── ヘテロ原子置換基ヘルパー ──────────────────────────────────────────

def _make_oxy_name(alkyl: str) -> str:
    """
    alkyl 名から oxy 接頭辞を構築する。
    例: 'methyl' → 'methoxy', 'ethyl' → 'ethoxy',
        '1-methylethyl' → '(1-methylethoxy)',
        'propanoyl' → 'propanoyloxy' (acyl groups keep the 'yl')
    """
    if alkyl.endswith("yl"):
        base = alkyl[:-2]
        # acyl groups end in 'o' before 'yl' (propanoyl, butanoyl, benzoyl)
        # → keep the 'yl' in the oxy name: "propanoyloxy" not "propanooxy"
        if base.endswith("o"):
            oxy = alkyl + "oxy"
        else:
            oxy = base + "oxy"
        # ロカント(数字)またはハイフンを含む場合は括弧で囲む
        if any(c.isdigit() for c in base) or "-" in base:
            return f"({oxy})"
        return oxy
    return f"({alkyl})oxy"


def _make_sulfanyl_name(alkyl: str) -> str:
    """
    alkyl 名から sulfanyl 接頭辞を構築する。
    例: 'methyl' → 'methylsulfanyl', '1-methylethyl' → '(1-methylethyl)sulfanyl'
    """
    if any(c.isdigit() for c in alkyl) or "-" in alkyl:
        return f"({alkyl})sulfanyl"
    return alkyl + "sulfanyl"
