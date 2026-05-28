"""
官能基の検出と優先順位付け。

IUPAC 2013 Blue Book P-65 の seniority order に従い、
分子グラフから官能基を検出して優先順に返す。
"""

from __future__ import annotations

from dataclasses import dataclass

from .constants import FUNCTIONAL_GROUP_PRIORITY
from .molecule_analyzer import MoleculeGraph, get_bond_order, get_atom


@dataclass
class FunctionalGroup:
    group_type: str        # 'carboxylic_acid', 'alcohol', ...
    atom_indices: list[int]  # 官能基を構成する原子インデックス
    priority: int          # FUNCTIONAL_GROUP_PRIORITY の値


def detect_groups(graph: MoleculeGraph) -> list[FunctionalGroup]:
    """
    MoleculeGraph を走査し、全官能基を検出して優先順位の高い順に返す。
    最初の要素が principal characteristic group になる。

    ハロゲン (F, Cl, Br, I) は置換基として扱い、ここには含めない。
    芳香族化合物は現在未対応 (NotImplementedError)。
    """
    groups: list[FunctionalGroup] = []

    # 芳香族・環状分子の場合はアルケン/アルキン検出をスキップ
    # (ring_handler.py が環状命名を担当)
    has_aromatic = any(a.is_aromatic for a in graph.atoms)

    # --- 官能基パターン検出 ---
    # 炭素原子のリストを走査
    for atom in graph.atoms:
        if atom.symbol != "C":
            continue
        idx = atom.idx

        # 酸無水物: C(=O)-O-C(=O) (priority=88) — ester より先に判定
        if _is_anhydride(graph, idx):
            o_link: int | None = None
            c2: int | None = None
            for nb_idx in graph.adjacency[idx]:
                nb = get_atom(graph, nb_idx)
                if nb.symbol == "O" and get_bond_order(graph, idx, nb_idx) == 1.0:
                    for o_nb_idx in graph.adjacency[nb_idx]:
                        if o_nb_idx == idx:
                            continue
                        o_nb = get_atom(graph, o_nb_idx)
                        if o_nb.symbol == "C" and _get_double_bonded_oxygen(graph, o_nb_idx) is not None:
                            o_link = nb_idx
                            c2 = o_nb_idx
                            break
                if o_link is not None:
                    break
            indices = [idx]
            if o_link is not None:
                indices.append(o_link)
            if c2 is not None:
                indices.append(c2)
            groups.append(FunctionalGroup(
                group_type="anhydride",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["anhydride"],
            ))
            continue

        # 炭酸エステル: RO-C(=O)-OR (priority=91, Phase 39) — エステルより先に判定
        if _is_carbonate(graph, idx):
            o_idx = _get_double_bonded_oxygen(graph, idx)
            indices = [idx] + ([o_idx] if o_idx is not None else [])
            for nb_idx in graph.adjacency[idx]:
                nb = get_atom(graph, nb_idx)
                if nb.symbol == "O" and get_bond_order(graph, idx, nb_idx) == 1.0:
                    indices.append(nb_idx)
            groups.append(FunctionalGroup(
                group_type="carbonate",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["carbonate"],
            ))
            continue

        # カルバメート: N-C(=O)-O-R (priority=89, Phase 42) — ester/amide より先に判定
        if _is_carbamate(graph, idx):
            o_idx = _get_double_bonded_oxygen(graph, idx)
            indices = [idx] + ([o_idx] if o_idx is not None else [])
            for nb_idx in graph.adjacency[idx]:
                nb = get_atom(graph, nb_idx)
                if nb.symbol in ("O", "N") and get_bond_order(graph, idx, nb_idx) == 1.0:
                    indices.append(nb_idx)
            groups.append(FunctionalGroup(
                group_type="carbamate",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["carbamate"],
            ))
            continue

        # クロロホルメート: Cl-C(=O)-O-R (priority=86, Phase 60) — esterより先に判定
        if _is_chloroformate(graph, idx):
            o_idx = _get_double_bonded_oxygen(graph, idx)
            indices = [idx] + ([o_idx] if o_idx is not None else [])
            for nb_idx in graph.adjacency[idx]:
                nb = get_atom(graph, nb_idx)
                if nb.symbol == "O" and get_bond_order(graph, idx, nb_idx) == 1.0:
                    indices.append(nb_idx)
                    break
            groups.append(FunctionalGroup(
                group_type="chloroformate",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["chloroformate"],
            ))
            continue

        # ペルオキシ酸: C(=O)-O-O-H (priority=99, Phase 77) — ester より先に判定
        if _is_peroxyacid(graph, idx):
            o_double = _get_double_bonded_oxygen(graph, idx)
            indices = [idx] + ([o_double] if o_double is not None else [])
            for nb_idx in graph.adjacency[idx]:
                nb = get_atom(graph, nb_idx)
                if nb.symbol == "O" and get_bond_order(graph, idx, nb_idx) == 1.0:
                    indices.append(nb_idx)
                    for o2_nb in graph.adjacency[nb_idx]:
                        if o2_nb != idx and get_atom(graph, o2_nb).symbol == "O":
                            indices.append(o2_nb)
                    break
            groups.append(FunctionalGroup(
                group_type="peroxyacid",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["peroxyacid"],
            ))
            continue

        # エステル: C(=O)-O-R (priority=90)
        if _is_ester(graph, idx):
            o_idx = _get_double_bonded_oxygen(graph, idx)
            indices = [idx] + ([o_idx] if o_idx is not None else [])
            # エステル酸素 (単結合 O-R) も追加
            for nb_idx in graph.adjacency[idx]:
                nb = get_atom(graph, nb_idx)
                if nb.symbol == "O" and get_bond_order(graph, idx, nb_idx) == 1.0:
                    indices.append(nb_idx)
                    break
            groups.append(FunctionalGroup(
                group_type="ester",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["ester"],
            ))
            continue

        # 酸ハライド: C(=O)-X (priority=85)
        if _is_acid_halide(graph, idx):
            o_idx = _get_double_bonded_oxygen(graph, idx)
            indices = [idx] + ([o_idx] if o_idx is not None else [])
            groups.append(FunctionalGroup(
                group_type="acid_halide",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["acid_halide"],
            ))
            continue

        # カルバミン酸: RnN-C(=O)-OH (保留名 carbamic acid, priority=99, Phase 71)
        if _is_carbamic_acid(graph, idx):
            o_idx = _get_double_bonded_oxygen(graph, idx)
            # N を直接取得 (H なし N,N-二置換も対応)
            n_idx_ca = next(
                (nb for nb in graph.adjacency[idx]
                 if get_atom(graph, nb).symbol == "N"
                 and get_bond_order(graph, idx, nb) == 1.0),
                None,
            )
            oh_idx = _get_carbamic_oh(graph, idx)
            indices = [idx]
            if o_idx is not None:
                indices.append(o_idx)
            if n_idx_ca is not None:
                indices.append(n_idx_ca)
            if oh_idx is not None:
                indices.append(oh_idx)
            groups.append(FunctionalGroup(
                group_type="carbamic_acid",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["carbamic_acid"],
            ))
            continue

        # チオアミド: C(=S)-NR₂ (priority=93, Phase 41) — アミン検出より先に処理
        if _is_thioamide(graph, idx):
            s_idx = _get_double_bonded_sulfur(graph, idx)
            n_idx_ta = _get_thioamide_nitrogen(graph, idx)
            indices = [idx]
            if s_idx is not None:
                indices.append(s_idx)
            if n_idx_ta is not None:
                indices.append(n_idx_ta)
            groups.append(FunctionalGroup(
                group_type="thioamide",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["thioamide"],
            ))
            continue

        # チオカルボン酸: Phase 149 — チオアルデヒドより先に判定
        if not atom.in_ring:
            _o2_tc = [nb for nb in graph.adjacency[idx]
                      if get_atom(graph, nb).symbol == "O"
                      and get_bond_order(graph, idx, nb) == 2.0]
            _s2_tc = [nb for nb in graph.adjacency[idx]
                      if get_atom(graph, nb).symbol == "S"
                      and get_bond_order(graph, idx, nb) == 2.0]
            _s1_tc = [nb for nb in graph.adjacency[idx]
                      if get_atom(graph, nb).symbol == "S"
                      and get_bond_order(graph, idx, nb) == 1.0
                      and (get_atom(graph, nb).num_hs >= 1
                           or any(get_atom(graph, hh).symbol == "H"
                                  for hh in graph.adjacency[nb]))]
            _o1_tc = [nb for nb in graph.adjacency[idx]
                      if get_atom(graph, nb).symbol == "O"
                      and get_bond_order(graph, idx, nb) == 1.0
                      and (get_atom(graph, nb).num_hs >= 1
                           or any(get_atom(graph, hh).symbol == "H"
                                  for hh in graph.adjacency[nb]))]
            if _o2_tc and _s1_tc:
                # C(=O)-SH: thioic S-acid
                groups.append(FunctionalGroup(
                    group_type="thioic_s_acid",
                    atom_indices=[idx] + _o2_tc + _s1_tc,
                    priority=FUNCTIONAL_GROUP_PRIORITY["thioic_s_acid"],
                ))
                continue
            if _s2_tc and _o1_tc:
                # C(=S)-OH: thioic O-acid
                groups.append(FunctionalGroup(
                    group_type="thioic_o_acid",
                    atom_indices=[idx] + _s2_tc + _o1_tc,
                    priority=FUNCTIONAL_GROUP_PRIORITY["thioic_o_acid"],
                ))
                continue
            if _s2_tc and _s1_tc:
                # C(=S)-SH: dithioic acid
                groups.append(FunctionalGroup(
                    group_type="dithioic_acid",
                    atom_indices=[idx] + _s2_tc + _s1_tc,
                    priority=FUNCTIONAL_GROUP_PRIORITY["dithioic_acid"],
                ))
                continue

        # チオアルデヒド / チオケトン: C=S (Phase 121)
        # N=C=S (イソチオシアネート) や C=S-N (チオアミド) は除外
        if not atom.in_ring:
            s_idx = _get_double_bonded_sulfur(graph, idx)
            if s_idx is not None and not get_atom(graph, s_idx).in_ring:
                # N または O への二重結合がないことを確認 (isothiocyanate 等を除外)
                has_other_double = any(
                    get_bond_order(graph, idx, nb) == 2.0 and nb != s_idx
                    for nb in graph.adjacency[idx]
                    if get_atom(graph, nb).symbol in ("N", "O")
                )
                if not has_other_double:
                    c_single = [nb for nb in graph.adjacency[idx]
                                if nb != s_idx and get_atom(graph, nb).symbol == "C"
                                and get_bond_order(graph, idx, nb) == 1.0]
                    if len(c_single) <= 1:
                        groups.append(FunctionalGroup(
                            group_type="thioaldehyde",
                            atom_indices=[idx, s_idx],
                            priority=FUNCTIONAL_GROUP_PRIORITY["thioaldehyde"],
                        ))
                        continue
                    elif len(c_single) == 2:
                        groups.append(FunctionalGroup(
                            group_type="thioketone",
                            atom_indices=[idx, s_idx],
                            priority=FUNCTIONAL_GROUP_PRIORITY["thioketone"],
                        ))
                        continue

        # ヒドラジド: C(=O)-NH-NH₂ (priority=94, Phase 75) — アミドより先に判定
        if _is_hydrazide(graph, idx):
            o_idx = _get_double_bonded_oxygen(graph, idx)
            n_idx = _get_amide_nitrogen(graph, idx)
            indices = [idx]
            if o_idx is not None:
                indices.append(o_idx)
            if n_idx is not None:
                indices.append(n_idx)
            groups.append(FunctionalGroup(
                group_type="hydrazide",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["hydrazide"],
            ))
            continue

        # アミド: C(=O)-NH₂ (第一級、priority=95)
        if _is_amide(graph, idx):
            o_idx = _get_double_bonded_oxygen(graph, idx)
            n_idx = _get_amide_nitrogen(graph, idx)
            indices = [idx]
            if o_idx is not None:
                indices.append(o_idx)
            if n_idx is not None:
                indices.append(n_idx)
            groups.append(FunctionalGroup(
                group_type="amide",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["amide"],
            ))
            continue

        # カルボキシレートアニオン: C(=O)[O-]  (Phase 146)
        if _is_carboxylate(graph, idx):
            o_indices = _get_carboxylate_oxygens(graph, idx)
            groups.append(FunctionalGroup(
                group_type="carboxylate",
                atom_indices=[idx] + o_indices,
                priority=FUNCTIONAL_GROUP_PRIORITY.get("carboxylate", 97),
            ))
            continue

        # カルボン酸: C(=O)O-H
        # 条件: C に O が二重結合、かつ別の O-H（or O-） が単結合
        if _is_carboxylic_acid(graph, idx):
            o_indices = _get_carbonyl_oxygens(graph, idx)
            groups.append(FunctionalGroup(
                group_type="carboxylic_acid",
                atom_indices=[idx] + o_indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["carboxylic_acid"],
            ))
            continue  # 同一炭素に他の官能基は付かない

        # シアン酸エステル: O-C≡N (Phase 69) — ニトリルより先に判定
        if _is_cyanate(graph, idx):
            n_idx_cy = _get_triple_bonded_nitrogen(graph, idx)
            o_idx_cy = _get_cyanate_oxygen(graph, idx)
            indices = [idx]
            if n_idx_cy is not None:
                indices.append(n_idx_cy)
            if o_idx_cy is not None:
                indices.append(o_idx_cy)
            groups.append(FunctionalGroup(
                group_type="cyanate",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["cyanate"],
            ))
            continue

        # チオシアン酸エステル: S-C≡N (Phase 69) — ニトリルより先に判定
        if _is_thiocyanate(graph, idx):
            n_idx_tc = _get_triple_bonded_nitrogen(graph, idx)
            s_idx_tc = _get_thiocyanate_sulfur(graph, idx)
            indices = [idx]
            if n_idx_tc is not None:
                indices.append(n_idx_tc)
            if s_idx_tc is not None:
                indices.append(s_idx_tc)
            groups.append(FunctionalGroup(
                group_type="thiocyanate",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["thiocyanate"],
            ))
            continue

        # ニトリル: C≡N (末端)
        if _is_nitrile(graph, idx):
            n_idx = _get_triple_bonded_nitrogen(graph, idx)
            indices = [idx] + ([n_idx] if n_idx is not None else [])
            groups.append(FunctionalGroup(
                group_type="nitrile",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["nitrile"],
            ))
            continue

        # アルデヒド: C(=O)H (末端カルボニル、H付き)
        if _is_aldehyde(graph, idx):
            o_idx = _get_double_bonded_oxygen(graph, idx)
            indices = [idx] + ([o_idx] if o_idx is not None else [])
            groups.append(FunctionalGroup(
                group_type="aldehyde",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["aldehyde"],
            ))
            continue

        # ケトン: C(=O)C (カルボニル、両側が炭素)
        if _is_ketone(graph, idx):
            o_idx = _get_double_bonded_oxygen(graph, idx)
            indices = [idx] + ([o_idx] if o_idx is not None else [])
            groups.append(FunctionalGroup(
                group_type="ketone",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["ketone"],
            ))
            continue

        # カルボジイミド: R-N=C=N-R (Phase 73)
        if _is_carbodiimide(graph, idx):
            n_double = [nb for nb in graph.adjacency[idx]
                        if get_atom(graph, nb).symbol == "N"
                        and get_bond_order(graph, idx, nb) == 2.0]
            groups.append(FunctionalGroup(
                group_type="carbodiimide",
                atom_indices=[idx] + n_double,
                priority=FUNCTIONAL_GROUP_PRIORITY["carbodiimide"],
            ))
            continue

        # アミジン: C(=N-H)(N-H) (Phase 70) — imine より先に検出
        if _is_amidine(graph, idx):
            n_imine_idx, n_amine_idx = _get_amidine_nitrogens(graph, idx)
            indices = [idx]
            if n_imine_idx is not None:
                indices.append(n_imine_idx)
            if n_amine_idx is not None:
                indices.append(n_amine_idx)
            groups.append(FunctionalGroup(
                group_type="amidine",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["amidine"],
            ))
            continue

        # イミン: C=N-H (第一級イミン)
        if _is_imine(graph, idx):
            n_idx_imine = _get_imine_nitrogen(graph, idx)
            indices = [idx] + ([n_idx_imine] if n_idx_imine is not None else [])
            groups.append(FunctionalGroup(
                group_type="imine",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["imine"],
            ))
            continue

        # オキシム: C=N-OH (ケトオキシムまたはアルドキシム)
        if _is_ketoxime(graph, idx):
            n_idx_ox = _get_oxime_nitrogen(graph, idx)
            indices = [idx] + ([n_idx_ox] if n_idx_ox is not None else [])
            groups.append(FunctionalGroup(
                group_type="ketoxime",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["ketoxime"],
            ))
            continue

        if _is_aldoxime(graph, idx):
            n_idx_ox = _get_oxime_nitrogen(graph, idx)
            indices = [idx] + ([n_idx_ox] if n_idx_ox is not None else [])
            groups.append(FunctionalGroup(
                group_type="aldoxime",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["aldoxime"],
            ))
            continue

        # セミカルバゾン: C=N-N-C(=O)-NH₂ (Phase 80) — hydrazone より先に判定
        _sc = _is_semicarbazone_or_thio(graph, idx)
        if _sc is not None:
            n_idx_hz = _get_hydrazone_nitrogen(graph, idx)
            indices = [idx] + ([n_idx_hz] if n_idx_hz is not None else [])
            groups.append(FunctionalGroup(
                group_type=_sc,
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY[_sc],
            ))
            continue

        # ヒドラゾン: C=N-NH₂ / C=N-NHR (Phase 43)
        if _is_kethydrazone(graph, idx):
            n_idx_hz = _get_hydrazone_nitrogen(graph, idx)
            indices = [idx] + ([n_idx_hz] if n_idx_hz is not None else [])
            groups.append(FunctionalGroup(
                group_type="kethydrazone",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["kethydrazone"],
            ))
            continue

        if _is_aldhydrazone(graph, idx):
            n_idx_hz = _get_hydrazone_nitrogen(graph, idx)
            indices = [idx] + ([n_idx_hz] if n_idx_hz is not None else [])
            groups.append(FunctionalGroup(
                group_type="aldhydrazone",
                atom_indices=indices,
                priority=FUNCTIONAL_GROUP_PRIORITY["aldhydrazone"],
            ))
            continue

    # アルコール: C-O-H (C に隣接する OH)
    # ※ カルボン酸・アルデヒド・ケトンでない炭素に付く OH
    carbonyl_carbons = {g.atom_indices[0] for g in groups}
    for atom in graph.atoms:
        if atom.symbol != "O":
            continue
        o_idx = atom.idx
        neighbors = graph.adjacency[o_idx]
        c_neighbors = [n for n in neighbors if get_atom(graph, n).symbol == "C"]
        h_neighbors = [n for n in neighbors if get_atom(graph, n).symbol == "H"]
        if c_neighbors and h_neighbors:
            c_idx = c_neighbors[0]
            if c_idx not in carbonyl_carbons:
                groups.append(FunctionalGroup(
                    group_type="alcohol",
                    atom_indices=[c_idx, o_idx],
                    priority=FUNCTIONAL_GROUP_PRIORITY["alcohol"],
                ))

    # ヒドロペルオキシド: C-O1-O2-H (Phase 44)
    # Phase 45: C が carbonyl C (C=O) の場合はペルオキシ酸なので除外
    seen_peroxy: set[tuple[int, int]] = set()
    for atom in graph.atoms:
        if atom.symbol != "O":
            continue
        o1_idx = atom.idx
        # O1 に C 隣接 AND O 隣接があるか
        o1_c_nbrs = [nb for nb in graph.adjacency[o1_idx] if get_atom(graph, nb).symbol == "C"]
        o1_o_nbrs = [nb for nb in graph.adjacency[o1_idx] if get_atom(graph, nb).symbol == "O"]
        if not o1_c_nbrs or not o1_o_nbrs:
            continue
        o2_idx = o1_o_nbrs[0]
        # O2 に H があること
        if not any(get_atom(graph, nb).symbol == "H" for nb in graph.adjacency[o2_idx]):
            continue
        # O1 に隣接する C が carbonyl C (=O あり) の場合はペルオキシ酸 → 除外
        c_idx = o1_c_nbrs[0]
        if _get_double_bonded_oxygen(graph, c_idx) is not None:
            continue
        pair = (min(o1_idx, o2_idx), max(o1_idx, o2_idx))
        if pair in seen_peroxy:
            continue
        seen_peroxy.add(pair)
        groups.append(FunctionalGroup(
            group_type="hydroperoxide",
            atom_indices=[c_idx, o1_idx, o2_idx],
            priority=FUNCTIONAL_GROUP_PRIORITY["hydroperoxide"],
        ))

    # 有機ペルオキシド: C-O-O-C (H なし) (Phase 57)
    seen_peroxide: set[tuple[int, int]] = set()
    for atom in graph.atoms:
        if atom.symbol != "O":
            continue
        o1_idx = atom.idx
        o1_c_nbrs = [nb for nb in graph.adjacency[o1_idx] if get_atom(graph, nb).symbol == "C"]
        o1_o_nbrs = [nb for nb in graph.adjacency[o1_idx] if get_atom(graph, nb).symbol == "O"]
        if not o1_c_nbrs or not o1_o_nbrs:
            continue
        o2_idx = o1_o_nbrs[0]
        # O2 に C 隣接 AND H なし (hydroperoxide は O2 に H あり → 除外済み)
        o2_c_nbrs = [nb for nb in graph.adjacency[o2_idx]
                     if nb != o1_idx and get_atom(graph, nb).symbol == "C"]
        if not o2_c_nbrs:
            continue
        # H なし確認 (両 O に H がないこと)
        if (any(get_atom(graph, nb).symbol == "H" for nb in graph.adjacency[o1_idx])
                or any(get_atom(graph, nb).symbol == "H" for nb in graph.adjacency[o2_idx])):
            continue
        # carbonyl C 除外 (ペルオキシ酸)
        c1_idx = o1_c_nbrs[0]
        c2_idx = o2_c_nbrs[0]
        if (_get_double_bonded_oxygen(graph, c1_idx) is not None
                or _get_double_bonded_oxygen(graph, c2_idx) is not None):
            continue
        pair = (min(o1_idx, o2_idx), max(o1_idx, o2_idx))
        if pair in seen_peroxide:
            continue
        seen_peroxide.add(pair)
        groups.append(FunctionalGroup(
            group_type="peroxide",
            atom_indices=[c1_idx, o1_idx, o2_idx, c2_idx],
            priority=FUNCTIONAL_GROUP_PRIORITY["peroxide"],
        ))

    # チオール / スルホキシド / スルホン / スルホンアミド: S 原子を走査
    for atom in graph.atoms:
        if atom.symbol != "S":
            continue
        if atom.in_ring:
            continue
        s_idx = atom.idx
        neighbors = graph.adjacency[s_idx]
        h_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "H"]
        c_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "C"]
        n_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "N"]
        o_double = [nb for nb in neighbors
                    if get_atom(graph, nb).symbol == "O" and get_bond_order(graph, s_idx, nb) == 2.0]

        # スルホン酸: C-S(=O)₂-OH
        o_single = [nb for nb in neighbors
                    if get_atom(graph, nb).symbol == "O"
                    and get_bond_order(graph, s_idx, nb) == 1.0]
        o_single_oh = [nb for nb in o_single
                       if any(get_atom(graph, onh).symbol == "H"
                              for onh in graph.adjacency[nb])]

        # スルホニルクロライド: C-S(=O)₂-Cl (Phase 59)
        cl_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "Cl"]
        if (len(o_double) == 2 and len(c_neighbors) == 1 and len(cl_neighbors) == 1
                and len(n_neighbors) == 0 and not o_single_oh):
            groups.append(FunctionalGroup(
                group_type="sulfonyl_chloride",
                atom_indices=[s_idx] + o_double + c_neighbors + cl_neighbors,
                priority=FUNCTIONAL_GROUP_PRIORITY["sulfonyl_chloride"],
            ))

        elif (len(o_double) == 2 and len(c_neighbors) == 1 and len(o_single_oh) >= 1
                and len(n_neighbors) == 0):
            groups.append(FunctionalGroup(
                group_type="sulfonic_acid",
                atom_indices=[s_idx] + o_double + c_neighbors + o_single_oh,
                priority=FUNCTIONAL_GROUP_PRIORITY["sulfonic_acid"],
            ))
        elif (len(o_double) == 1 and len(c_neighbors) == 1 and len(o_single_oh) >= 1
              and len(n_neighbors) == 0):
            # スルフィン酸: C-S(=O)-OH (Phase 38)
            groups.append(FunctionalGroup(
                group_type="sulfinic_acid",
                atom_indices=[s_idx] + o_double + c_neighbors + o_single_oh,
                priority=FUNCTIONAL_GROUP_PRIORITY["sulfinic_acid"],
            ))
        elif len(h_neighbors) >= 1 and len(c_neighbors) == 1 and not o_double:
            # チオール: C-SH
            c_idx = c_neighbors[0]
            groups.append(FunctionalGroup(
                group_type="thiol",
                atom_indices=[c_idx, s_idx],
                priority=FUNCTIONAL_GROUP_PRIORITY["thiol"],
            ))
        elif (len(o_double) == 2 and len(n_neighbors) == 2 and len(c_neighbors) == 0):
            # スルファミド: H2N-S(=O)₂-NH2 (Phase 67)
            groups.append(FunctionalGroup(
                group_type="sulfamide",
                atom_indices=[s_idx] + o_double + n_neighbors,
                priority=FUNCTIONAL_GROUP_PRIORITY["sulfamide"],
            ))
        elif (len(o_double) == 2 and len(n_neighbors) == 1 and len(c_neighbors) == 0
              and len(o_single_oh) >= 1):
            # スルファミン酸: H2N-S(=O)₂-OH (Phase 67)
            groups.append(FunctionalGroup(
                group_type="sulfamic_acid",
                atom_indices=[s_idx] + o_double + n_neighbors + o_single_oh,
                priority=FUNCTIONAL_GROUP_PRIORITY["sulfamic_acid"],
            ))
        elif len(o_double) == 2 and len(n_neighbors) >= 1 and len(c_neighbors) == 1:
            # スルホンアミド: C-S(=O)₂-N
            groups.append(FunctionalGroup(
                group_type="sulfonamide",
                atom_indices=[s_idx] + o_double + c_neighbors + n_neighbors,
                priority=FUNCTIONAL_GROUP_PRIORITY["sulfonamide"],
            ))
        elif (len(o_double) == 2 and len(c_neighbors) == 1
              and not o_single_oh and len(n_neighbors) == 0):
            # スルホン酸エステル: C-S(=O)₂-O-C (Phase 82) — sulfone より先に
            o_ester = [nb for nb in o_single if nb not in o_single_oh
                       and any(get_atom(graph, occ).symbol == "C"
                               for occ in graph.adjacency[nb] if occ != s_idx)]
            if o_ester:
                groups.append(FunctionalGroup(
                    group_type="sulfonate_ester",
                    atom_indices=[s_idx] + o_double + c_neighbors + [o_ester[0]],
                    priority=FUNCTIONAL_GROUP_PRIORITY["sulfonate_ester"],
                ))
            else:
                # 通常スルホン (両側 C)
                pass  # len(c_neighbors)==1 では sulfone にならない
        elif (len(o_double) == 1 and len(c_neighbors) == 1
              and not o_single_oh and len(n_neighbors) == 0):
            # スルフィン酸エステル: C-S(=O)-O-C (Phase 82) — sulfoxide より先に
            o_ester = [nb for nb in o_single if nb not in o_single_oh
                       and any(get_atom(graph, occ).symbol == "C"
                               for occ in graph.adjacency[nb] if occ != s_idx)]
            if o_ester:
                groups.append(FunctionalGroup(
                    group_type="sulfinate_ester",
                    atom_indices=[s_idx] + o_double + c_neighbors + [o_ester[0]],
                    priority=FUNCTIONAL_GROUP_PRIORITY["sulfinate_ester"],
                ))
            else:
                pass  # will not match sulfoxide (len(c_neighbors)==1)
        elif len(o_double) == 2 and len(c_neighbors) == 2:
            # スルホン: C-S(=O)₂-C
            groups.append(FunctionalGroup(
                group_type="sulfone",
                atom_indices=[s_idx] + o_double + c_neighbors,
                priority=FUNCTIONAL_GROUP_PRIORITY["sulfone"],
            ))
        elif len(o_double) == 1 and len(c_neighbors) == 2:
            # スルホキシド: C-S(=O)-C
            groups.append(FunctionalGroup(
                group_type="sulfoxide",
                atom_indices=[s_idx] + o_double + c_neighbors,
                priority=FUNCTIONAL_GROUP_PRIORITY["sulfoxide"],
            ))
        elif (len(h_neighbors) == 0 and len(c_neighbors) == 2
              and not o_double and len(n_neighbors) == 0):
            # チオエステル vs チオエーテル
            carbonyl_cs = [c for c in c_neighbors if _has_double_bonded_oxygen(graph, c)]
            if carbonyl_cs:
                # チオエステル: R-C(=O)-S-R' (Phase 55)
                alkyl_cs = [c for c in c_neighbors if c not in carbonyl_cs]
                groups.append(FunctionalGroup(
                    group_type="thioester",
                    atom_indices=[s_idx, carbonyl_cs[0]] + alkyl_cs,
                    priority=FUNCTIONAL_GROUP_PRIORITY["thioester"],
                ))
            else:
                # チオエーテル: C-S-C
                groups.append(FunctionalGroup(
                    group_type="sulfide",
                    atom_indices=[s_idx] + c_neighbors,
                    priority=FUNCTIONAL_GROUP_PRIORITY["sulfide"],
                ))
        elif (len(h_neighbors) == 0 and len(c_neighbors) == 1
              and not o_double and len(n_neighbors) == 0):
            # ジスルフィド: C-S-S-C (Phase 56)
            s_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "S"]
            if len(s_neighbors) == 1:
                s2_idx = s_neighbors[0]
                s2_c = [nb for nb in graph.adjacency[s2_idx]
                        if nb != s_idx and get_atom(graph, nb).symbol == "C"]
                if s2_c and s_idx < s2_idx:  # 重複登録を防ぐ
                    groups.append(FunctionalGroup(
                        group_type="disulfide",
                        atom_indices=[s_idx, s2_idx] + c_neighbors + s2_c,
                        priority=FUNCTIONAL_GROUP_PRIORITY["disulfide"],
                    ))

    # アミン検出: 環外 N を走査（第一級・第二級・第三級）
    # アミド (C=O に結合した N) は除外する
    for atom in graph.atoms:
        if atom.symbol != "N":
            continue
        if atom.in_ring:
            continue
        if atom.formal_charge == 1:
            # N+ で O 隣接なし (ammonium) は別途処理; N-oxide/nitro はそのまま
            n_idx_chk = atom.idx
            o_nb_chk = [nb for nb in graph.adjacency[n_idx_chk]
                        if get_atom(graph, nb).symbol == "O"]
            n_nb_chk = [nb for nb in graph.adjacency[n_idx_chk]
                        if get_atom(graph, nb).symbol == "N"]
            if not o_nb_chk and not n_nb_chk:
                continue  # 純アンモニウムは ammonium 検出へ
        n_idx = atom.idx
        neighbors = graph.adjacency[n_idx]
        h_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "H"]
        c_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "C"]

        # イソシアネート / イソチオシアネート: R-N=C=O/S (Phase 54/68)
        c_dbl = [nb for nb in c_neighbors if get_bond_order(graph, n_idx, nb) == 2.0]
        c_sgl = [nb for nb in c_neighbors if get_bond_order(graph, n_idx, nb) == 1.0]
        if len(c_dbl) == 1 and len(c_sgl) == 1 and not h_neighbors:
            iso_c = c_dbl[0]
            for iso_c_nb in graph.adjacency[iso_c]:
                iso_c_nb_sym = get_atom(graph, iso_c_nb).symbol
                if get_bond_order(graph, iso_c, iso_c_nb) == 2.0:
                    if iso_c_nb_sym == "O":
                        groups.append(FunctionalGroup(
                            group_type="isocyanate",
                            atom_indices=[c_sgl[0], n_idx, iso_c],
                            priority=FUNCTIONAL_GROUP_PRIORITY["isocyanate"],
                        ))
                        break
                    elif iso_c_nb_sym == "S":
                        groups.append(FunctionalGroup(
                            group_type="isothiocyanate",
                            atom_indices=[c_sgl[0], n_idx, iso_c, iso_c_nb],
                            priority=FUNCTIONAL_GROUP_PRIORITY["isothiocyanate"],
                        ))
                        break
            continue

        # アジド: R-N=[N+]=[N-] (N に単結合 C + 二重結合 N) (Phase 53)
        n_dbl = [nb for nb in graph.adjacency[n_idx]
                 if get_atom(graph, nb).symbol == "N"
                 and get_bond_order(graph, n_idx, nb) == 2.0]
        if len(n_dbl) == 1 and len(c_neighbors) == 1 and not h_neighbors:
            n2_idx = n_dbl[0]
            # N2 の隣に N3 があるか (N2=N3 または N2-N3)
            n3_candidates = [nb for nb in graph.adjacency[n2_idx]
                             if nb != n_idx and get_atom(graph, nb).symbol == "N"]
            if n3_candidates:
                groups.append(FunctionalGroup(
                    group_type="azide",
                    atom_indices=[c_neighbors[0], n_idx, n2_idx, n3_candidates[0]],
                    priority=FUNCTIONAL_GROUP_PRIORITY["azide"],
                ))
                continue

        # ニトロソ: R-N=O (N に単結合 C + 二重結合 O のみ, H なし) (Phase 52)
        o_all = [nb for nb in graph.adjacency[n_idx]
                 if get_atom(graph, nb).symbol == "O"]
        o_dbl = [nb for nb in o_all
                 if get_bond_order(graph, n_idx, nb) == 2.0]
        if len(o_all) == 1 and len(o_dbl) == 1 and len(c_neighbors) == 1 and not h_neighbors:
            groups.append(FunctionalGroup(
                group_type="nitroso",
                atom_indices=[c_neighbors[0], n_idx, o_dbl[0]],
                priority=FUNCTIONAL_GROUP_PRIORITY["nitroso"],
            ))
            continue

        if len(h_neighbors) >= 2 and len(c_neighbors) == 1:
            # 第一級アミン: NH₂ に C が 1 個
            c_idx = c_neighbors[0]
            if _has_double_bonded_oxygen(graph, c_idx):
                continue
            if _has_double_bonded_sulfur(graph, c_idx):  # チオアミド除外
                continue
            groups.append(FunctionalGroup(
                group_type="amine",
                atom_indices=[c_idx, n_idx],
                priority=FUNCTIONAL_GROUP_PRIORITY["amine"],
            ))
        elif len(h_neighbors) == 1 and len(c_neighbors) == 2:
            # 第二級アミン: NH に C が 2 個
            if any(_has_double_bonded_oxygen(graph, c) for c in c_neighbors):
                continue  # 二級アミドは除外
            if any(_has_double_bonded_sulfur(graph, c) for c in c_neighbors):
                continue  # 二級チオアミドは除外
            groups.append(FunctionalGroup(
                group_type="amine",
                atom_indices=[n_idx] + c_neighbors,
                priority=FUNCTIONAL_GROUP_PRIORITY["amine"],
            ))
        elif len(h_neighbors) == 0 and len(c_neighbors) == 3:
            # 第三級アミン: N に C が 3 個
            if any(_has_double_bonded_oxygen(graph, c) for c in c_neighbors):
                continue  # アミド様結合は除外
            if any(_has_double_bonded_sulfur(graph, c) for c in c_neighbors):
                continue  # 三級チオアミドは除外
            groups.append(FunctionalGroup(
                group_type="amine",
                atom_indices=[n_idx] + c_neighbors,
                priority=FUNCTIONAL_GROUP_PRIORITY["amine"],
            ))

    # アルケン: C=C (非芳香族、非環状のみ)
    seen_double = set()
    for atom in graph.atoms:
        if atom.symbol != "C" or atom.in_ring:
            continue
        for neighbor_idx in graph.adjacency[atom.idx]:
            neighbor = get_atom(graph, neighbor_idx)
            if neighbor.symbol != "C" or neighbor.in_ring:
                continue
            bo = get_bond_order(graph, atom.idx, neighbor_idx)
            if bo == 2.0:
                pair = (min(atom.idx, neighbor_idx), max(atom.idx, neighbor_idx))
                if pair not in seen_double:
                    seen_double.add(pair)
                    groups.append(FunctionalGroup(
                        group_type="alkene",
                        atom_indices=list(pair),
                        priority=FUNCTIONAL_GROUP_PRIORITY["alkene"],
                    ))

    # アルキン: C≡C (非環状のみ)
    seen_triple = set()
    for atom in graph.atoms:
        if atom.symbol != "C" or atom.in_ring:
            continue
        for neighbor_idx in graph.adjacency[atom.idx]:
            neighbor = get_atom(graph, neighbor_idx)
            if neighbor.symbol != "C" or neighbor.in_ring:
                continue
            bo = get_bond_order(graph, atom.idx, neighbor_idx)
            if bo == 3.0:
                pair = (min(atom.idx, neighbor_idx), max(atom.idx, neighbor_idx))
                if pair not in seen_triple:
                    seen_triple.add(pair)
                    groups.append(FunctionalGroup(
                        group_type="alkyne",
                        atom_indices=list(pair),
                        priority=FUNCTIONAL_GROUP_PRIORITY["alkyne"],
                    ))

    # ─── Phase 143: リン化合物検出 ────────────────────────────────────────
    for atom in graph.atoms:
        if atom.symbol != "P" or atom.in_ring:
            continue
        p_idx = atom.idx
        neighbors = graph.adjacency[p_idx]
        c_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "C"]
        o_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "O"]
        o_double = [nb for nb in o_neighbors if get_bond_order(graph, p_idx, nb) == 2.0]
        o_single = [nb for nb in o_neighbors if get_bond_order(graph, p_idx, nb) == 1.0]
        o_single_oh = [nb for nb in o_single
                       if any(get_atom(graph, onh).symbol == "H"
                              for onh in graph.adjacency[nb])]

        if len(o_double) == 1 and len(o_single_oh) >= 2 and len(c_neighbors) == 1:
            # ホスホン酸: R-P(=O)(OH)2
            groups.append(FunctionalGroup(
                group_type="phosphonic_acid",
                atom_indices=[p_idx] + c_neighbors + o_double + o_single_oh,
                priority=FUNCTIONAL_GROUP_PRIORITY["phosphonic_acid"],
            ))
        elif len(o_double) == 1 and len(o_single_oh) >= 1 and len(c_neighbors) == 2:
            # ホスフィン酸: R2P(=O)(OH)
            groups.append(FunctionalGroup(
                group_type="phosphinic_acid",
                atom_indices=[p_idx] + c_neighbors + o_double + o_single_oh,
                priority=FUNCTIONAL_GROUP_PRIORITY["phosphinic_acid"],
            ))
        elif c_neighbors and not o_neighbors:
            # ホスファン: R_n-PH_{3-n}
            groups.append(FunctionalGroup(
                group_type="phosphane",
                atom_indices=[p_idx] + c_neighbors,
                priority=FUNCTIONAL_GROUP_PRIORITY["phosphane"],
            ))
        elif len(o_double) == 1 and not o_single_oh:
            # O-アルキル結合のある P=O 化合物 (エステル類)
            o_ester = [nb for nb in o_single
                       if any(get_atom(graph, occ).symbol == "C"
                              for occ in graph.adjacency[nb] if occ != p_idx)]
            if o_ester and len(c_neighbors) == 0:
                # ホスフェートエステル: (RO)3P=O
                groups.append(FunctionalGroup(
                    group_type="phosphate_ester",
                    atom_indices=[p_idx] + o_double + o_ester,
                    priority=FUNCTIONAL_GROUP_PRIORITY.get("phosphate_ester", 87),
                ))
            elif o_ester and len(c_neighbors) == 1:
                # ホスホネートエステル: R-P(=O)(OR)2
                groups.append(FunctionalGroup(
                    group_type="phosphonate_ester",
                    atom_indices=[p_idx] + c_neighbors + o_double + o_ester,
                    priority=FUNCTIONAL_GROUP_PRIORITY.get("phosphonate_ester", 86),
                ))
            elif o_ester and len(c_neighbors) == 2:
                # ホスフィネートエステル: R2-P(=O)(OR)
                groups.append(FunctionalGroup(
                    group_type="phosphinate_ester",
                    atom_indices=[p_idx] + c_neighbors + o_double + o_ester,
                    priority=FUNCTIONAL_GROUP_PRIORITY.get("phosphinate_ester", 85),
                ))

    # ─── Phase 143: ホウ素化合物検出 ──────────────────────────────────────
    for atom in graph.atoms:
        if atom.symbol != "B" or atom.in_ring:
            continue
        b_idx = atom.idx
        neighbors = graph.adjacency[b_idx]
        c_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "C"]
        o_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "O"]
        o_single_oh = [nb for nb in o_neighbors
                       if get_bond_order(graph, b_idx, nb) == 1.0
                       and any(get_atom(graph, onh).symbol == "H"
                               for onh in graph.adjacency[nb])]

        if len(o_single_oh) >= 2 and len(c_neighbors) == 1:
            # ボロン酸: R-B(OH)2
            groups.append(FunctionalGroup(
                group_type="boronic_acid",
                atom_indices=[b_idx] + c_neighbors + o_single_oh,
                priority=FUNCTIONAL_GROUP_PRIORITY["boronic_acid"],
            ))
        elif len(o_single_oh) == 1 and len(c_neighbors) == 2:
            # ボリン酸: R2B(OH)
            groups.append(FunctionalGroup(
                group_type="borinic_acid",
                atom_indices=[b_idx] + c_neighbors + o_single_oh,
                priority=FUNCTIONAL_GROUP_PRIORITY["borinic_acid"],
            ))
        elif c_neighbors and not o_neighbors:
            # ボラン: R_n-BH_{3-n}
            groups.append(FunctionalGroup(
                group_type="borane_org",
                atom_indices=[b_idx] + c_neighbors,
                priority=FUNCTIONAL_GROUP_PRIORITY["borane_org"],
            ))

    # ─── Phase 143: ケイ素化合物検出 ──────────────────────────────────────
    for atom in graph.atoms:
        if atom.symbol != "Si" or atom.in_ring:
            continue
        si_idx = atom.idx
        neighbors = graph.adjacency[si_idx]
        c_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "C"]
        if c_neighbors:
            groups.append(FunctionalGroup(
                group_type="silane_org",
                atom_indices=[si_idx] + c_neighbors,
                priority=FUNCTIONAL_GROUP_PRIORITY["silane_org"],
            ))

    # ─── Phase 146: アンモニウムイオン検出 ───────────────────────────────
    for atom in graph.atoms:
        if atom.symbol != "N" or atom.in_ring or atom.formal_charge != 1:
            continue
        n_idx = atom.idx
        neighbors = graph.adjacency[n_idx]
        # O 隣接あり → N-oxide or nitro group → ammonium ではない
        o_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "O"]
        if o_neighbors:
            continue
        # N=N 隣接あり → diazo/azo → ammonium ではない
        n_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "N"]
        if n_neighbors:
            continue
        c_neighbors = [nb for nb in neighbors if get_atom(graph, nb).symbol == "C"]
        if c_neighbors:
            groups.append(FunctionalGroup(
                group_type="ammonium",
                atom_indices=[n_idx] + c_neighbors,
                priority=FUNCTIONAL_GROUP_PRIORITY.get("ammonium", 11),
            ))

    # ─── Phase 150: ニトレート/ニトライトエステル検出 ─────────────────
    for atom in graph.atoms:
        if atom.symbol != "N" or atom.in_ring:
            continue
        n_idx = atom.idx
        n_nbs = graph.adjacency[n_idx]
        o_double = [nb for nb in n_nbs
                    if get_atom(graph, nb).symbol == "O"
                    and get_bond_order(graph, n_idx, nb) == 2.0]
        o_single = [nb for nb in n_nbs
                    if get_atom(graph, nb).symbol == "O"
                    and get_bond_order(graph, n_idx, nb) == 1.0]
        if atom.formal_charge == 1 and len(o_double) == 1 and len(o_single) == 2:
            # Nitrate ester: R-O-[N+](=O)[O-]
            # One of the single-O must be an alkoxy (has C neighbor), the other is [O-]
            o_alkyl = [ob for ob in o_single
                       if any(get_atom(graph, c).symbol == "C"
                              for c in graph.adjacency[ob])]
            o_neg = [ob for ob in o_single if get_atom(graph, ob).formal_charge == -1]
            if len(o_alkyl) == 1 and len(o_neg) >= 1:
                groups.append(FunctionalGroup(
                    group_type="nitrate_ester",
                    atom_indices=[n_idx, o_alkyl[0]] + o_double + o_neg[:1],
                    priority=FUNCTIONAL_GROUP_PRIORITY["nitrate_ester"],
                ))
        elif atom.formal_charge == 0 and len(o_double) == 1 and len(o_single) == 1:
            # Nitrite ester: R-O-N=O
            o_alkyl_n = [ob for ob in o_single
                         if any(get_atom(graph, c).symbol == "C"
                                for c in graph.adjacency[ob])]
            if len(o_alkyl_n) == 1:
                groups.append(FunctionalGroup(
                    group_type="nitrite_ester",
                    atom_indices=[n_idx, o_alkyl_n[0]] + o_double,
                    priority=FUNCTIONAL_GROUP_PRIORITY["nitrite_ester"],
                ))

    # アルカン: 上記なし
    if not groups:
        groups.append(FunctionalGroup(
            group_type="alkane",
            atom_indices=[],
            priority=FUNCTIONAL_GROUP_PRIORITY["alkane"],
        ))

    # 優先順位の高い順にソート
    groups.sort(key=lambda g: g.priority, reverse=True)

    # 同一最高優先度グループを diol/dione/dioic_acid 等に集約
    groups = aggregate_groups(groups)

    return groups


# ─── 内部ヘルパー ────────────────────────────────────────────

def _is_carboxylic_acid(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C(=O)O-H パターンを持つかチェック。"""
    has_double_o = False
    has_single_oh = False

    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol != "O":
            continue
        bo = get_bond_order(graph, c_idx, nb_idx)
        if bo == 2.0:
            has_double_o = True
        elif bo == 1.0:
            # この O に H が付いているか
            for o_nb in graph.adjacency[nb_idx]:
                if get_atom(graph, o_nb).symbol == "H":
                    has_single_oh = True
                    break

    return has_double_o and has_single_oh


def _is_aldehyde(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C(=O)H パターンを持つかチェック（末端アルデヒド）。"""
    has_double_o = _get_double_bonded_oxygen(graph, c_idx) is not None
    if not has_double_o:
        return False

    # 同じ C に H が直接付いているか（AddHs後）
    has_h = any(
        get_atom(graph, nb).symbol == "H"
        for nb in graph.adjacency[c_idx]
    )

    # カルボン酸の単結合O が存在してはいけない
    has_single_oh = False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "O" and get_bond_order(graph, c_idx, nb_idx) == 1.0:
            has_single_oh = True
            break

    return has_double_o and has_h and not has_single_oh


def _is_ketone(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C(=O)C パターン（ケトン）かチェック。累積 C=C=O (ketene 型) も含む。"""
    o_idx = _get_double_bonded_oxygen(graph, c_idx)
    if o_idx is None:
        return False

    # カルボン酸・アルデヒドでないこと
    if _is_carboxylic_acid(graph, c_idx):
        return False
    if _is_aldehyde(graph, c_idx):
        return False

    c_neighbors = [
        nb for nb in graph.adjacency[c_idx]
        if get_atom(graph, nb).symbol == "C"
    ]
    if len(c_neighbors) >= 2:
        return True

    # Phase 129: ketene 型 C=C=O: 隣接 C との結合が二重結合 (累積系)
    if len(c_neighbors) == 1:
        nb = c_neighbors[0]
        if get_bond_order(graph, c_idx, nb) == 2.0:
            return True

    return False


def _get_double_bonded_oxygen(graph: MoleculeGraph, c_idx: int) -> int | None:
    """C に二重結合で付いている O のインデックスを返す。なければ None。"""
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "O" and get_bond_order(graph, c_idx, nb_idx) == 2.0:
            return nb_idx
    return None


def _get_carbonyl_oxygens(graph: MoleculeGraph, c_idx: int) -> list[int]:
    """カルボン酸 C に付く両方の O インデックスを返す。"""
    result = []
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "O":
            result.append(nb_idx)
    return result


def _is_carboxylate(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C(=O)[O-] パターン (カルボキシレートアニオン) かチェック。

    炭酸水素塩 (bicarbonate: C(=O)([O-])OH) は除外。
    """
    has_double_o = False
    has_neg_o = False
    has_single_oh = False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol != "O":
            continue
        bo = get_bond_order(graph, c_idx, nb_idx)
        if bo == 2.0:
            has_double_o = True
        elif bo == 1.0 and nb.formal_charge == -1:
            has_neg_o = True
        elif bo == 1.0 and nb.num_hs >= 1:
            has_single_oh = True
    # bicarbonate has C=O + C-[O-] + C-OH → exclude by checking no OH
    return has_double_o and has_neg_o and not has_single_oh


def _get_carboxylate_oxygens(graph: MoleculeGraph, c_idx: int) -> list[int]:
    """カルボキシレート C に付く両方の O インデックスを返す。"""
    result = []
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "O":
            result.append(nb_idx)
    return result


def _is_nitrile(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C≡N パターン（末端ニトリル）かチェック。"""
    atom = get_atom(graph, c_idx)
    if atom.in_ring:
        return False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "N" and get_bond_order(graph, c_idx, nb_idx) == 3.0:
            # N が他の重原子と結合していないことを確認
            n_heavy = [
                n for n in graph.adjacency[nb_idx]
                if n != c_idx and get_atom(graph, n).symbol != "H"
            ]
            if not n_heavy:
                return True
    return False


def _get_triple_bonded_nitrogen(graph: MoleculeGraph, c_idx: int) -> int | None:
    """C≡N の N インデックスを返す。なければ None。"""
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "N" and get_bond_order(graph, c_idx, nb_idx) == 3.0:
            return nb_idx
    return None


def _has_double_bonded_oxygen(graph: MoleculeGraph, c_idx: int) -> bool:
    """C に C=O（カルボニル）があるかチェック（アミド判定用）。"""
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "O" and get_bond_order(graph, c_idx, nb_idx) == 2.0:
            return True
    return False


def _is_imine(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C=N-H または C=N-R パターン（一・二級イミン）かチェック。"""
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "N" and get_bond_order(graph, c_idx, nb_idx) == 2.0:
            has_h = any(get_atom(graph, n).symbol == "H" for n in graph.adjacency[nb_idx])
            if has_h:
                return True
            # N-置換イミン: N に C のみが付く (酸素・硫黄なし)
            n_heavy = [n for n in graph.adjacency[nb_idx]
                       if n != c_idx and get_atom(graph, n).symbol != "H"]
            n_c_only = all(get_atom(graph, n).symbol == "C" for n in n_heavy)
            if n_heavy and n_c_only:
                # イソシアネート (C=N=O) / カルボジイミド (C=N=C) 除外:
                # C が N 以外にも二重結合を持つ場合は真のイミンではない
                c_has_other_double = any(
                    nb2 != nb_idx and get_bond_order(graph, c_idx, nb2) == 2.0
                    for nb2 in graph.adjacency[c_idx]
                )
                if not c_has_other_double:
                    return True
    return False


def _get_imine_nitrogen(graph: MoleculeGraph, c_idx: int) -> int | None:
    """C=N-H または C=N-R の N インデックスを返す。なければ None。"""
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "N" and get_bond_order(graph, c_idx, nb_idx) == 2.0:
            has_h = any(get_atom(graph, n).symbol == "H" for n in graph.adjacency[nb_idx])
            if has_h:
                return nb_idx
            # N-置換イミン: C のみの隣接 (Phase 64) — イソシアネート除外
            n_heavy = [n for n in graph.adjacency[nb_idx]
                       if n != c_idx and get_atom(graph, n).symbol != "H"]
            if n_heavy and all(get_atom(graph, n).symbol == "C" for n in n_heavy):
                c_has_other_double = any(
                    nb2 != nb_idx and get_bond_order(graph, c_idx, nb2) == 2.0
                    for nb2 in graph.adjacency[c_idx]
                )
                if not c_has_other_double:
                    return nb_idx
    return None


def _is_ketoxime(graph: MoleculeGraph, c_idx: int) -> bool:
    """C=N-OH パターンで C の C 隣接が 2+ のケトオキシムかチェック。"""
    atom = get_atom(graph, c_idx)
    if atom.in_ring:
        return False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "N" and get_bond_order(graph, c_idx, nb_idx) == 2.0:
            for n_nb_idx in graph.adjacency[nb_idx]:
                n_nb = get_atom(graph, n_nb_idx)
                if n_nb.symbol == "O" and get_bond_order(graph, nb_idx, n_nb_idx) == 1.0:
                    if any(get_atom(graph, o_nb).symbol == "H" for o_nb in graph.adjacency[n_nb_idx]):
                        c_count = sum(1 for nb2 in graph.adjacency[c_idx]
                                      if get_atom(graph, nb2).symbol == "C")
                        return c_count >= 2
    return False


def _is_aldoxime(graph: MoleculeGraph, c_idx: int) -> bool:
    """C=N-OH パターンで C の C 隣接が 0 または 1 のアルドキシムかチェック。"""
    atom = get_atom(graph, c_idx)
    if atom.in_ring:
        return False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "N" and get_bond_order(graph, c_idx, nb_idx) == 2.0:
            for n_nb_idx in graph.adjacency[nb_idx]:
                n_nb = get_atom(graph, n_nb_idx)
                if n_nb.symbol == "O" and get_bond_order(graph, nb_idx, n_nb_idx) == 1.0:
                    if any(get_atom(graph, o_nb).symbol == "H" for o_nb in graph.adjacency[n_nb_idx]):
                        c_count = sum(1 for nb2 in graph.adjacency[c_idx]
                                      if get_atom(graph, nb2).symbol == "C")
                        return c_count < 2
    return False


def _get_oxime_nitrogen(graph: MoleculeGraph, c_idx: int) -> int | None:
    """オキシム C に二重結合した N のインデックスを返す。なければ None。"""
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "N" and get_bond_order(graph, c_idx, nb_idx) == 2.0:
            for n_nb_idx in graph.adjacency[nb_idx]:
                n_nb = get_atom(graph, n_nb_idx)
                if n_nb.symbol == "O" and get_bond_order(graph, nb_idx, n_nb_idx) == 1.0:
                    if any(get_atom(graph, o_nb).symbol == "H" for o_nb in graph.adjacency[n_nb_idx]):
                        return nb_idx
    return None


def _is_carbonate(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が RO-C(=O)-OR パターン（炭酸エステル）かチェック。
    炭酸エステル: 1 個の C=O と 2 個の C-O-R, C 直結なし。
    Phase 45: 環内 C は除外（環状炭酸エステルはラクトン系として処理）。"""
    if get_atom(graph, c_idx).in_ring:  # 環状炭酸エステルは対象外
        return False
    if _get_double_bonded_oxygen(graph, c_idx) is None:
        return False
    # カルボニル C に直結する炭素がないこと
    c_neighbors = [nb for nb in graph.adjacency[c_idx]
                   if get_atom(graph, nb).symbol == "C"]
    if c_neighbors:
        return False
    # 単結合 O が 2 個で各々に C が隣接すること
    single_o_with_c = []
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol != "O" or get_bond_order(graph, c_idx, nb_idx) != 1.0:
            continue
        if any(get_atom(graph, on).symbol == "C" for on in graph.adjacency[nb_idx]
               if on != c_idx):
            single_o_with_c.append(nb_idx)
    return len(single_o_with_c) == 2


def _is_anhydride(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C(=O)-O-C(=O) パターン（酸無水物）かチェック。環状は除く（Phase 127）。"""
    if get_atom(graph, c_idx).in_ring:
        return False
    if _get_double_bonded_oxygen(graph, c_idx) is None:
        return False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol != "O" or get_bond_order(graph, c_idx, nb_idx) != 1.0:
            continue
        for o_nb_idx in graph.adjacency[nb_idx]:
            if o_nb_idx == c_idx:
                continue
            o_nb = get_atom(graph, o_nb_idx)
            if o_nb.symbol == "C" and _get_double_bonded_oxygen(graph, o_nb_idx) is not None:
                return True
    return False


def _is_ester(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C(=O)-O-R パターン（エステル）かチェック。"""
    has_double_o = False
    has_single_o_c = False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol != "O":
            continue
        bo = get_bond_order(graph, c_idx, nb_idx)
        if bo == 2.0:
            has_double_o = True
        elif bo == 1.0:
            # この O に炭素が隣接する（O-R 部分）
            for o_nb in graph.adjacency[nb_idx]:
                if get_atom(graph, o_nb).symbol == "C" and o_nb != c_idx:
                    has_single_o_c = True
    return has_double_o and has_single_o_c


def _get_double_bonded_sulfur(graph: MoleculeGraph, c_idx: int) -> int | None:
    """C に二重結合している S のインデックスを返す。なければ None。"""
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "S" and get_bond_order(graph, c_idx, nb_idx) == 2.0:
            return nb_idx
    return None


def _has_double_bonded_sulfur(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C=S を持つか判定。"""
    return _get_double_bonded_sulfur(graph, c_idx) is not None


def _is_thioamide(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C(=S)-NR₂ パターン（チオアミド）かチェック。"""
    if _get_double_bonded_sulfur(graph, c_idx) is None:
        return False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "N" and not nb.in_ring:
            if get_bond_order(graph, c_idx, nb_idx) == 1.0:
                return True
    return False


def _get_thioamide_nitrogen(graph: MoleculeGraph, c_idx: int) -> int | None:
    """チオアミド C に結合した N インデックスを返す。"""
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "N" and not nb.in_ring:
            if get_bond_order(graph, c_idx, nb_idx) == 1.0:
                return nb_idx
    return None


def _is_carbamate(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が N-C(=O)-O-R パターン（カルバメート）かチェック。"""
    has_double_o = False
    has_single_o_c = False
    has_single_n = False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "O":
            bo = get_bond_order(graph, c_idx, nb_idx)
            if bo == 2.0:
                has_double_o = True
            elif bo == 1.0:
                for o_nb in graph.adjacency[nb_idx]:
                    if get_atom(graph, o_nb).symbol == "C" and o_nb != c_idx:
                        has_single_o_c = True
        elif nb.symbol == "N" and get_bond_order(graph, c_idx, nb_idx) == 1.0:
            has_single_n = True
    return has_double_o and has_single_o_c and has_single_n


def _get_imine_double_bonded_nitrogen(graph: MoleculeGraph, c_idx: int) -> int | None:
    """C に二重結合している N のインデックスを返す（ヒドラゾン検出用）。"""
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "N" and get_bond_order(graph, c_idx, nb_idx) == 2.0:
            return nb_idx
    return None


def _is_kethydrazone(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C(=N-NH₂) パターン（ケトヒドラゾン）かチェック。
    ケトン由来: C に 2 個以上の C 隣接または C+H。"""
    n1_idx = _get_imine_double_bonded_nitrogen(graph, c_idx)
    if n1_idx is None:
        return False
    # N1 の隣に別の N があること（N1-N2 は単結合; ジアゾ C=[N+]=[N-] は除外）
    for nb in graph.adjacency[n1_idx]:
        if nb == c_idx:
            continue
        if get_atom(graph, nb).symbol == "N":
            if get_bond_order(graph, n1_idx, nb) != 1.0:
                continue
            # C の C 隣接が 2 個以上 → ケトン由来
            c_nbrs = [n for n in graph.adjacency[c_idx]
                      if get_atom(graph, n).symbol == "C"]
            h_nbrs = [n for n in graph.adjacency[c_idx]
                      if get_atom(graph, n).symbol == "H"]
            return len(c_nbrs) >= 2 or (len(c_nbrs) == 1 and len(h_nbrs) == 0)
    return False


def _is_aldhydrazone(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C(=N-NH₂) パターン（アルドヒドラゾン）かチェック。
    アルデヒド由来: C に H が 1 個付く。"""
    n1_idx = _get_imine_double_bonded_nitrogen(graph, c_idx)
    if n1_idx is None:
        return False
    for nb in graph.adjacency[n1_idx]:
        if nb == c_idx:
            continue
        if get_atom(graph, nb).symbol == "N":
            if get_bond_order(graph, n1_idx, nb) != 1.0:
                continue
            h_nbrs = [n for n in graph.adjacency[c_idx]
                      if get_atom(graph, n).symbol == "H"]
            c_nbrs = [n for n in graph.adjacency[c_idx]
                      if get_atom(graph, n).symbol == "C"]
            return len(h_nbrs) >= 1 and len(c_nbrs) <= 1
    return False


def _get_hydrazone_nitrogen(graph: MoleculeGraph, c_idx: int) -> int | None:
    """ヒドラゾン C に二重結合している N インデックスを返す。"""
    return _get_imine_double_bonded_nitrogen(graph, c_idx)


def _is_chloroformate(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が X-C(=O)-O-R パターン（ハロホルメート）かチェック (Phase 60/65)。"""
    if _get_double_bonded_oxygen(graph, c_idx) is None:
        return False
    has_halide = False
    has_ester_o = False
    has_c_neighbor = False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol in ("F", "Cl", "Br", "I"):
            has_halide = True
        elif nb.symbol == "O" and get_bond_order(graph, c_idx, nb_idx) == 1.0:
            o_c_nbrs = [n for n in graph.adjacency[nb_idx]
                        if n != c_idx and get_atom(graph, n).symbol == "C"]
            if o_c_nbrs:
                has_ester_o = True
        elif nb.symbol == "C":
            has_c_neighbor = True
    return has_halide and has_ester_o and not has_c_neighbor


def _is_acid_halide(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C(=O)-X パターン（酸ハライド）かチェック。"""
    has_double_o = _get_double_bonded_oxygen(graph, c_idx) is not None
    if not has_double_o:
        return False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol in ("F", "Cl", "Br", "I"):
            return True
    return False


def _is_amide(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C(=O)-NR₂ パターン（一・二・三級アミド）かチェック。"""
    has_double_o = _get_double_bonded_oxygen(graph, c_idx) is not None
    if not has_double_o:
        return False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        # 環外 N に単結合 (ニトリルの三重結合は除外)
        if nb.symbol == "N" and not nb.in_ring:
            if get_bond_order(graph, c_idx, nb_idx) == 1.0:
                return True
    return False


def _get_amide_nitrogen(graph: MoleculeGraph, c_idx: int) -> int | None:
    """アミド C に結合した N インデックスを返す（全次数対応）。"""
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "N" and not nb.in_ring:
            if get_bond_order(graph, c_idx, nb_idx) == 1.0:
                return nb_idx
    return None


def aggregate_groups(groups: list[FunctionalGroup]) -> list[FunctionalGroup]:
    """
    最高優先度の同一 group_type が複数ある場合、diol/dione/dioic_acid 等に集約する。

    例: [alcohol, alcohol, alkane] → [diol, alkane]
         [carboxylic_acid, carboxylic_acid] → [dioic_acid]
    """
    if not groups:
        return groups

    top_type = groups[0].group_type
    top_priority = groups[0].priority
    top_groups = [g for g in groups if g.group_type == top_type]

    if len(top_groups) <= 1:
        return groups

    _multi_map: dict[tuple[str, int], str] = {
        ("carboxylic_acid", 2): "dioic_acid",
        ("alcohol", 2): "diol",
        ("alcohol", 3): "triol",
        ("ketone", 2): "dione",
        ("aldehyde", 2): "dial",
        ("ester", 2): "diester",
        ("amine", 2): "diamine",
        ("amine", 3): "triamine",
        ("acid_halide", 2): "diacid_halide",
        ("nitrile", 2): "dinitrile",
        ("carboxylate", 2): "dicarboxylate",
    }
    multi_type = _multi_map.get((top_type, len(top_groups)))
    if multi_type is None:
        return groups  # 3+ ketones 等の未対応ケースはそのまま返す

    seen: set[int] = set()
    merged_atoms: list[int] = []
    for g in top_groups:
        for ai in g.atom_indices:
            if ai not in seen:
                seen.add(ai)
                merged_atoms.append(ai)

    merged = FunctionalGroup(
        group_type=multi_type,
        atom_indices=merged_atoms,
        priority=top_priority,
    )
    others = [g for g in groups if g.group_type != top_type]
    return [merged] + others


def principal_group(groups: list[FunctionalGroup]) -> FunctionalGroup | None:
    """官能基リストの最優先グループを返す（alkane の場合は None）。"""
    if not groups:
        return None
    top = groups[0]
    if top.group_type == "alkane":
        return None
    return top


# ─── Phase 68–73: 新規ヘルパー関数 ──────────────────────────────────────

def _is_carbamic_acid(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が RnN-C(=O)-OH パターン（カルバミン酸）かチェック (n=0,1,2)。"""
    has_imine_o = False
    has_n = False
    has_oh = False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        bo = get_bond_order(graph, c_idx, nb_idx)
        if nb.symbol == "O" and bo == 2.0:
            has_imine_o = True
        elif nb.symbol == "O" and bo == 1.0:
            if any(get_atom(graph, n).symbol == "H" for n in graph.adjacency[nb_idx]):
                has_oh = True
        elif nb.symbol == "N" and bo == 1.0:
            has_n = True
    return has_imine_o and has_n and has_oh


def _get_carbamic_oh(graph: MoleculeGraph, c_idx: int) -> int | None:
    """カルバミン酸 C の OH 酸素インデックスを返す。"""
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "O" and get_bond_order(graph, c_idx, nb_idx) == 1.0:
            if any(get_atom(graph, n).symbol == "H" for n in graph.adjacency[nb_idx]):
                return nb_idx
    return None


def _is_amidine(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C(=N-H)(N-H) パターン（アミジン）かチェック。"""
    has_imine_n = False
    has_amine_n = False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol != "N":
            continue
        bo = get_bond_order(graph, c_idx, nb_idx)
        n_h = sum(1 for n in graph.adjacency[nb_idx] if get_atom(graph, n).symbol == "H")
        if bo == 2.0 and n_h >= 1:
            has_imine_n = True
        elif bo == 1.0 and n_h >= 1:
            has_amine_n = True
    return has_imine_n and has_amine_n


def _get_amidine_nitrogens(
    graph: MoleculeGraph, c_idx: int
) -> tuple[int | None, int | None]:
    """アミジン C の (=N-H, -N-H) インデックスのペアを返す。"""
    n_imine: int | None = None
    n_amine: int | None = None
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol != "N":
            continue
        bo = get_bond_order(graph, c_idx, nb_idx)
        n_h = sum(1 for n in graph.adjacency[nb_idx] if get_atom(graph, n).symbol == "H")
        if bo == 2.0 and n_h >= 1:
            n_imine = nb_idx
        elif bo == 1.0 and n_h >= 1:
            n_amine = nb_idx
    return n_imine, n_amine


def _is_cyanate(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が O-C≡N パターン（シアン酸エステル）かチェック。"""
    has_triple_n = False
    has_single_o = False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        bo = get_bond_order(graph, c_idx, nb_idx)
        if nb.symbol == "N" and bo == 3.0:
            n_heavy = [n for n in graph.adjacency[nb_idx]
                       if n != c_idx and get_atom(graph, n).symbol != "H"]
            if not n_heavy:
                has_triple_n = True
        elif nb.symbol == "O" and bo == 1.0:
            has_single_o = True
    return has_triple_n and has_single_o


def _get_cyanate_oxygen(graph: MoleculeGraph, c_idx: int) -> int | None:
    """シアン酸エステル C の O インデックスを返す。"""
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "O" and get_bond_order(graph, c_idx, nb_idx) == 1.0:
            return nb_idx
    return None


def _is_thiocyanate(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が S-C≡N パターン（チオシアン酸エステル）かチェック。"""
    has_triple_n = False
    has_single_s = False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        bo = get_bond_order(graph, c_idx, nb_idx)
        if nb.symbol == "N" and bo == 3.0:
            n_heavy = [n for n in graph.adjacency[nb_idx]
                       if n != c_idx and get_atom(graph, n).symbol != "H"]
            if not n_heavy:
                has_triple_n = True
        elif nb.symbol == "S" and bo == 1.0:
            has_single_s = True
    return has_triple_n and has_single_s


def _get_thiocyanate_sulfur(graph: MoleculeGraph, c_idx: int) -> int | None:
    """チオシアン酸エステル C の S インデックスを返す。"""
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "S" and get_bond_order(graph, c_idx, nb_idx) == 1.0:
            return nb_idx
    return None


def _is_carbodiimide(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が R-N=C=N-R パターン（カルボジイミド）かチェック。"""
    if get_atom(graph, c_idx).in_ring:
        return False
    n_double = [nb for nb in graph.adjacency[c_idx]
                if get_atom(graph, nb).symbol == "N"
                and get_bond_order(graph, c_idx, nb) == 2.0]
    if len(n_double) != 2:
        return False
    # 各 N が外部 C (single bond) を持つことを確認
    for n_idx in n_double:
        c_sgl = [nb for nb in graph.adjacency[n_idx]
                 if nb != c_idx and get_atom(graph, nb).symbol == "C"
                 and get_bond_order(graph, n_idx, nb) == 1.0]
        if not c_sgl:
            return False
    return True


def _is_semicarbazone_or_thio(
    graph: MoleculeGraph, c_idx: int
) -> str | None:
    """C が C=N-N-C(=O)-NH₂ (セミカルバゾン) または C=N-N-C(=S)-NH₂ (チオセミカルバゾン) かチェック。
    C1 (imine 側) から検出。戻り値は group_type 文字列 or None。
    """
    # C1 に N への二重結合があるか
    n1_idx = _get_imine_double_bonded_nitrogen(graph, c_idx)
    if n1_idx is None:
        return None

    # N1 の隣の N2 (単結合)
    n2_candidates = [nb for nb in graph.adjacency[n1_idx]
                     if nb != c_idx and get_atom(graph, nb).symbol == "N"
                     and get_bond_order(graph, n1_idx, nb) == 1.0]
    if not n2_candidates:
        return None
    n2_idx = n2_candidates[0]

    # N2 の隣の C2 (単結合)
    c2_candidates = [nb for nb in graph.adjacency[n2_idx]
                     if nb != n1_idx and get_atom(graph, nb).symbol == "C"
                     and get_bond_order(graph, n2_idx, nb) == 1.0]
    if not c2_candidates:
        return None
    c2_idx = c2_candidates[0]

    # C2 に =O または =S があり、かつ NH₂ または NH がある
    has_carbonyl = _get_double_bonded_oxygen(graph, c2_idx) is not None
    has_thio = _get_double_bonded_sulfur(graph, c2_idx) is not None
    if not (has_carbonyl or has_thio):
        return None
    has_n_h = any(
        get_atom(graph, nb).symbol == "N"
        and get_bond_order(graph, c2_idx, nb) == 1.0
        and any(get_atom(graph, nnh).symbol == "H" for nnh in graph.adjacency[nb])
        for nb in graph.adjacency[c2_idx]
    )
    if not has_n_h:
        return None

    # C1 が H を持つ (アルデヒド型) か確認
    h_count = sum(1 for nb in graph.adjacency[c_idx] if get_atom(graph, nb).symbol == "H")
    c_count = sum(1 for nb in graph.adjacency[c_idx]
                  if nb != n1_idx and get_atom(graph, nb).symbol == "C")

    if has_carbonyl:
        return "aldsemicarbazone" if (h_count >= 1 and c_count <= 1) else "semicarbazone"
    else:
        return "aldthiosemicarbazone" if (h_count >= 1 and c_count <= 1) else "thiosemicarbazone"


def _is_peroxyacid(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C(=O)-O-O-H パターン（ペルオキシ酸）かチェック。
    ester より先に呼ぶこと。
    """
    has_double_o = _get_double_bonded_oxygen(graph, c_idx) is not None
    if not has_double_o:
        return False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol != "O" or get_bond_order(graph, c_idx, nb_idx) != 1.0:
            continue
        # この O1 に O2 が隣接する
        o2_candidates = [n for n in graph.adjacency[nb_idx]
                         if n != c_idx and get_atom(graph, n).symbol == "O"]
        if not o2_candidates:
            continue
        o2_idx = o2_candidates[0]
        # O2 に H があること
        if any(get_atom(graph, n).symbol == "H" for n in graph.adjacency[o2_idx]):
            return True
    return False


def _is_hydrazide(graph: MoleculeGraph, c_idx: int) -> bool:
    """C が C(=O)-N-N パターン（ヒドラジド）かチェック。
    amide より先に呼ぶこと。
    アジド (N=N double bond) とセミカルバゾン (N2 が C=N を持つ) は除外。
    """
    has_double_o = _get_double_bonded_oxygen(graph, c_idx) is not None
    if not has_double_o:
        return False
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol != "N" or get_bond_order(graph, c_idx, nb_idx) != 1.0:
            continue
        # N1 (amide-side N). 隣の N2 を探す。
        for n2_idx in graph.adjacency[nb_idx]:
            if n2_idx == c_idx:
                continue
            n2 = get_atom(graph, n2_idx)
            if n2.symbol != "N":
                continue
            # N1-N2 が単結合であること (azide は N=N double bond → 除外)
            if get_bond_order(graph, nb_idx, n2_idx) != 1.0:
                continue
            # N2 が C への二重結合を持たないこと (semicarbazone は N2=C → 除外)
            n2_dbl_c = any(
                get_atom(graph, nb3).symbol == "C"
                and get_bond_order(graph, n2_idx, nb3) == 2.0
                for nb3 in graph.adjacency[n2_idx] if nb3 != nb_idx
            )
            if n2_dbl_c:
                continue
            return True
    return False
