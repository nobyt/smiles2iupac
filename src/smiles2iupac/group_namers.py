"""
group_namers.py: _name_* ハンドラ関数と PGRP_DISPATCH テーブル。

__init__.py から抽出した命名ハンドラ群。
"""

from __future__ import annotations

from .chain_finder import (
    collect_acid_chain as _collect_acid_chain,
    chain_through_pivot as _chain_through_pivot,
    chain_multiple_bonds as _chain_multiple_bonds,
)

def _name_diester(graph, pgrp, get_atom) -> str:
    """
    ジエステル命名: di{alkyl} {acid}dioate 等 (Phase 51)
    例: CCOC(=O)CC(=O)OCC → diethyl propanedioate
    """
    from .molecule_analyzer import get_bond_order
    from .substituent import _name_carbon_substituent
    from .constants import CHAIN_PREFIX, MULTIPLIER

    # atom_indices: [c1, ..., c2, ...] の順 (2つのカルボニル C を含む)
    # カルボニル C を特定
    carbonyl_cs: list[int] = []
    for ai in pgrp.atom_indices:
        atom = get_atom(graph, ai)
        if atom.symbol != "C":
            continue
        for nb_idx in graph.adjacency[ai]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "O" and get_bond_order(graph, ai, nb_idx) == 2.0:
                carbonyl_cs.append(ai)
                break
    if len(carbonyl_cs) < 2:
        return "diester"

    c1, c2 = carbonyl_cs[0], carbonyl_cs[1]

    # 各エステル O → アルキル基
    def _get_ester_alkyl(carbonyl_c: int) -> str:
        for nb_idx in graph.adjacency[carbonyl_c]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "O" and get_bond_order(graph, carbonyl_c, nb_idx) == 1.0:
                alkyl_cs = [
                    n for n in graph.adjacency[nb_idx]
                    if n != carbonyl_c and get_atom(graph, n).symbol == "C"
                ]
                if alkyl_cs:
                    return _name_carbon_substituent(graph, alkyl_cs[0], {carbonyl_c, nb_idx})
                return "methyl"
        return "methyl"

    alkyl1 = _get_ester_alkyl(c1)
    alkyl2 = _get_ester_alkyl(c2)

    # 酸鎖 (c1 〜 c2 を結ぶ炭素鎖) を BFS で収集
    ester_os: set[int] = set()
    for ci in (c1, c2):
        for nb_idx in graph.adjacency[ci]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "O" and get_bond_order(graph, ci, nb_idx) == 1.0:
                ester_os.add(nb_idx)
                for nb2 in graph.adjacency[nb_idx]:
                    if nb2 != ci and get_atom(graph, nb2).symbol == "C":
                        ester_os.add(nb2)

    # Phase 268/320: 芳香族/シクロアルカン環に両カルボニルが直結の場合
    ring1, rc1, arom1 = _ring_for_exo_c(graph, c1, get_atom)
    ring2, rc2, arom2 = _ring_for_exo_c(graph, c2, get_atom)
    if ring1 is not None and ring2 is not None and ring1 == ring2:
        ring_list = sorted(ring1)
        from .ring_handler import _assign_ring_locants
        ring_chain = _assign_ring_locants(graph, ring_list, arom1, "arene" if arom1 else "alkane",
                                          [rc1, rc2])
        loc1 = ring_chain.locant_map.get(rc1, 1)
        loc2 = ring_chain.locant_map.get(rc2, 2)
        locs = sorted([loc1, loc2])
        n = len(ring1)
        ring_base = "benzene" if (arom1 and n == 6) else f"cyclo{CHAIN_PREFIX.get(n, f'C{n}')}ane"
        acid_name = f"{ring_base}-{locs[0]},{locs[1]}-dicarboxylate"
        if alkyl1 == alkyl2:
            mult = MULTIPLIER.get(2, "di")
            return f"{mult}{alkyl1} {acid_name}"
        return f"{min(alkyl1, alkyl2)} {max(alkyl1, alkyl2)} {acid_name}"

    chain_fwd_de = _collect_acid_chain(graph, c1, ester_os, get_atom)
    chain_rev_de = _collect_acid_chain(graph, c2, ester_os, get_atom)
    ene_fwd_de, yne_fwd_de = _chain_multiple_bonds(graph, chain_fwd_de)
    ene_rev_de, yne_rev_de = _chain_multiple_bonds(graph, chain_rev_de)
    mb_fwd_de = sorted(ene_fwd_de + yne_fwd_de)
    mb_rev_de = sorted(ene_rev_de + yne_rev_de)
    if mb_rev_de and (not mb_fwd_de or mb_rev_de < mb_fwd_de):
        acid_carbons = chain_rev_de
        _ene_de, _yne_de = ene_rev_de, yne_rev_de
    else:
        acid_carbons = chain_fwd_de
        _ene_de, _yne_de = ene_fwd_de, yne_fwd_de
    n_acid = len(acid_carbons)

    stem = CHAIN_PREFIX.get(n_acid, f"C{n_acid}")
    if _ene_de or _yne_de:
        from .name_assembler import _format_multiple_bonds as _fmt_de
        acid_name = f"{stem}{_fmt_de(_ene_de, _yne_de)}edioate"
    else:
        _DIACID_RETAINED = {
            2: "oxalate", 3: "malonate", 4: "succinate", 5: "glutarate",
            6: "adipate", 7: "pimelate", 8: "suberate", 9: "azelate", 10: "sebacate",
        }
        acid_name = _DIACID_RETAINED.get(n_acid, f"{stem}anedioate")

    stereo_pfx_de = ""
    if _ene_de:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_de = PrincipalChain(atom_indices=acid_carbons,
                                locant_map={c: i + 1 for i, c in enumerate(acid_carbons)})
        _stereo_de = assign_stereochemistry(graph, _pc_de)
        if _stereo_de:
            stereo_pfx_de = "(" + ",".join(d.strip("()") for d in _stereo_de) + ")-"

    # alkyl 部分の組み立て
    if alkyl1 == alkyl2:
        mult = MULTIPLIER.get(2, "di")
        return f"{stereo_pfx_de}{mult}{alkyl1} {acid_name}"
    return f"{stereo_pfx_de}{min(alkyl1, alkyl2)} {max(alkyl1, alkyl2)} {acid_name}"


def _name_dioic_acid(graph, pgrp, get_atom):
    """環状ジカルボン酸命名: cycloalkane-X,Y-dicarboxylic acid (Phase 311)"""
    from .molecule_analyzer import get_bond_order
    from .constants import CHAIN_PREFIX

    pgrp_atoms = set(pgrp.atom_indices)
    cooh_cs = []
    for ai in pgrp.atom_indices:
        if get_atom(graph, ai).symbol != "C":
            continue
        for nb_idx in graph.adjacency[ai]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "O" and get_bond_order(graph, ai, nb_idx) == 2.0:
                cooh_cs.append(ai)
                break

    if len(cooh_cs) < 2:
        return None

    c1, c2 = cooh_cs[0], cooh_cs[1]

    def _ring_attachment(c_idx: int):
        for nb_idx in graph.adjacency[c_idx]:
            if nb_idx in pgrp_atoms:
                continue
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "C" and nb.in_ring:
                ring = next(
                    (rt for rt in (graph.ring_atom_sets or []) if nb_idx in rt), None
                )
                return nb_idx, ring
        return None, None

    ring_c1, ring1 = _ring_attachment(c1)
    ring_c2, ring2 = _ring_attachment(c2)

    if ring1 is None or ring2 is None or set(ring1) != set(ring2):
        return None

    ring_set = set(ring1)
    ring_list = sorted(ring_set)
    ring_size = len(ring_list)

    is_benzene = (ring_size == 6
                  and all(get_atom(graph, a).symbol == "C" and get_atom(graph, a).is_aromatic
                          for a in ring_list))

    from .ring_handler import _assign_ring_locants
    ring_chain = _assign_ring_locants(graph, ring_list, is_benzene,
                                      "arene" if is_benzene else "cycloalkane",
                                      [ring_c1, ring_c2])
    loc1 = ring_chain.locant_map.get(ring_c1, 1)
    loc2 = ring_chain.locant_map.get(ring_c2, 2)
    locs = sorted([loc1, loc2])

    if is_benzene:
        return f"benzene-{locs[0]},{locs[1]}-dicarboxylic acid"

    # Phase 523: heteroaromatic rings (furan, thiophene, pyridine, ...)
    is_aromatic_ring = all(get_atom(graph, a).is_aromatic for a in ring_list)
    has_heteroatom = any(get_atom(graph, a).symbol != "C" for a in ring_list)
    if is_aromatic_ring and has_heteroatom:
        # Get ring atoms in connectivity order via BFS through ring bonds
        ring_set_bfs = set(ring_list)
        ordered: list[int] = [ring_list[0]]
        seen: set[int] = {ring_list[0]}
        while len(ordered) < len(ring_list):
            cur = ordered[-1]
            moved = False
            for nb in graph.adjacency[cur]:
                if nb in ring_set_bfs and nb not in seen:
                    ordered.append(nb)
                    seen.add(nb)
                    moved = True
                    break
            if not moved:
                break
        from .heterocycle_handler import _find_best_start, _canonical_sig, _RETAINED_NAMES
        rotation = _find_best_start(ordered, graph)
        sig = _canonical_sig(rotation, graph)
        entry = _RETAINED_NAMES.get((True, sig))
        if entry is not None:
            het_name, has_ind_h = entry
            ind_h = "1H-" if has_ind_h else ""
            # Try both ring directions; pick the one giving the lowest locants.
            rev_rotation = [rotation[0]] + list(reversed(rotation[1:]))
            best_locs: list[int] | None = None
            for rot in (rotation, rev_rotation):
                loc_map_het = {atom: i + 1 for i, atom in enumerate(rot)}
                h_l1 = loc_map_het.get(ring_c1)
                h_l2 = loc_map_het.get(ring_c2)
                if h_l1 is not None and h_l2 is not None:
                    pair = sorted([h_l1, h_l2])
                    if best_locs is None or pair < best_locs:
                        best_locs = pair
            if best_locs is not None:
                return f"{ind_h}{het_name}-{best_locs[0]},{best_locs[1]}-dicarboxylic acid"

    stem = CHAIN_PREFIX.get(ring_size, f"C{ring_size}")
    return f"cyclo{stem}ane-{locs[0]},{locs[1]}-dicarboxylic acid"


def _name_dicarboxylate(graph, pgrp, get_atom) -> str:
    """ジカルボキシレートジアニオン: butanedioate = succinate  (Phase 146)"""
    from .constants import CHAIN_PREFIX
    from .molecule_analyzer import get_bond_order
    # 両端の carboxylate C を収集（直接 O 隣接でカルボニルを持つ C）
    c_idxs = []
    for ai in pgrp.atom_indices:
        if get_atom(graph, ai).symbol == "C":
            has_co = any(get_atom(graph, nb).symbol == "O" and
                         get_bond_order(graph, ai, nb) == 2.0
                         for nb in graph.adjacency[ai])
            if has_co:
                c_idxs.append(ai)
    if not c_idxs:
        return "dioate"
    c1 = c_idxs[0]
    c2 = c_idxs[1] if len(c_idxs) >= 2 else c1
    chain_fwd_dc = _collect_acid_chain(graph, c1, set(), get_atom)
    chain_rev_dc = _collect_acid_chain(graph, c2, set(), get_atom)
    ene_fwd_dc, yne_fwd_dc = _chain_multiple_bonds(graph, chain_fwd_dc)
    ene_rev_dc, yne_rev_dc = _chain_multiple_bonds(graph, chain_rev_dc)
    mb_fwd_dc = sorted(ene_fwd_dc + yne_fwd_dc)
    mb_rev_dc = sorted(ene_rev_dc + yne_rev_dc)
    if mb_rev_dc and (not mb_fwd_dc or mb_rev_dc < mb_fwd_dc):
        chain = chain_rev_dc
        _ene_dc, _yne_dc = ene_rev_dc, yne_rev_dc
    else:
        chain = chain_fwd_dc
        _ene_dc, _yne_dc = ene_fwd_dc, yne_fwd_dc
    n = len(chain)
    if _ene_dc or _yne_dc:
        from .name_assembler import _format_multiple_bonds as _fmt_dc
        stem = CHAIN_PREFIX.get(n, f"C{n}")
        base = f"{stem}{_fmt_dc(_ene_dc, _yne_dc)}edioate"
        if _ene_dc:
            from .stereochemistry import assign_stereochemistry
            from .chain_finder import PrincipalChain
            _pc_dc = PrincipalChain(atom_indices=chain,
                                    locant_map={c: i + 1 for i, c in enumerate(chain)})
            _stereo_dc = assign_stereochemistry(graph, _pc_dc)
            if _stereo_dc:
                _comb_dc = ",".join(d.strip("()") for d in _stereo_dc)
                return f"({_comb_dc})-{base}"
        return base
    # 保留ジアニオン名 (飽和のみ)
    retained = {2: "oxalate", 3: "malonate", 4: "succinate", 5: "glutarate",
                6: "adipate", 7: "pimelate", 8: "suberate", 9: "azelate", 10: "sebacate"}
    if n in retained:
        return retained[n]
    stem = CHAIN_PREFIX.get(n, f"C{n}")
    return f"{stem}anedioate"


def _name_carboxylate(graph, pgrp, get_atom) -> str:
    """カルボキシレートアニオン: ethanoate, formate, acetate  (Phase 146)"""
    from .constants import CHAIN_PREFIX
    carbonyl_c = pgrp.atom_indices[0]
    # 芳香族環に直接結合したカルボキシレート → benzoate (Phase 265)
    o_idxs = {nb for nb in graph.adjacency[carbonyl_c] if get_atom(graph, nb).symbol == "O"}
    for nb_idx in graph.adjacency[carbonyl_c]:
        if nb_idx in o_idxs:
            continue
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "C" and nb.is_aromatic and nb.in_ring:
            ring_atoms = next(
                (set(rt) for rt in (graph.ring_atom_sets or []) if nb_idx in rt), set()
            )
            if (len(ring_atoms) == 6
                    and all(get_atom(graph, a).symbol == "C" for a in ring_atoms)):
                return "benzoate"
    acid_chain = _collect_acid_chain(graph, carbonyl_c, o_idxs, get_atom)
    n = len(acid_chain)
    if n == 1:
        return "formate"
    if n == 2:
        return "acetate"
    stem = CHAIN_PREFIX.get(n, f"C{n}")
    ene, yne = _chain_multiple_bonds(graph, acid_chain)
    if ene or yne:
        from .name_assembler import _format_multiple_bonds as _fmt
        base_name = f"{stem}{_fmt(ene, yne)}oate"
        if ene:
            from .stereochemistry import assign_stereochemistry
            from .chain_finder import PrincipalChain
            _pc_cox = PrincipalChain(atom_indices=acid_chain,
                                     locant_map={c: i + 1 for i, c in enumerate(acid_chain)})
            _stereo_cox = assign_stereochemistry(graph, _pc_cox)
            if _stereo_cox:
                _comb_cox = ",".join(d.strip("()") for d in _stereo_cox)
                return f"({_comb_cox})-{base_name}"
        return base_name
    return f"{stem}anoate"


def _name_thioic_acid(graph, pgrp, get_atom) -> str:
    """チオカルボン酸命名 (Phase 149)"""
    from .constants import CHAIN_PREFIX
    gtype = pgrp.group_type
    carbonyl_c = pgrp.atom_indices[0]
    chalcogen_idxs = {nb for nb in graph.adjacency[carbonyl_c]
                      if get_atom(graph, nb).symbol in ("O", "S")}
    # 芳香族環に直接結合したカルボニル → benzenecarbothioic (Phase 265)
    for nb_idx in graph.adjacency[carbonyl_c]:
        if nb_idx in chalcogen_idxs:
            continue
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "C" and nb.is_aromatic and nb.in_ring:
            ring_atoms = next(
                (set(rt) for rt in (graph.ring_atom_sets or []) if nb_idx in rt), set()
            )
            if (len(ring_atoms) == 6
                    and all(get_atom(graph, a).symbol == "C" for a in ring_atoms)):
                if gtype == "thioic_s_acid":
                    return "benzenecarbothioic S-acid"
                elif gtype == "thioic_o_acid":
                    return "benzenecarbothioic O-acid"
                else:
                    return "benzenecarbodithioic acid"
            # Heteroaromatic ring
            if (any(get_atom(graph, a).symbol != "C" for a in ring_atoms)
                    and all(get_atom(graph, a).is_aromatic for a in ring_atoms)):
                _apfx_th = _aryl_sulfonyl_prefix(graph, nb_idx, carbonyl_c, get_atom)
                if _apfx_th is not None:
                    if gtype == "thioic_s_acid":
                        return f"{_apfx_th}carbothioic S-acid"
                    elif gtype == "thioic_o_acid":
                        return f"{_apfx_th}carbothioic O-acid"
                    else:
                        return f"{_apfx_th}carbodithioic acid"
    acid_chain = _collect_acid_chain(graph, carbonyl_c, chalcogen_idxs, get_atom)
    n = len(acid_chain)
    stem = CHAIN_PREFIX.get(n, f"C{n}")
    ene, yne = _chain_multiple_bonds(graph, acid_chain)
    if ene or yne:
        from .name_assembler import _format_multiple_bonds as _fmt
        mb = _fmt(ene, yne)
    else:
        mb = ""
    stereo_pfx = ""
    if ene:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_thi = PrincipalChain(atom_indices=acid_chain,
                                 locant_map={c: i + 1 for i, c in enumerate(acid_chain)})
        _stereo_thi = assign_stereochemistry(graph, _pc_thi)
        if _stereo_thi:
            _comb_thi = ",".join(d.strip("()") for d in _stereo_thi)
            stereo_pfx = f"({_comb_thi})-"
    if gtype == "thioic_s_acid":
        if not mb:
            return f"{stereo_pfx}{stem}anethioic S-acid"
        return f"{stereo_pfx}{stem}{mb}ethioic S-acid"
    elif gtype == "thioic_o_acid":
        if not mb:
            return f"{stereo_pfx}{stem}anethioic O-acid"
        return f"{stereo_pfx}{stem}{mb}ethioic O-acid"
    else:  # dithioic_acid
        if not mb:
            return f"{stereo_pfx}{stem}anedithioic acid"
        return f"{stereo_pfx}{stem}{mb}edithioic acid"


def _name_nitrate_ester(graph, pgrp, get_atom) -> str:
    """ニトレートエステル命名: {alkyl} nitrate (Phase 150)"""
    from .substituent import _name_carbon_substituent
    from .molecule_analyzer import get_bond_order
    n_idx = pgrp.atom_indices[0]
    ester_o = pgrp.atom_indices[1]
    alkyl_cs = [nb for nb in graph.adjacency[ester_o]
                if nb != n_idx and get_atom(graph, nb).symbol == "C"]
    if alkyl_cs:
        alkyl_name = _name_carbon_substituent(graph, alkyl_cs[0], {ester_o, n_idx})
    else:
        alkyl_name = "methyl"
    return f"{alkyl_name} nitrate"


def _name_nitrite_ester(graph, pgrp, get_atom) -> str:
    """ニトライトエステル命名: {alkyl} nitrite (Phase 150)"""
    from .substituent import _name_carbon_substituent
    n_idx = pgrp.atom_indices[0]
    ester_o = pgrp.atom_indices[1]
    alkyl_cs = [nb for nb in graph.adjacency[ester_o]
                if nb != n_idx and get_atom(graph, nb).symbol == "C"]
    if alkyl_cs:
        alkyl_name = _name_carbon_substituent(graph, alkyl_cs[0], {ester_o, n_idx})
    else:
        alkyl_name = "methyl"
    return f"{alkyl_name} nitrite"


def _name_imidic_acid(graph, pgrp, get_atom) -> str:
    """
    イミド酸命名: {stem}imidic acid / N-{sub}{stem}imidic acid (E/Z 対応)
    例: CC(=N)O  → ethanimidic acid
        CC(=NC)O → N-methylethanimidic acid
    """
    from .constants import CHAIN_PREFIX, MULTIPLIER
    from .substituent import _name_carbon_substituent
    from collections import Counter as _Ctr_ia

    c_idx = pgrp.atom_indices[0]
    o_idxs = {nb for nb in graph.adjacency[c_idx] if get_atom(graph, nb).symbol == "O"}
    n_idxs = {nb for nb in graph.adjacency[c_idx] if get_atom(graph, nb).symbol == "N"}

    # Ring-attached imidic acid: benzene → "benzenecarboximidic acid",
    # heteroaromatic → "pyridine-2-carboximidic acid" etc.
    for _rn_ia in graph.adjacency[c_idx]:
        _rna_ia = get_atom(graph, _rn_ia)
        if not (_rna_ia.symbol == "C" and _rna_ia.in_ring and _rna_ia.is_aromatic):
            continue
        _apfx_ia = _aryl_sulfonyl_prefix(graph, _rn_ia, c_idx, get_atom)
        if _apfx_ia is None:
            continue
        _n_pfx_ia = ""
        for _ni_ia in n_idxs:
            _nsubs_ia = [nb for nb in graph.adjacency[_ni_ia]
                         if nb != c_idx and get_atom(graph, nb).symbol == "C"]
            if _nsubs_ia:
                _sn_ia = [_name_carbon_substituent(graph, nb, {_ni_ia}) for nb in _nsubs_ia]
                _cnt_ia = _Ctr_ia(_sn_ia)
                _pts_ia: list[str] = []
                for _s in sorted(_cnt_ia):
                    _c = _cnt_ia[_s]
                    _ss = f"({_s})" if _s.startswith("(") else _s
                    if _c == 1:
                        _pts_ia.append(f"N-{_ss}")
                    else:
                        _pts_ia.append(f"N,N-{MULTIPLIER.get(_c, str(_c))}{_ss}")
                _n_pfx_ia = "-".join(_pts_ia)
            break
        return f"{_n_pfx_ia}{_apfx_ia}carboximidic acid"

    acid_chain = _collect_acid_chain(graph, c_idx, o_idxs | n_idxs, get_atom)
    n_c = len(acid_chain)
    stem = CHAIN_PREFIX.get(n_c, f"C{n_c}")
    _ene_ia, _yne_ia = _chain_multiple_bonds(graph, acid_chain)
    if _ene_ia or _yne_ia:
        from .name_assembler import _format_multiple_bonds as _fmt_ia
        base_name = f"{stem}{_fmt_ia(_ene_ia, _yne_ia)}imidic acid"
    else:
        base_name = f"{stem}animidic acid"

    stereo_pfx_ia = ""
    if _ene_ia:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_ia = PrincipalChain(atom_indices=acid_chain,
                                locant_map={c: i + 1 for i, c in enumerate(acid_chain)})
        _stereo_ia = assign_stereochemistry(graph, _pc_ia)
        if _stereo_ia:
            stereo_pfx_ia = "(" + ",".join(d.strip("()") for d in _stereo_ia) + ")-"

    # N-置換基
    n_prefix_ia = ""
    for n_idx in n_idxs:
        n_subs_ia = [nb for nb in graph.adjacency[n_idx]
                     if nb != c_idx and get_atom(graph, nb).symbol == "C"]
        if n_subs_ia:
            from .constants import MULTIPLIER
            sub_names_ia = [_name_carbon_substituent(graph, nb, {n_idx}) for nb in n_subs_ia]
            cnt_map_ia = _Ctr_ia(sub_names_ia)
            parts_ia = []
            for sub in sorted(cnt_map_ia):
                cnt = cnt_map_ia[sub]
                sub_str_ia = f"({sub})" if sub.startswith("(") else sub
                if cnt == 1:
                    parts_ia.append(f"N-{sub_str_ia}")
                else:
                    mult = MULTIPLIER.get(cnt, f"{cnt}")
                    parts_ia.append(f"N,N-{mult}{sub_str_ia}")
            n_prefix_ia = "-".join(parts_ia)
        break  # only one N expected

    return f"{stereo_pfx_ia}{n_prefix_ia}{base_name}"


def _name_imidate_ester(graph, pgrp, get_atom) -> str | None:
    """
    イミダートエステル命名: {alkyl} {stem}imidate / {alkyl} N-{sub}{stem}imidate
    例: CC(=N)OCC    → ethyl ethanimidate
        CCOC(=NC)CC  → ethyl N-methylpropanimidate
    """
    from .constants import CHAIN_PREFIX
    from .substituent import _name_carbon_substituent

    c_idx = pgrp.atom_indices[0]
    n_idx = next((ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "N"), None)
    # ester O を特定
    o_idx = next((ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "O"), None)
    if o_idx is None:
        return None

    # Phase 508: C=N-...-O が全て同一環内にある場合は環命名系に委譲
    if (get_atom(graph, c_idx).in_ring
            and n_idx is not None and get_atom(graph, n_idx).in_ring
            and get_atom(graph, o_idx).in_ring):
        return None

    # alkyl (O 側)
    ester_c = next(
        (nb for nb in graph.adjacency[o_idx]
         if nb != c_idx and get_atom(graph, nb).symbol == "C"),
        None,
    )
    alkyl_name = _name_carbon_substituent(graph, ester_c, {o_idx}) if ester_c else "methyl"

    # acid chain (C(=N) 側、O を除く)
    acid_chain = _collect_acid_chain(graph, c_idx, {o_idx}, get_atom)
    n_c = len(acid_chain)
    stem = CHAIN_PREFIX.get(n_c, f"C{n_c}")
    _ene_im, _yne_im = _chain_multiple_bonds(graph, acid_chain)
    if _ene_im or _yne_im:
        from .name_assembler import _format_multiple_bonds as _fmt_im
        acid_name = f"{stem}{_fmt_im(_ene_im, _yne_im)}imidate"
    else:
        acid_name = f"{stem}animidate"

    stereo_pfx_im = ""
    if _ene_im:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_im = PrincipalChain(atom_indices=acid_chain,
                                locant_map={c: i + 1 for i, c in enumerate(acid_chain)})
        _stereo_im = assign_stereochemistry(graph, _pc_im)
        if _stereo_im:
            stereo_pfx_im = "(" + ",".join(d.strip("()") for d in _stereo_im) + ")-"

    # N-置換基
    n_prefix = ""
    if n_idx is not None:
        n_subs = [nb for nb in graph.adjacency[n_idx]
                  if nb != c_idx and get_atom(graph, nb).symbol == "C"]
        if n_subs:
            from .constants import MULTIPLIER
            from .name_assembler import _needs_bis_tris as _nbp_im
            from collections import Counter as _Ctr
            sub_names = [_name_carbon_substituent(graph, nb, {n_idx}) for nb in n_subs]
            cnt_map = _Ctr(sub_names)
            parts = []
            for sub in sorted(cnt_map):
                cnt = cnt_map[sub]
                sub_str = f"({sub})" if _nbp_im(sub) else sub
                if cnt == 1:
                    parts.append(f"N-{sub_str}")
                else:
                    mult = MULTIPLIER.get(cnt, f"{cnt}")
                    parts.append(f"N,N-{mult}{sub_str}")
            n_prefix = "-".join(parts)

    # Ring-attached imidate (benzene or heteroaromatic)
    for _rn_imid in graph.adjacency[c_idx]:
        _rna_imid = get_atom(graph, _rn_imid)
        if not (_rna_imid.symbol == "C" and _rna_imid.in_ring and _rna_imid.is_aromatic):
            continue
        _apfx_imid = _aryl_sulfonyl_prefix(graph, _rn_imid, c_idx, get_atom)
        if _apfx_imid is None:
            continue
        return f"{alkyl_name} {n_prefix}{_apfx_imid}carboximidate"

    return f"{alkyl_name} {stereo_pfx_im}{n_prefix}{acid_name}"


def _name_ester(graph, pgrp, get_atom) -> str:
    """
    エステル命名: (alkyl) (stem)oate
    例: CCOC(=O)C → ethyl ethanoate
    """
    from .molecule_analyzer import carbon_indices
    from .substituent import _name_carbon_substituent
    from .constants import CHAIN_PREFIX

    # pgrp.atom_indices[0] = 酸カルボニル C
    carbonyl_c = pgrp.atom_indices[0]

    # エステル O (単結合、R 側) と アルキル C を特定
    ester_o: int | None = None
    for nb_idx in graph.adjacency[carbonyl_c]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "O":
            from .molecule_analyzer import get_bond_order
            if get_bond_order(graph, carbonyl_c, nb_idx) == 1.0:
                ester_o = nb_idx
                break

    # アルキル基 (ester O の隣の C 群)
    alkyl_name = "methyl"  # デフォルト
    if ester_o is not None:
        alkyl_cs = [
            nb for nb in graph.adjacency[ester_o]
            if nb != carbonyl_c and get_atom(graph, nb).symbol == "C"
        ]
        if alkyl_cs:
            excluded = {carbonyl_c, ester_o}
            alkyl_name = _name_carbon_substituent(graph, alkyl_cs[0], excluded)
            # methyl→methyl, ethyl→ethyl, propyl→propyl ...
            # _name_carbon_substituent returns "methyl", "ethyl" etc. ✓

    # 酸鎖 (カルボニル C 側) — ester_o を除外した炭素鎖長
    ester_o_set = {ester_o} if ester_o else set()
    if ester_o:
        ester_o_set.update(
            nb for nb in graph.adjacency[ester_o]
            if nb != carbonyl_c and get_atom(graph, nb).symbol == "C"
        )

    # 酸側: カルボニル C の隣が芳香環 C の場合 → benzoate (Phase 47)
    acid_ring_c: int | None = None
    for nb_idx in graph.adjacency[carbonyl_c]:
        if nb_idx in ester_o_set:
            continue
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "C" and nb.is_aromatic and nb.in_ring:
            acid_ring_c = nb_idx
            break

    if acid_ring_c is not None:
        # Phase 155: 芳香族エステル → 環種別に命名
        # BFS で連結する環原子を収集（縮合環対策: 単一環のみ対応）
        ring_c_set: set[int] = set()
        _q_r: list[int] = [acid_ring_c]
        while _q_r:
            _a_r = _q_r.pop()
            if _a_r in ring_c_set or not get_atom(graph, _a_r).in_ring:
                continue
            ring_c_set.add(_a_r)
            _q_r.extend(graph.adjacency[_a_r])

        _ring_sz = len(ring_c_set)
        _all_arom = all(get_atom(graph, a).is_aromatic for a in ring_c_set)

        if _ring_sz in (5, 6) and _all_arom:
            _pure_benz = all(get_atom(graph, a).symbol == "C" for a in ring_c_set)
            # 環内隣接グラフ
            _radj: dict[int, list[int]] = {
                a: [nb for nb in graph.adjacency[a] if nb in ring_c_set]
                for a in ring_c_set
            }

            def _traverse_ring(start: int, nxt: int) -> list[int] | None:
                order = [start, nxt]
                prev, cur = start, nxt
                while len(order) < _ring_sz:
                    nexts = [nb for nb in _radj[cur] if nb != prev]
                    if not nexts:
                        return None
                    prev, cur = cur, nexts[0]
                    order.append(cur)
                return order if len(order) == _ring_sz else None

            _nbrs = _radj.get(acid_ring_c, [])

            if _pure_benz and len(_nbrs) >= 2:
                # ベンゼン環エステル → {prefix}benzoate
                best_order_b: list[int] | None = None
                best_slocs_b: list[int] | None = None
                for _nb0 in _nbrs:
                    _ord = _traverse_ring(acid_ring_c, _nb0)
                    if _ord is None:
                        continue
                    _lm = {a: i + 1 for i, a in enumerate(_ord)}
                    _slocs = sorted(
                        _lm[a] for a in ring_c_set
                        if a != acid_ring_c
                        and any(
                            nb not in ring_c_set
                            and get_atom(graph, nb).symbol != "H"
                            and nb != carbonyl_c
                            for nb in graph.adjacency[a]
                        )
                    )
                    if best_slocs_b is None or _slocs < best_slocs_b:
                        best_slocs_b = _slocs
                        best_order_b = _ord
                if best_order_b is not None:
                    _lm_b = {a: i + 1 for i, a in enumerate(best_order_b)}
                    from .substituent import name_substituent as _nsub_b
                    from .name_assembler import _build_prefix as _bpfx_b
                    _excl_b = ring_c_set | {carbonyl_c}
                    _rsubs_b: list[tuple[int, str]] = []
                    for _ai_b in ring_c_set:
                        if _ai_b == acid_ring_c:
                            continue
                        for _nb_b in graph.adjacency[_ai_b]:
                            if _nb_b in _excl_b or get_atom(graph, _nb_b).symbol == "H":
                                continue
                            _sn_b = _nsub_b(graph, _nb_b, _excl_b)
                            if _sn_b:
                                _rsubs_b.append((_lm_b[_ai_b], _sn_b))
                    if not _rsubs_b:
                        return f"{alkyl_name} benzoate"
                    return f"{alkyl_name} {_bpfx_b(_rsubs_b)}benzoate"

            elif not _pure_benz and len(_nbrs) >= 2:
                # ヘテロ芳香族環エステル → {ring}-{locant}-carboxylate
                from .heterocycle_handler import (
                    _find_best_start as _fbs_e,
                    _RETAINED_NAMES as _RN_e,
                    _canonical_sig as _csig_e,
                    _is_aromatic_ring as _iar_e,
                )
                _ring_list_e = list(ring_c_set)
                _rot_e = _fbs_e(_ring_list_e, graph)
                _lm_e = {a: i + 1 for i, a in enumerate(_rot_e)}
                _rloc_e = _lm_e.get(acid_ring_c)
                if _rloc_e is not None:
                    _is_arom_e = _iar_e(_ring_list_e, graph)
                    _csig_e_val = _csig_e(_rot_e, graph)
                    _retained_e = _RN_e.get((_is_arom_e, _csig_e_val))
                    if _retained_e:
                        _rbase_e, _nh_e = _retained_e
                        if _nh_e:
                            _rbase_e = f"1H-{_rbase_e}"
                        return f"{alkyl_name} {_rbase_e}-{_rloc_e}-carboxylate"

    # Phase 188: 非芳香族環に結合したカルボニル → cycloalkanecarboxylate 型
    for _nb188 in graph.adjacency[carbonyl_c]:
        if _nb188 in ester_o_set:
            continue
        _nb188a = get_atom(graph, _nb188)
        if _nb188a.symbol == "C" and _nb188a.in_ring and not _nb188a.is_aromatic:
            _ring_set188 = next(
                (rt for rt in (graph.ring_atom_sets or []) if _nb188 in rt), None
            )
            if (_ring_set188 is not None
                    and all(get_atom(graph, a).symbol == "C" for a in _ring_set188)):
                from .ring_handler import (
                    _assign_ring_locants,
                    collect_ring_substituents,
                    assemble_ring_name,
                )
                _ring_list188 = list(_ring_set188)
                _ring_chain188 = _assign_ring_locants(
                    graph, _ring_list188, False, "alkane", [carbonyl_c]
                )
                _ring_subs188 = collect_ring_substituents(
                    graph, _ring_chain188, [carbonyl_c]
                )
                _ring_base188 = assemble_ring_name(
                    _ring_chain188, _ring_subs188, "alkane", None, []
                )
                return f"{alkyl_name} {_ring_base188}carboxylate"
            break

    # 酸側の炭素鎖を DFS で収集（ester_o 方向は除外）
    acid_carbons = _collect_acid_chain(graph, carbonyl_c, ester_o_set, get_atom)
    n_acid = len(acid_carbons)

    stem = CHAIN_PREFIX.get(n_acid, f"C{n_acid}")
    _ene_e, _yne_e = _chain_multiple_bonds(graph, acid_carbons)
    if _ene_e or _yne_e:
        from .name_assembler import _format_multiple_bonds as _fmt_e
        acid_name = f"{stem}{_fmt_e(_ene_e, _yne_e)}oate"
    elif n_acid == 1:
        # Phase 120: ギ酸エステルは "formate" (IUPAC 2013 PIN)
        acid_name = "formate"
    elif n_acid == 2:
        # Phase 117: 酢酸エステルは "acetate" (IUPAC 2013 PIN)
        acid_name = "acetate"
    else:
        acid_name = f"{stem}anoate"

    # 鎖上の置換基を収集 (Phase 72)
    from .substituent import collect_substituents
    from .name_assembler import _build_prefix
    locant_map = {c: i + 1 for i, c in enumerate(acid_carbons)}
    principal_atoms = set(pgrp.atom_indices) | ester_o_set
    for nb_idx in graph.adjacency[carbonyl_c]:
        if get_atom(graph, nb_idx).symbol == "O":
            principal_atoms.add(nb_idx)
    subs = collect_substituents(graph, acid_carbons, locant_map, list(principal_atoms))
    if subs:
        prefix = _build_prefix(subs)
        acid_name = f"{prefix}{acid_name}"

    # E/Z 立体化学 (Phase 316)
    stereo_prefix = ""
    if _ene_e:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_e = PrincipalChain(atom_indices=acid_carbons,
                               locant_map={c: i + 1 for i, c in enumerate(acid_carbons)})
        _stereo_e = assign_stereochemistry(graph, _pc_e)
        if _stereo_e:
            _combined_e = ",".join(d.strip("()") for d in _stereo_e)
            stereo_prefix = f"({_combined_e})-"

    return f"{stereo_prefix}{alkyl_name} {acid_name}"


def _name_acid_halide(graph, pgrp, get_atom) -> str:
    """
    酸ハライド命名: (stem)oyl (halide)
    例: CC(=O)Cl → ethanoyl chloride
        ClCC(=O)Cl → 2-chloroethanoyl chloride  (Phase 72)
    """
    from .constants import CHAIN_PREFIX, HALOGEN_PREFIX
    from .substituent import collect_substituents
    from .name_assembler import _build_prefix

    carbonyl_c = pgrp.atom_indices[0]

    # ハロゲン特定
    halide_name = "chloride"
    halogen_set: set[int] = set()
    for nb_idx in graph.adjacency[carbonyl_c]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol in HALOGEN_PREFIX:
            halide_name = {
                "F": "fluoride", "Cl": "chloride",
                "Br": "bromide", "I": "iodide",
            }.get(nb.symbol, "halide")
            halogen_set.add(nb_idx)
            break

    # 芳香環に直接結合したカルボニル → benzoyl 型 / ring-N-carbonyl 型 (Phase 107, 525)
    for nb_idx in graph.adjacency[carbonyl_c]:
        if nb_idx in halogen_set:
            continue
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "C" and nb.is_aromatic and nb.in_ring:
            ring_atoms = next(
                (set(rt) for rt in (graph.ring_atom_sets or []) if nb_idx in rt), set()
            )
            # 純粋ベンゼン環 → benzoyl
            if (len(ring_atoms) == 6
                    and all(get_atom(graph, a).symbol == "C" for a in ring_atoms)):
                return f"benzoyl {halide_name}"
            # ヘテロ芳香環 → ring-N-carbonyl halide
            has_het = any(get_atom(graph, a).symbol != "C" for a in ring_atoms)
            if has_het and all(get_atom(graph, a).is_aromatic for a in ring_atoms):
                aryl_pfx = _aryl_sulfonyl_prefix(graph, nb_idx, carbonyl_c, get_atom)
                if aryl_pfx is not None:
                    return f"{aryl_pfx}carbonyl {halide_name}"
            # 他の縮合芳香環等は後続の通常パスへ
            break

    # Phase 173: 非芳香族環に結合したカルボニル → cyclohexanecarbonyl chloride 型
    for nb_idx in graph.adjacency[carbonyl_c]:
        if nb_idx in halogen_set:
            continue
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "C" and nb.in_ring and not nb.is_aromatic:
            ring_atoms_set = next(
                (rt for rt in (graph.ring_atom_sets or []) if nb_idx in rt), None
            )
            if ring_atoms_set is not None and all(
                get_atom(graph, a).symbol == "C" for a in ring_atoms_set
            ):
                from .ring_handler import (
                    _assign_ring_locants,
                    collect_ring_substituents,
                    assemble_ring_name,
                )
                from .constants import CHAIN_PREFIX
                ring_list = list(ring_atoms_set)
                ring_chain = _assign_ring_locants(
                    graph, ring_list, False, "alkane", [carbonyl_c]
                )
                substituents = collect_ring_substituents(
                    graph, ring_chain, [carbonyl_c]
                )
                stem = CHAIN_PREFIX.get(len(ring_list), f"C{len(ring_list)}")
                ring_base = assemble_ring_name(ring_chain, substituents, "alkane", None, [])
                return f"{ring_base}carbonyl {halide_name}"

    # 酸側炭素鎖
    acid_carbons = _collect_acid_chain(graph, carbonyl_c, halogen_set, get_atom)
    n_acid = len(acid_carbons)
    stem = CHAIN_PREFIX.get(n_acid, f"C{n_acid}")
    _ene_h, _yne_h = _chain_multiple_bonds(graph, acid_carbons)
    if _ene_h or _yne_h:
        from .name_assembler import _format_multiple_bonds as _fmt_h
        base = f"{stem}{_fmt_h(_ene_h, _yne_h)}oyl {halide_name}"
    else:
        base = f"{stem}anoyl {halide_name}"

    # Phase 198: カルボニル C に 2 個のハロゲンが付く場合 → carbonyl dihalide
    if n_acid == 1:
        _halide_names = {"F": "fluoride", "Cl": "chloride", "Br": "bromide", "I": "iodide"}
        _all_hal = [
            nb for nb in graph.adjacency[carbonyl_c]
            if get_atom(graph, nb).symbol in _halide_names
        ]
        if len(_all_hal) == 2:
            hal_syms = sorted(
                _halide_names[get_atom(graph, h).symbol] for h in _all_hal
            )
            if hal_syms[0] == hal_syms[1]:
                return f"carbonyl di{hal_syms[0]}"
            return f"carbonyl {hal_syms[0]} {hal_syms[1]}"

    # 鎖上の置換基を収集 (Phase 72)
    locant_map = {c: i + 1 for i, c in enumerate(acid_carbons)}
    # 主官能基原子: carbonyl_c の =O, ハロゲン, carbonyl_c 自身
    principal_atoms = set(pgrp.atom_indices) | halogen_set
    for nb_idx in graph.adjacency[carbonyl_c]:
        if get_atom(graph, nb_idx).symbol == "O":
            principal_atoms.add(nb_idx)
    subs = collect_substituents(graph, acid_carbons, locant_map, list(principal_atoms))
    # Phase 392: for a 1-carbon chain the only locant is trivially 1 → omit it (IUPAC P-14.5.2)
    if n_acid == 1:
        subs = [(None, name) for _, name in subs]

    # E/Z 立体化学 (Phase 316)
    _stereo_prefix_ah = ""
    if _ene_h:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_ah = PrincipalChain(atom_indices=acid_carbons,
                                locant_map={c: i + 1 for i, c in enumerate(acid_carbons)})
        _stereo_ah = assign_stereochemistry(graph, _pc_ah)
        if _stereo_ah:
            _comb_ah = ",".join(d.strip("()") for d in _stereo_ah)
            _stereo_prefix_ah = f"({_comb_ah})-"

    if not subs:
        # Phase 120/254: 置換基なし 1C/2C の保留名 (IUPAC 2013 PIN)
        if n_acid == 1:
            return f"formyl {halide_name}"
        if n_acid == 2:
            return f"acetyl {halide_name}"
        return f"{_stereo_prefix_ah}{base}"
    prefix = _build_prefix(subs)
    return f"{_stereo_prefix_ah}{prefix}{base}"


def _name_diacid_halide(graph, pgrp, get_atom) -> str:
    """
    ジアシルハライド命名: {stem}anedioyl di{halide}
    例: ClC(=O)C(=O)Cl → ethanedioyl dichloride
        ClC(=O)CC(=O)Cl → propanedioyl dichloride
    """
    from .constants import CHAIN_PREFIX

    # atom_indices = [c1, o1, c2, o2] from aggregate_groups
    # Find the two carbonyl carbons (C with halogen neighbor)
    halide_map = {"F": "fluoride", "Cl": "chloride", "Br": "bromide", "I": "iodide"}
    halogen_syms = set(halide_map.keys())

    carbonyl_cs = []
    halide_names = []
    for ai in pgrp.atom_indices:
        atom = get_atom(graph, ai)
        if atom.symbol != "C":
            continue
        for nb_idx in graph.adjacency[ai]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol in halogen_syms:
                carbonyl_cs.append(ai)
                halide_names.append(halide_map[nb.symbol])
                break

    if len(carbonyl_cs) < 2:
        return "diacid halide"

    c1, c2 = carbonyl_cs[0], carbonyl_cs[1]
    chain_fwd = _collect_acid_chain(graph, c1, set(), get_atom)
    chain_rev = _collect_acid_chain(graph, c2, set(), get_atom)
    ene_fwd, yne_fwd = _chain_multiple_bonds(graph, chain_fwd)
    ene_rev, yne_rev = _chain_multiple_bonds(graph, chain_rev)
    mb_fwd = sorted(ene_fwd + yne_fwd)
    mb_rev = sorted(ene_rev + yne_rev)
    if mb_rev and (not mb_fwd or mb_rev < mb_fwd):
        acid_carbons = chain_rev
        _ene_dah, _yne_dah = ene_rev, yne_rev
    else:
        acid_carbons = chain_fwd
        _ene_dah, _yne_dah = ene_fwd, yne_fwd

    n = len(acid_carbons)
    stem = CHAIN_PREFIX.get(n, f"C{n}")
    if _ene_dah or _yne_dah:
        from .name_assembler import _format_multiple_bonds as _fmt_dah
        chain_part = f"{stem}{_fmt_dah(_ene_dah, _yne_dah)}edioyl"
    else:
        chain_part = f"{stem}anedioyl"

    stereo_pfx_dah = ""
    if _ene_dah:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_dah = PrincipalChain(atom_indices=acid_carbons,
                                 locant_map={c: i + 1 for i, c in enumerate(acid_carbons)})
        _stereo_dah = assign_stereochemistry(graph, _pc_dah)
        if _stereo_dah:
            stereo_pfx_dah = "(" + ",".join(d.strip("()") for d in _stereo_dah) + ")-"

    halide_names_sorted = sorted(halide_names)
    if halide_names_sorted[0] == halide_names_sorted[1]:
        halide_str = f"di{halide_names_sorted[0]}"
    else:
        halide_str = f"{halide_names_sorted[0]} {halide_names_sorted[1]}"

    return f"{stereo_pfx_dah}{chain_part} {halide_str}"


def _name_sulfonate_sulfinate_ester(graph, pgrp, get_atom) -> str:
    """
    スルホン酸エステル / スルフィン酸エステル命名 (Phase 82)
    例: CS(=O)OC → methyl methanesulfinate
        CS(=O)(=O)OC → methyl methanesulfonate
    """
    from .constants import CHAIN_PREFIX
    from .substituent import _name_carbon_substituent

    s_idx = pgrp.atom_indices[0]
    kind = "sulfonate" if pgrp.group_type == "sulfonate_ester" else "sulfinate"

    # S に結合した C (acid 側) と O-C (ester 側) を特定
    c_acid = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    o_ester = next((nb for nb in pgrp.atom_indices if get_atom(graph, nb).symbol == "O"
                    and any(get_atom(graph, occ).symbol == "C"
                            for occ in graph.adjacency[nb] if occ != s_idx)), None)

    # acid 側の鎖
    stereo_pfx_ssest = ""
    if c_acid:
        aryl = _aryl_sulfonyl_prefix(graph, c_acid[0], s_idx, get_atom)
        if aryl is not None:
            acid_name = f"{aryl}{kind}"
        else:
            acid_chain_se, locant_se = _chain_through_pivot(graph, c_acid[0], {s_idx}, get_atom)
            n = len(acid_chain_se)
            acid_stem = CHAIN_PREFIX.get(n, f"C{n}")
            _ene_ssest, _yne_ssest = _chain_multiple_bonds(graph, acid_chain_se)
            if _ene_ssest:
                from .stereochemistry import assign_stereochemistry
                from .chain_finder import PrincipalChain
                _pc_ssest = PrincipalChain(
                    atom_indices=acid_chain_se,
                    locant_map={c: i + 1 for i, c in enumerate(acid_chain_se)})
                _stereo_ssest = assign_stereochemistry(graph, _pc_ssest)
                if _stereo_ssest:
                    stereo_pfx_ssest = "(" + ",".join(d.strip("()") for d in _stereo_ssest) + ")-"
            if _ene_ssest or _yne_ssest:
                from .name_assembler import _format_multiple_bonds as _fmt_ssest
                mb_str_ssest = _fmt_ssest(_ene_ssest, _yne_ssest)
                if n >= 3:
                    acid_name = f"{acid_stem}{mb_str_ssest}e-{locant_se}-{kind}"
                else:
                    acid_name = f"{acid_stem}{mb_str_ssest}e{kind}"
            else:
                if n >= 3:
                    acid_name = f"{acid_stem}ane-{locant_se}-{kind}"
                else:
                    acid_name = f"{acid_stem}ane{kind}"
    else:
        acid_name = kind

    # ester 側の alkyl
    if o_ester is not None:
        ester_c = next((nb for nb in graph.adjacency[o_ester]
                        if nb != s_idx and get_atom(graph, nb).symbol == "C"), None)
        ester_alkyl = (_name_carbon_substituent(graph, ester_c, {o_ester})
                       if ester_c is not None else "methyl")
    else:
        ester_alkyl = "methyl"

    return f"{stereo_pfx_ssest}{ester_alkyl} {acid_name}"


def _name_sulfoxide_sulfone(graph, pgrp, get_atom) -> str:
    """
    スルホキシド・スルホン命名: IUPAC 2013 P-65.3.1 dialkyl sulfone/sulfoxide 形式
    例: CS(=O)C       → dimethyl sulfoxide
        CS(=O)(=O)C   → dimethyl sulfone
        CS(=O)CC      → ethyl methyl sulfoxide
        CS(=O)(=O)CC  → ethyl methyl sulfone
    """
    import re as _re_sox
    from .substituent import _name_carbon_substituent, name_substituent

    s_idx = pgrp.atom_indices[0]
    type_word = "sulfoxide" if pgrp.group_type == "sulfoxide" else "sulfone"

    c_neighbors = [nb for nb in graph.adjacency[s_idx]
                   if get_atom(graph, nb).symbol == "C"]
    if len(c_neighbors) < 2:
        return type_word

    c1, c2 = c_neighbors[0], c_neighbors[1]

    def _group_name(c_idx: int) -> str:
        atom = get_atom(graph, c_idx)
        if atom.in_ring and atom.is_aromatic:
            return name_substituent(graph, c_idx, {s_idx}) or "phenyl"
        return _name_carbon_substituent(graph, c_idx, {s_idx})

    name1 = _group_name(c1)
    name2 = _group_name(c2)

    def _needs_parens(nm: str) -> bool:
        return bool(_re_sox.search(r"[0-9]", nm)) or nm.startswith("(")

    def _alpha_base(nm: str) -> str:
        s = _re_sox.sub(r"^\(", "", nm)
        s = _re_sox.sub(r"^(di|tri|tetra|bis|tris)", "", s)
        return s.lower()

    names = sorted([name1, name2], key=_alpha_base)

    if names[0] == names[1]:
        nm = names[0]
        prefix = f"bis({nm})" if _needs_parens(nm) else f"di{nm}"
    else:
        parts = [f"({nm})" if _needs_parens(nm) else nm for nm in names]
        prefix = " ".join(parts)

    return f"{prefix} {type_word}"


def _name_selenoxide_selenone(graph, pgrp, get_atom) -> str:
    """セレノキシド・セレノン: dialkyl selenoxide / selenone (IUPAC 2013, Phase 519)"""
    import re as _re_seo
    from .substituent import _name_carbon_substituent, name_substituent

    se_idx = pgrp.atom_indices[0]
    type_word = "selenoxide" if pgrp.group_type == "selenoxide" else "selenone"

    c_neighbors = [nb for nb in graph.adjacency[se_idx]
                   if get_atom(graph, nb).symbol == "C"]
    if len(c_neighbors) < 2:
        return type_word

    def _group_name(c_idx: int) -> str:
        atom = get_atom(graph, c_idx)
        if atom.in_ring and atom.is_aromatic:
            return name_substituent(graph, c_idx, {se_idx}) or "phenyl"
        return _name_carbon_substituent(graph, c_idx, {se_idx})

    name1 = _group_name(c_neighbors[0])
    name2 = _group_name(c_neighbors[1])

    def _needs_parens(nm: str) -> bool:
        return bool(_re_seo.search(r"[0-9]", nm)) or nm.startswith("(")

    def _alpha_base(nm: str) -> str:
        s = _re_seo.sub(r"^\(", "", nm)
        s = _re_seo.sub(r"^(di|tri|tetra|bis|tris)", "", s)
        return s.lower()

    names = sorted([name1, name2], key=_alpha_base)

    if names[0] == names[1]:
        nm = names[0]
        prefix = f"bis({nm})" if _needs_parens(nm) else f"di{nm}"
    else:
        parts = [f"({nm})" if _needs_parens(nm) else nm for nm in names]
        prefix = " ".join(parts)

    return f"{prefix} {type_word}"


def _name_sulfonamide(graph, pgrp, get_atom) -> str:
    """
    スルホンアミド命名: {stem}anesulfonamide
    例: CS(=O)(=O)N → methanesulfonamide
        CS(=O)(=O)NC → N-methylmethanesulfonamide
    """
    from .constants import CHAIN_PREFIX, MULTIPLIER
    from .substituent import _name_carbon_substituent
    from collections import Counter

    s_idx = pgrp.atom_indices[0]

    # S に付いた C 鎖（酸側）
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_s:
        return "sulfonamide"
    aryl = _aryl_sulfonyl_prefix(graph, c_on_s[0], s_idx, get_atom)
    stereo_pfx_snam = ""
    if aryl is not None:
        # N 上置換基は後続の N-subs ロジックで追記するため base だけ作る
        base = f"{aryl}sulfonamide"
    else:
        chain_sn, locant_sn = _chain_through_pivot(graph, c_on_s[0], {s_idx}, get_atom)
        stem = CHAIN_PREFIX.get(len(chain_sn), f"C{len(chain_sn)}")
        _ene_sn, _yne_sn = _chain_multiple_bonds(graph, chain_sn)
        if _ene_sn:
            from .stereochemistry import assign_stereochemistry
            from .chain_finder import PrincipalChain
            _pc_snam = PrincipalChain(atom_indices=chain_sn,
                                      locant_map={c: i + 1 for i, c in enumerate(chain_sn)})
            _stereo_snam = assign_stereochemistry(graph, _pc_snam)
            if _stereo_snam:
                stereo_pfx_snam = "(" + ",".join(d.strip("()") for d in _stereo_snam) + ")-"
        if _ene_sn or _yne_sn:
            from .name_assembler import _format_multiple_bonds as _fmt_sn
            mb_str = _fmt_sn(_ene_sn, _yne_sn)
            if len(chain_sn) >= 3:
                base = f"{stem}{mb_str}e-{locant_sn}-sulfonamide"
            else:
                base = f"{stem}{mb_str}esulfonamide"
        else:
            if len(chain_sn) >= 3:
                base = f"{stem}ane-{locant_sn}-sulfonamide"
            else:
                base = f"{stem}anesulfonamide"

    # N 上の置換基
    n_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "N"]
    if not n_on_s:
        return f"{stereo_pfx_snam}{base}"

    n_idx = n_on_s[0]
    n_atom = get_atom(graph, n_idx)

    # Phase 394: N が環内にある場合 → {N_loc}-({sulfonyl_stem}sulfonyl){ring_base}
    if n_atom.in_ring:
        ring_set = next(
            (rt for rt in (graph.ring_atom_sets or []) if n_idx in rt), None
        )
        if ring_set is not None:
            from .heterocycle_handler import (
                _find_best_start as _fbs_sn,
                _build_locant_map as _blm_sn,
                _is_aromatic_ring as _iar_sn,
                _canonical_sig as _csig_sn,
                _RETAINED_NAMES as _RN_sn,
                _match_hantzsch_widman as _mhw_sn,
            )
            ring_list = list(ring_set)
            rot = _fbs_sn(ring_list, graph)
            lmap = _blm_sn(rot)
            n_locant = lmap.get(n_idx)
            is_arom = _iar_sn(ring_list, graph)
            csig = _csig_sn(rot, graph)
            retained = _RN_sn.get((is_arom, csig))
            if retained is not None:
                rbase, _ = retained
            else:
                rbase = _mhw_sn(ring_list, graph)
            if rbase is not None:
                sulfonyl_stem = base.replace("sulfonamide", "sulfonyl")
                loc_str = f"{n_locant}-" if n_locant is not None else ""
                return f"{loc_str}({stereo_pfx_snam}{sulfonyl_stem}){rbase}"

    c_on_n = [nb for nb in graph.adjacency[n_idx]
               if get_atom(graph, nb).symbol == "C"]
    if not c_on_n:
        return f"{stereo_pfx_snam}{base}"

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in c_on_n]
    sub_counts = Counter(n_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if sub.startswith("(") else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub_str}")

    prefix = "-".join(prefix_parts)
    return f"{stereo_pfx_snam}{prefix}{base}"


def _name_sulfonyl_azide(graph, pgrp, get_atom) -> str:
    """スルホニルアジド命名 (Phase 355): R-S(=O)₂-N₃ → {stem}anesulfonyl azide"""
    from .constants import CHAIN_PREFIX
    s_idx = pgrp.atom_indices[0]
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_s:
        return "sulfonyl azide"
    aryl_sa = _aryl_sulfonyl_prefix(graph, c_on_s[0], s_idx, get_atom)
    if aryl_sa is not None:
        return f"{aryl_sa}sulfonyl azide"
    chain_sa, locant_sa = _chain_through_pivot(graph, c_on_s[0], {s_idx}, get_atom)
    stem = CHAIN_PREFIX.get(len(chain_sa), f"C{len(chain_sa)}")
    _ene_sa, _yne_sa = _chain_multiple_bonds(graph, chain_sa)
    if _ene_sa or _yne_sa:
        from .name_assembler import _format_multiple_bonds as _fmt_sa
        mb_str = _fmt_sa(_ene_sa, _yne_sa)
        if len(chain_sa) >= 3:
            return f"{stem}{mb_str}e-{locant_sa}-sulfonyl azide"
        return f"{stem}{mb_str}esulfonyl azide"
    if len(chain_sa) >= 3:
        return f"{stem}ane-{locant_sa}-sulfonyl azide"
    return f"{stem}anesulfonyl azide"


def _name_sulfonohydrazide(graph, pgrp, get_atom) -> str:
    """スルホノヒドラジド命名 (Phase 353): R-S(=O)₂-NHNH₂ → {stem}anesulfonohydrazide"""
    from .constants import CHAIN_PREFIX, MULTIPLIER
    from .substituent import _name_carbon_substituent
    from collections import Counter

    s_idx = pgrp.atom_indices[0]
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_s:
        return "sulfonohydrazide"
    aryl = _aryl_sulfonyl_prefix(graph, c_on_s[0], s_idx, get_atom)
    if aryl is not None:
        base = f"{aryl}sulfonohydrazide"
    else:
        chain_sh, locant_sh = _chain_through_pivot(graph, c_on_s[0], {s_idx}, get_atom)
        stem = CHAIN_PREFIX.get(len(chain_sh), f"C{len(chain_sh)}")
        _ene_sh, _yne_sh = _chain_multiple_bonds(graph, chain_sh)
        if _ene_sh or _yne_sh:
            from .name_assembler import _format_multiple_bonds as _fmt_sh
            mb_str = _fmt_sh(_ene_sh, _yne_sh)
            if len(chain_sh) >= 3:
                base = f"{stem}{mb_str}e-{locant_sh}-sulfonohydrazide"
            else:
                base = f"{stem}{mb_str}esulfonohydrazide"
        else:
            if len(chain_sh) >= 3:
                base = f"{stem}ane-{locant_sh}-sulfonohydrazide"
            else:
                base = f"{stem}anesulfonohydrazide"

    # N' substituents on the terminal N2
    n1_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "N"]
    if not n1_on_s:
        return base
    n1_idx = n1_on_s[0]
    from .functional_group import get_bond_order as _gbo_sh
    n2_nbrs = [nb for nb in graph.adjacency[n1_idx]
               if nb != s_idx and get_atom(graph, nb).symbol == "N"
               and _gbo_sh(graph, n1_idx, nb) == 1.0]
    if not n2_nbrs:
        return base
    n2_idx = n2_nbrs[0]
    c_on_n2 = [nb for nb in graph.adjacency[n2_idx]
                if nb != n1_idx and get_atom(graph, nb).symbol == "C"]
    if not c_on_n2:
        return base
    n2_subs = [_name_carbon_substituent(graph, c, {n2_idx}) for c in c_on_n2]
    sub_counts = Counter(n2_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if sub.startswith("(") else sub
        if cnt == 1:
            prefix_parts.append(f"N'-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N',N'-{mult}{sub_str}")
    prefix = "-".join(prefix_parts)
    return f"{prefix}{base}"


def _name_sulfamic_acid(graph, pgrp, get_atom) -> str:
    """スルファミン酸 (N-置換含む): (N-R)sulfamic acid (Phase 303)"""
    from .constants import MULTIPLIER
    from .substituent import _name_carbon_substituent
    from collections import Counter

    s_idx = pgrp.atom_indices[0]
    n_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "N"]
    if not n_on_s:
        return "sulfamic acid"

    n_idx = n_on_s[0]
    c_on_n = [nb for nb in graph.adjacency[n_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_n:
        return "sulfamic acid"

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in c_on_n]
    sub_counts = Counter(n_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if sub.startswith("(") else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub_str}")

    prefix = "-".join(prefix_parts)
    return f"{prefix}sulfamic acid"


def _name_n_substituted_sulfamide(graph, pgrp, get_atom) -> str:
    """N-置換スルファミド: N-alkylsulfamide, N,N'-dialkylsulfamide 等 (Phase 342)"""
    from .constants import MULTIPLIER
    from .substituent import _name_carbon_substituent
    from collections import Counter

    s_idx = pgrp.atom_indices[0]
    n_atoms = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "N"]
    if len(n_atoms) < 2:
        return "sulfamide"

    # Collect substituents from both N atoms with their N-position label
    # Following IUPAC convention: N and N' for the two nitrogens
    all_sub_pairs: list[tuple[str, str]] = []  # (N-label, sub_name)
    for i, n_idx in enumerate(n_atoms):
        c_subs = [nb for nb in graph.adjacency[n_idx] if get_atom(graph, nb).symbol == "C"]
        for c in c_subs:
            sub_name = _name_carbon_substituent(graph, c, {n_idx, s_idx})
            label = "N" if i == 0 else "N'"
            all_sub_pairs.append((label, sub_name))

    if not all_sub_pairs:
        return "sulfamide"

    # Assign N/N' labels in alphabetical order of substituent name
    sub_names_sorted = sorted(sub for _, sub in all_sub_pairs)
    sub_counts_sf = Counter(sub_names_sorted)

    # Re-label: alphabetically first unique substituent gets "N", second gets "N'"
    n_labels = ["N", "N'"]
    unique_subs = sorted(sub_counts_sf.keys())
    prefix_parts_sf: list[str] = []
    label_idx = 0
    for sub in unique_subs:
        cnt = sub_counts_sf[sub]
        sub_str = f"({sub})" if sub.startswith("(") else sub
        if cnt == 1:
            lbl = n_labels[label_idx] if label_idx < len(n_labels) else f"N{label_idx + 1}"
            prefix_parts_sf.append(f"{lbl}-{sub_str}")
            label_idx += 1
        else:
            labels_used = [n_labels[i] if i < len(n_labels) else f"N{i + 1}"
                           for i in range(label_idx, label_idx + cnt)]
            labels_str = ",".join(labels_used)
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts_sf.append(f"{labels_str}-{mult}{sub_str}")
            label_idx += cnt

    prefix_sf = "-".join(prefix_parts_sf)
    return f"{prefix_sf}sulfamide"


def _name_sulfinylhydrazide(graph, pgrp, get_atom) -> str:
    """スルフィニルヒドラジド命名 (Phase 354): R-S(=O)-NHNH₂ → {stem}anesulfinylhydrazide"""
    from .constants import CHAIN_PREFIX, MULTIPLIER
    from .substituent import _name_carbon_substituent
    from .functional_group import get_bond_order as _gbo_sfh
    from collections import Counter

    s_idx = pgrp.atom_indices[0]
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_s:
        return "sulfinylhydrazide"
    aryl_sfh = _aryl_sulfonyl_prefix(graph, c_on_s[0], s_idx, get_atom)
    if aryl_sfh is not None:
        base = f"{aryl_sfh}sulfinylhydrazide"
    else:
        chain_sfh, locant_sfh = _chain_through_pivot(graph, c_on_s[0], {s_idx}, get_atom)
        stem = CHAIN_PREFIX.get(len(chain_sfh), f"C{len(chain_sfh)}")
        _ene_sfh, _yne_sfh = _chain_multiple_bonds(graph, chain_sfh)
        if _ene_sfh or _yne_sfh:
            from .name_assembler import _format_multiple_bonds as _fmt_sfh
            mb_str = _fmt_sfh(_ene_sfh, _yne_sfh)
            if len(chain_sfh) >= 3:
                base = f"{stem}{mb_str}e-{locant_sfh}-sulfinylhydrazide"
            else:
                base = f"{stem}{mb_str}esulfinylhydrazide"
        else:
            if len(chain_sfh) >= 3:
                base = f"{stem}ane-{locant_sfh}-sulfinylhydrazide"
            else:
                base = f"{stem}anesulfinylhydrazide"

    n1_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "N"]
    if not n1_on_s:
        return base
    n1_idx = n1_on_s[0]
    n2_nbrs = [nb for nb in graph.adjacency[n1_idx]
               if nb != s_idx and get_atom(graph, nb).symbol == "N"
               and _gbo_sfh(graph, n1_idx, nb) == 1.0]
    if not n2_nbrs:
        return base
    n2_idx = n2_nbrs[0]
    c_on_n2 = [nb for nb in graph.adjacency[n2_idx]
                if nb != n1_idx and get_atom(graph, nb).symbol == "C"]
    if not c_on_n2:
        return base
    n2_subs = [_name_carbon_substituent(graph, c, {n2_idx}) for c in c_on_n2]
    sub_counts = Counter(n2_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if sub.startswith("(") else sub
        if cnt == 1:
            prefix_parts.append(f"N'-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N',N'-{mult}{sub_str}")
    return "-".join(prefix_parts) + base


def _name_sulfinamide(graph, pgrp, get_atom) -> str:
    """スルフィナミド: {stem}anesulfinamide  (Phase 200)
    例: CS(=O)N → methanesulfinamide
        CS(=O)NC → N-methylmethanesulfinamide
    """
    from .constants import CHAIN_PREFIX, MULTIPLIER
    from .substituent import _name_carbon_substituent
    from collections import Counter

    s_idx = pgrp.atom_indices[0]
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_s:
        return "sulfinamide"

    chain_sn, locant_sn = _chain_through_pivot(graph, c_on_s[0], {s_idx}, get_atom)
    stem = CHAIN_PREFIX.get(len(chain_sn), f"C{len(chain_sn)}")
    _ene_sn, _yne_sn = _chain_multiple_bonds(graph, chain_sn)
    stereo_pfx_sfin = ""
    if _ene_sn:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_sfin = PrincipalChain(atom_indices=chain_sn,
                                   locant_map={c: i + 1 for i, c in enumerate(chain_sn)})
        _stereo_sfin = assign_stereochemistry(graph, _pc_sfin)
        if _stereo_sfin:
            stereo_pfx_sfin = "(" + ",".join(d.strip("()") for d in _stereo_sfin) + ")-"
    if _ene_sn or _yne_sn:
        from .name_assembler import _format_multiple_bonds as _fmt_sn
        mb_str = _fmt_sn(_ene_sn, _yne_sn)
        if len(chain_sn) >= 3:
            base = f"{stem}{mb_str}e-{locant_sn}-sulfinamide"
        else:
            base = f"{stem}{mb_str}esulfinamide"
    else:
        if len(chain_sn) >= 3:
            base = f"{stem}ane-{locant_sn}-sulfinamide"
        else:
            base = f"{stem}anesulfinamide"

    n_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "N"]
    if not n_on_s:
        return f"{stereo_pfx_sfin}{base}"

    n_idx = n_on_s[0]
    c_on_n = [nb for nb in graph.adjacency[n_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_n:
        return f"{stereo_pfx_sfin}{base}"

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in c_on_n]
    sub_counts = Counter(n_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if sub.startswith("(") else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub_str}")

    prefix = "-".join(prefix_parts)
    return f"{stereo_pfx_sfin}{prefix}{base}"


def _aryl_sulfonyl_prefix(graph, c_start, s_idx, get_atom) -> str | None:
    """
    c_start が芳香環の C であれば環の接頭辞を返す。
    ベンゼン: "benzene" (置換基つきなら "4-methylbenzene" 等)
    ヘテロ芳香族: "pyridine-2-" / "thiophene-2-" 等 (末尾ハイフンつき)
    """
    atom = get_atom(graph, c_start)
    if not (atom.is_aromatic and atom.in_ring):
        return None
    ring_atoms = next(
        (rt for rt in (graph.ring_atom_sets or []) if c_start in rt), None
    )
    if ring_atoms is None:
        return None

    # Benzene: existing path
    if (len(ring_atoms) == 6
            and all(get_atom(graph, a).symbol == "C" for a in ring_atoms)):
        from .ring_handler import _assign_ring_locants, collect_ring_substituents, assemble_ring_name
        ring_chain = _assign_ring_locants(graph, list(ring_atoms), True, "alkane", [s_idx])
        substituents = collect_ring_substituents(graph, ring_chain, [s_idx])
        return assemble_ring_name(ring_chain, substituents, "alkane", None, [])

    # Heteroaromatic ring: use _RETAINED_NAMES
    has_heteroatom = any(get_atom(graph, a).symbol != "C" for a in ring_atoms)
    if has_heteroatom and all(get_atom(graph, a).is_aromatic for a in ring_atoms):
        from .heterocycle_handler import _find_best_start, _canonical_sig, _RETAINED_NAMES
        ring_list = sorted(ring_atoms)
        ring_set = set(ring_list)
        ordered: list[int] = [ring_list[0]]
        seen: set[int] = {ring_list[0]}
        while len(ordered) < len(ring_list):
            cur = ordered[-1]
            for nb in graph.adjacency[cur]:
                if nb in ring_set and nb not in seen:
                    ordered.append(nb)
                    seen.add(nb)
                    break
            else:
                break
        rotation = _find_best_start(ordered, graph)
        sig = _canonical_sig(rotation, graph)
        entry = _RETAINED_NAMES.get((True, sig))
        if entry is not None:
            het_name, has_ind_h = entry
            ind_h = "1H-" if has_ind_h else ""
            rev_rotation = [rotation[0]] + list(reversed(rotation[1:]))
            best_loc: int | None = None
            for rot in (rotation, rev_rotation):
                loc = {a: i + 1 for i, a in enumerate(rot)}.get(c_start)
                if loc is not None and (best_loc is None or loc < best_loc):
                    best_loc = loc
            if best_loc is not None:
                return f"{ind_h}{het_name}-{best_loc}-"

    return None


def _name_sulfonic_acid(graph, pgrp, get_atom) -> str:
    """
    スルホン酸命名: {stem}anesulfonic acid
    例: CS(=O)(=O)O → methanesulfonic acid
        CCCS(=O)(=O)O → propane-1-sulfonic acid
        CC(S(=O)(=O)O)C → propane-2-sulfonic acid
    """
    from .constants import CHAIN_PREFIX
    s_idx = pgrp.atom_indices[0]
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_s:
        return "sulfonic acid"
    aryl = _aryl_sulfonyl_prefix(graph, c_on_s[0], s_idx, get_atom)
    if aryl is not None:
        return f"{aryl}sulfonic acid"
    chain, locant = _chain_through_pivot(graph, c_on_s[0], {s_idx}, get_atom)
    stem = CHAIN_PREFIX.get(len(chain), f"C{len(chain)}")
    _ene_so, _yne_so = _chain_multiple_bonds(graph, chain)
    stereo_pfx_so = ""
    if _ene_so:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_so = PrincipalChain(atom_indices=chain,
                                locant_map={c: i + 1 for i, c in enumerate(chain)})
        _stereo_so = assign_stereochemistry(graph, _pc_so)
        if _stereo_so:
            stereo_pfx_so = "(" + ",".join(d.strip("()") for d in _stereo_so) + ")-"
    if _ene_so or _yne_so:
        from .name_assembler import _format_multiple_bonds as _fmt_so
        mb_str = _fmt_so(_ene_so, _yne_so)
        if len(chain) >= 3:
            return f"{stereo_pfx_so}{stem}{mb_str}e-{locant}-sulfonic acid"
        return f"{stereo_pfx_so}{stem}{mb_str}esulfonic acid"
    if len(chain) >= 3:
        return f"{stem}ane-{locant}-sulfonic acid"
    return f"{stem}anesulfonic acid"


def _name_sulfonate_anion(graph, pgrp, get_atom) -> str:
    """スルホン酸アニオン: {stem}anesulfonate (Phase 386)"""
    from .constants import CHAIN_PREFIX
    s_idx = pgrp.atom_indices[0]
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_s:
        return "sulfonate"
    aryl = _aryl_sulfonyl_prefix(graph, c_on_s[0], s_idx, get_atom)
    if aryl is not None:
        return f"{aryl}sulfonate"
    chain, locant = _chain_through_pivot(graph, c_on_s[0], {s_idx}, get_atom)
    stem = CHAIN_PREFIX.get(len(chain), f"C{len(chain)}")
    _ene_so, _yne_so = _chain_multiple_bonds(graph, chain)
    if _ene_so or _yne_so:
        from .name_assembler import _format_multiple_bonds as _fmt_so
        mb_str = _fmt_so(_ene_so, _yne_so)
        if len(chain) >= 3:
            return f"{stem}{mb_str}e-{locant}-sulfonate"
        return f"{stem}{mb_str}esulfonate"
    if len(chain) >= 3:
        return f"{stem}ane-{locant}-sulfonate"
    return f"{stem}anesulfonate"


def _name_disulfonic_acid(graph, pgrp, get_atom) -> str:
    """ジスルホン酸: ethane-1,2-disulfonic acid etc."""
    from .constants import CHAIN_PREFIX
    s_atoms = [ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "S"]
    c_pivots = []
    for s_idx in s_atoms:
        c_nbs = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
        if c_nbs:
            c_pivots.append((c_nbs[0], s_idx))
    if len(c_pivots) < 2:
        return "disulfonic acid"
    c1, s1 = c_pivots[0]
    chain, _ = _chain_through_pivot(graph, c1, {s1}, get_atom)
    lm_fwd = {c: i + 1 for i, c in enumerate(chain)}
    lm_rev = {c: len(chain) - i for i, c in enumerate(chain)}
    locs_fwd = sorted(lm_fwd.get(c, 0) for c, _ in c_pivots)
    locs_rev = sorted(lm_rev.get(c, 0) for c, _ in c_pivots)
    locs = locs_fwd if locs_fwd <= locs_rev else locs_rev
    stem = CHAIN_PREFIX.get(len(chain), f"C{len(chain)}")
    loc_str = ",".join(str(l) for l in locs)
    return f"{stem}ane-{loc_str}-disulfonic acid"


def _name_disulfonamide(graph, pgrp, get_atom) -> str:
    """ジスルホンアミド: ethane-1,2-disulfonamide etc."""
    from .constants import CHAIN_PREFIX
    s_atoms = [ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "S"]
    c_pivots = []
    for s_idx in s_atoms:
        c_nbs = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
        if c_nbs:
            c_pivots.append((c_nbs[0], s_idx))
    if len(c_pivots) < 2:
        return "disulfonamide"
    c1, s1 = c_pivots[0]
    chain, _ = _chain_through_pivot(graph, c1, {s1}, get_atom)
    lm_fwd = {c: i + 1 for i, c in enumerate(chain)}
    lm_rev = {c: len(chain) - i for i, c in enumerate(chain)}
    locs_fwd = sorted(lm_fwd.get(c, 0) for c, _ in c_pivots)
    locs_rev = sorted(lm_rev.get(c, 0) for c, _ in c_pivots)
    locs = locs_fwd if locs_fwd <= locs_rev else locs_rev
    stem = CHAIN_PREFIX.get(len(chain), f"C{len(chain)}")
    loc_str = ",".join(str(l) for l in locs)
    return f"{stem}ane-{loc_str}-disulfonamide"


def _name_diimine(graph, pgrp, get_atom) -> str:
    """ジイミン命名: alkane-X,Y-diimine (Phase 301)"""
    from .constants import CHAIN_PREFIX
    from .molecule_analyzer import get_bond_order
    n_atoms = [ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "N"]
    c_pivots = []
    for n_idx in n_atoms:
        for nb in graph.adjacency[n_idx]:
            if (get_atom(graph, nb).symbol == "C"
                    and get_bond_order(graph, n_idx, nb) == 2.0):
                c_pivots.append((nb, n_idx))
                break
    if len(c_pivots) < 2:
        return "diimine"
    c1, n1 = c_pivots[0]
    chain, _ = _chain_through_pivot(graph, c1, {n1}, get_atom)
    lm_fwd = {c: i + 1 for i, c in enumerate(chain)}
    lm_rev = {c: len(chain) - i for i, c in enumerate(chain)}
    locs_fwd = sorted(lm_fwd.get(c, 0) for c, _ in c_pivots)
    locs_rev = sorted(lm_rev.get(c, 0) for c, _ in c_pivots)
    locs = locs_fwd if locs_fwd <= locs_rev else locs_rev
    stem = CHAIN_PREFIX.get(len(chain), f"C{len(chain)}")
    loc_str = ",".join(str(l) for l in locs)
    return f"{stem}ane-{loc_str}-diimine"


def _name_dial(graph, pgrp, get_atom):
    """ベンゼン/シクロアルカン環上の2つのアルデヒド基 → {ring}-X,Y-dicarbaldehyde (Phase 302/318)"""
    from .molecule_analyzer import get_bond_order
    from .constants import CHAIN_PREFIX

    formyl_cs = []
    for ai in pgrp.atom_indices:
        if get_atom(graph, ai).symbol != "C":
            continue
        for nb_idx in graph.adjacency[ai]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "O" and get_bond_order(graph, ai, nb_idx) == 2.0:
                formyl_cs.append(ai)
                break

    if len(formyl_cs) < 2:
        return None

    c1, c2 = formyl_cs[0], formyl_cs[1]

    ring1, rc1, arom1 = _ring_for_exo_c(graph, c1, get_atom)
    ring2, rc2, arom2 = _ring_for_exo_c(graph, c2, get_atom)

    if ring1 is not None and ring2 is not None and ring1 == ring2:
        from .ring_handler import _assign_ring_locants
        ring_list = sorted(ring1)
        ring_chain = _assign_ring_locants(graph, ring_list, arom1, "arene" if arom1 else "alkane",
                                          [rc1, rc2])
        loc1 = ring_chain.locant_map.get(rc1, 1)
        loc2 = ring_chain.locant_map.get(rc2, 2)
        locs = sorted([loc1, loc2])
        n = len(ring1)
        if arom1 and n == 6:
            ring_base = "benzene"
        else:
            stem = CHAIN_PREFIX.get(n, f"C{n}")
            ring_base = f"cyclo{stem}ane"
        return f"{ring_base}-{locs[0]},{locs[1]}-dicarbaldehyde"

    return None


def _name_sulfinic_acid(graph, pgrp, get_atom) -> str:
    """
    スルフィン酸命名: {stem}anesulfinic acid (Phase 38)
    例: CS(=O)O → methanesulfinic acid
        CCCS(=O)O → propane-1-sulfinic acid
    """
    from .constants import CHAIN_PREFIX
    s_idx = pgrp.atom_indices[0]
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_s:
        return "sulfinic acid"
    aryl = _aryl_sulfonyl_prefix(graph, c_on_s[0], s_idx, get_atom)
    if aryl is not None:
        return f"{aryl}sulfinic acid"
    chain, locant = _chain_through_pivot(graph, c_on_s[0], {s_idx}, get_atom)
    stem = CHAIN_PREFIX.get(len(chain), f"C{len(chain)}")
    _ene_si, _yne_si = _chain_multiple_bonds(graph, chain)
    stereo_pfx_si = ""
    if _ene_si:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_si = PrincipalChain(atom_indices=chain,
                                locant_map={c: i + 1 for i, c in enumerate(chain)})
        _stereo_si = assign_stereochemistry(graph, _pc_si)
        if _stereo_si:
            stereo_pfx_si = "(" + ",".join(d.strip("()") for d in _stereo_si) + ")-"
    if _ene_si or _yne_si:
        from .name_assembler import _format_multiple_bonds as _fmt_si
        mb_str = _fmt_si(_ene_si, _yne_si)
        if len(chain) >= 3:
            return f"{stereo_pfx_si}{stem}{mb_str}e-{locant}-sulfinic acid"
        return f"{stereo_pfx_si}{stem}{mb_str}esulfinic acid"
    if len(chain) >= 3:
        return f"{stem}ane-{locant}-sulfinic acid"
    return f"{stem}anesulfinic acid"


def _name_sulfenic_acid(graph, pgrp, get_atom) -> str:
    """
    スルフェン酸命名: {stem}anesulfenic acid (Phase 166)
    例: CSO → methanesulfenic acid
        CCSO → ethanesulfenic acid
    """
    from .constants import CHAIN_PREFIX
    s_idx = pgrp.atom_indices[0]
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_s:
        return "sulfenic acid"
    aryl = _aryl_sulfonyl_prefix(graph, c_on_s[0], s_idx, get_atom)
    if aryl is not None:
        return f"{aryl}sulfenic acid"
    chain, locant = _chain_through_pivot(graph, c_on_s[0], {s_idx}, get_atom)
    stem = CHAIN_PREFIX.get(len(chain), f"C{len(chain)}")
    _ene_se, _yne_se = _chain_multiple_bonds(graph, chain)
    stereo_pfx_se = ""
    if _ene_se:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_se = PrincipalChain(atom_indices=chain,
                                locant_map={c: i + 1 for i, c in enumerate(chain)})
        _stereo_se = assign_stereochemistry(graph, _pc_se)
        if _stereo_se:
            stereo_pfx_se = "(" + ",".join(d.strip("()") for d in _stereo_se) + ")-"
    if _ene_se or _yne_se:
        from .name_assembler import _format_multiple_bonds as _fmt_se
        mb_str = _fmt_se(_ene_se, _yne_se)
        if len(chain) >= 3:
            return f"{stereo_pfx_se}{stem}{mb_str}e-{locant}-sulfenic acid"
        return f"{stereo_pfx_se}{stem}{mb_str}esulfenic acid"
    if len(chain) >= 3:
        return f"{stem}ane-{locant}-sulfenic acid"
    return f"{stem}anesulfenic acid"


def _name_sulfenyl_halide(graph, pgrp, get_atom) -> str:
    """スルフェニルハライド命名: {stem}anesulfenyl {halide} (Phase 375)"""
    from .constants import CHAIN_PREFIX
    _HALIDE_NAMES = {"Cl": "chloride", "F": "fluoride", "Br": "bromide", "I": "iodide"}
    s_idx = pgrp.atom_indices[0]
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    hal = next((get_atom(graph, nb).symbol for nb in graph.adjacency[s_idx]
                if get_atom(graph, nb).symbol in _HALIDE_NAMES), "Cl")
    hal_name = _HALIDE_NAMES.get(hal, "chloride")
    if not c_on_s:
        return f"sulfenyl {hal_name}"
    aryl = _aryl_sulfonyl_prefix(graph, c_on_s[0], s_idx, get_atom)
    if aryl is not None:
        return f"{aryl}sulfenyl {hal_name}"
    chain, locant = _chain_through_pivot(graph, c_on_s[0], {s_idx}, get_atom)
    stem = CHAIN_PREFIX.get(len(chain), f"C{len(chain)}")
    _ene_sfh, _yne_sfh = _chain_multiple_bonds(graph, chain)
    if _ene_sfh or _yne_sfh:
        from .name_assembler import _format_multiple_bonds as _fmt_sfh
        mb_str = _fmt_sfh(_ene_sfh, _yne_sfh)
        if len(chain) >= 3:
            return f"{stem}{mb_str}e-{locant}-sulfenyl {hal_name}"
        return f"{stem}{mb_str}esulfenyl {hal_name}"
    if len(chain) >= 3:
        return f"{stem}ane-{locant}-sulfenyl {hal_name}"
    return f"{stem}anesulfenyl {hal_name}"


def _name_sulfenate_ester(graph, pgrp, get_atom) -> str:
    """スルフェン酸エステル命名: {alkyl} {stem}anesulfenate (Phase 375)"""
    from .constants import CHAIN_PREFIX
    from .substituent import _name_carbon_substituent
    s_idx = pgrp.atom_indices[0]
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    o_idx = next((nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "O"), None)
    if o_idx is None or not c_on_s:
        return "sulfenate"
    # alkyl (O 側)
    ester_c = next(
        (nb for nb in graph.adjacency[o_idx]
         if nb != s_idx and get_atom(graph, nb).symbol == "C"),
        None,
    )
    alkyl = _name_carbon_substituent(graph, ester_c, {o_idx}) if ester_c else "methyl"
    aryl = _aryl_sulfonyl_prefix(graph, c_on_s[0], s_idx, get_atom)
    if aryl is not None:
        return f"{alkyl} {aryl}sulfenate"
    chain, locant = _chain_through_pivot(graph, c_on_s[0], {s_idx}, get_atom)
    stem = CHAIN_PREFIX.get(len(chain), f"C{len(chain)}")
    _ene_sfe, _yne_sfe = _chain_multiple_bonds(graph, chain)
    if _ene_sfe or _yne_sfe:
        from .name_assembler import _format_multiple_bonds as _fmt_sfe
        mb_str = _fmt_sfe(_ene_sfe, _yne_sfe)
        if len(chain) >= 3:
            return f"{alkyl} {stem}{mb_str}e-{locant}-sulfenate"
        return f"{alkyl} {stem}{mb_str}esulfenate"
    if len(chain) >= 3:
        return f"{alkyl} {stem}ane-{locant}-sulfenate"
    return f"{alkyl} {stem}anesulfenate"


def _name_sulfonyl_chloride(graph, pgrp, get_atom) -> str:
    """スルホニルハライド命名: {stem}anesulfonyl {halide} (Phase 59/110/177)"""
    from .constants import CHAIN_PREFIX
    _HALIDE_NAMES = {"Cl": "chloride", "F": "fluoride", "Br": "bromide", "I": "iodide"}
    s_idx = pgrp.atom_indices[0]
    # ハロゲン特定
    halide_name = "chloride"
    for nb in graph.adjacency[s_idx]:
        sym = get_atom(graph, nb).symbol
        if sym in _HALIDE_NAMES:
            halide_name = _HALIDE_NAMES[sym]
            break
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_s:
        return f"sulfonyl {halide_name}"
    aryl = _aryl_sulfonyl_prefix(graph, c_on_s[0], s_idx, get_atom)
    if aryl is not None:
        return f"{aryl}sulfonyl {halide_name}"
    chain_sc, locant_sc = _chain_through_pivot(graph, c_on_s[0], {s_idx}, get_atom)
    stem = CHAIN_PREFIX.get(len(chain_sc), f"C{len(chain_sc)}")
    _ene_sc, _yne_sc = _chain_multiple_bonds(graph, chain_sc)
    stereo_pfx_sc = ""
    if _ene_sc:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_sc = PrincipalChain(atom_indices=chain_sc,
                                locant_map={c: i + 1 for i, c in enumerate(chain_sc)})
        _stereo_sc = assign_stereochemistry(graph, _pc_sc)
        if _stereo_sc:
            stereo_pfx_sc = "(" + ",".join(d.strip("()") for d in _stereo_sc) + ")-"
    if _ene_sc or _yne_sc:
        from .name_assembler import _format_multiple_bonds as _fmt_sc
        mb_str = _fmt_sc(_ene_sc, _yne_sc)
        if len(chain_sc) >= 3:
            return f"{stereo_pfx_sc}{stem}{mb_str}e-{locant_sc}-sulfonyl {halide_name}"
        return f"{stereo_pfx_sc}{stem}{mb_str}esulfonyl {halide_name}"
    if len(chain_sc) >= 3:
        return f"{stem}ane-{locant_sc}-sulfonyl {halide_name}"
    return f"{stem}anesulfonyl {halide_name}"


def _name_chloroformate(graph, pgrp, get_atom) -> str:
    """ハロホルメート命名: {alkyl} carbono{halo}idate (Phase 60/65)"""
    from .substituent import _name_carbon_substituent
    _HALO_SUFFIX = {"F": "fluoridic acid", "Cl": "chloridate", "Br": "bromate", "I": "iodate"}
    _HALO_NAME = {"F": "carbonofluoridate", "Cl": "carbonochloridate",
                  "Br": "carbonobromate", "I": "carbonoiodate"}
    carbonyl_c = pgrp.atom_indices[0]
    # ハロゲン特定
    halide_suffix = "carbonochloridate"
    for nb in graph.adjacency[carbonyl_c]:
        sym = get_atom(graph, nb).symbol
        if sym in _HALO_NAME:
            halide_suffix = _HALO_NAME[sym]
            break
    ester_o_idx = next(
        (ai for ai in pgrp.atom_indices[2:]
         if get_atom(graph, ai).symbol == "O"),
        None,
    )
    if ester_o_idx is None:
        return halide_suffix
    alkyl_c = next(
        (nb for nb in graph.adjacency[ester_o_idx]
         if nb != carbonyl_c and get_atom(graph, nb).symbol == "C"),
        None,
    )
    if alkyl_c is None:
        return halide_suffix
    alkyl_name = _name_carbon_substituent(graph, alkyl_c, {ester_o_idx})
    return f"{alkyl_name} {halide_suffix}"


def _name_sulfide(graph, pgrp, get_atom) -> str:
    """
    チオエーテル命名: IUPAC 2013 P-63.6.1 dialkyl sulfide 形式
    例: CSC → dimethyl sulfide
        CSCC → ethyl methyl sulfide
        CSc1ccccc1 → methyl phenyl sulfide
    """
    import re as _re_sf
    from .substituent import _name_carbon_substituent, name_substituent

    s_idx = pgrp.atom_indices[0]
    c_neighbors = [nb for nb in graph.adjacency[s_idx]
                   if get_atom(graph, nb).symbol == "C"]
    if len(c_neighbors) < 2:
        return "sulfide"

    c1, c2 = c_neighbors[0], c_neighbors[1]

    def _group_name(c_idx: int) -> str:
        atom = get_atom(graph, c_idx)
        if atom.in_ring and atom.is_aromatic:
            return name_substituent(graph, c_idx, {s_idx}) or "phenyl"
        return _name_carbon_substituent(graph, c_idx, {s_idx})

    name1 = _group_name(c1)
    name2 = _group_name(c2)

    def _needs_parens(nm: str) -> bool:
        return bool(_re_sf.search(r"[0-9]", nm)) or nm.startswith("(")

    def _alpha_base(nm: str) -> str:
        s = _re_sf.sub(r"^\(", "", nm)
        s = _re_sf.sub(r"^(di|tri|tetra|bis|tris)", "", s)
        return s.lower()

    names = sorted([name1, name2], key=_alpha_base)

    if names[0] == names[1]:
        nm = names[0]
        prefix = f"bis({nm})" if _needs_parens(nm) else f"di{nm}"
    else:
        parts = [f"({nm})" if _needs_parens(nm) else nm for nm in names]
        prefix = " ".join(parts)

    return f"{prefix} sulfide"


def _name_sulfinyl_chloride(graph, pgrp, get_atom) -> str:
    """スルフィニルハライド命名: {stem}anesulfinyl {halide} (Phase 224)"""
    from .constants import CHAIN_PREFIX
    _HALIDE_NAMES = {"Cl": "chloride", "F": "fluoride", "Br": "bromide", "I": "iodide"}
    s_idx = pgrp.atom_indices[0]
    halide_name = "chloride"
    for nb in graph.adjacency[s_idx]:
        sym = get_atom(graph, nb).symbol
        if sym in _HALIDE_NAMES:
            halide_name = _HALIDE_NAMES[sym]
            break
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_s:
        return f"sulfinyl {halide_name}"
    aryl = _aryl_sulfonyl_prefix(graph, c_on_s[0], s_idx, get_atom)
    if aryl is not None:
        return f"{aryl}sulfinyl {halide_name}"
    chain_sc, locant_sc = _chain_through_pivot(graph, c_on_s[0], {s_idx}, get_atom)
    stem = CHAIN_PREFIX.get(len(chain_sc), f"C{len(chain_sc)}")
    _ene_sc, _yne_sc = _chain_multiple_bonds(graph, chain_sc)
    if _ene_sc or _yne_sc:
        from .name_assembler import _format_multiple_bonds as _fmt_sc
        mb_str = _fmt_sc(_ene_sc, _yne_sc)
        if len(chain_sc) >= 3:
            return f"{stem}{mb_str}e-{locant_sc}-sulfinyl {halide_name}"
        return f"{stem}{mb_str}esulfinyl {halide_name}"
    if len(chain_sc) >= 3:
        return f"{stem}ane-{locant_sc}-sulfinyl {halide_name}"
    return f"{stem}anesulfinyl {halide_name}"


def _name_carbonate(graph, pgrp, get_atom) -> str:
    """
    炭酸エステル命名: {alkyl1} {alkyl2} carbonate / {alkyl} hydrogen carbonate (Phase 39/250)
    例: COC(=O)OC → dimethyl carbonate
        COC(=O)O  → methyl hydrogen carbonate
    """
    from .constants import MULTIPLIER
    from .substituent import _name_carbon_substituent
    from collections import Counter

    carbonyl_c = pgrp.atom_indices[0]
    ester_os = [nb for nb in graph.adjacency[carbonyl_c]
                if get_atom(graph, nb).symbol == "O"
                and any(get_atom(graph, on).symbol == "C"
                        for on in graph.adjacency[nb] if on != carbonyl_c)]
    oh_os = [nb for nb in graph.adjacency[carbonyl_c]
             if get_atom(graph, nb).symbol == "O"
             and any(get_atom(graph, on).symbol == "H" for on in graph.adjacency[nb])]
    alkyl_names = []
    for o_idx in ester_os:
        c_nbrs = [nb for nb in graph.adjacency[o_idx]
                  if nb != carbonyl_c and get_atom(graph, nb).symbol == "C"]
        if c_nbrs:
            alkyl_names.append(_name_carbon_substituent(graph, c_nbrs[0], {o_idx, carbonyl_c}))

    if not alkyl_names:
        return "carbonate"

    alkyl_names.sort()
    counts = Counter(alkyl_names)
    parts = []
    for name in sorted(counts):
        cnt = counts[name]
        if cnt == 1:
            parts.append(name)
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            parts.append(f"{mult}{name}")
    # half-ester: one free OH → "R hydrogen carbonate"
    if oh_os and len(ester_os) == 1:
        parts.append("hydrogen")
    return " ".join(parts) + " carbonate"


def _name_carbonothioate(graph, pgrp, get_atom) -> str:
    """
    チオ炭酸エステル命名 (Phase 348): dialkyl carbonothioate (RO)2C=S
    例: COC(=S)OC → dimethyl carbonothioate
    """
    from .constants import MULTIPLIER
    from .substituent import _name_carbon_substituent
    from .functional_group import get_bond_order
    from collections import Counter

    central_c = pgrp.atom_indices[0]
    o_idxs = [ai for ai in pgrp.atom_indices[1:]
              if get_atom(graph, ai).symbol == "O"]
    alkyl_names = []
    for o_idx in o_idxs:
        c_nbrs = [nb for nb in graph.adjacency[o_idx]
                  if nb != central_c and get_atom(graph, nb).symbol == "C"]
        if c_nbrs:
            alkyl_names.append(_name_carbon_substituent(graph, c_nbrs[0], {o_idx, central_c}))
    if not alkyl_names:
        return "carbonothioate"
    alkyl_names.sort()
    counts = Counter(alkyl_names)
    parts = []
    for name in sorted(counts):
        cnt = counts[name]
        if cnt == 1:
            parts.append(name)
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            parts.append(f"{mult}{name}")
    return " ".join(parts) + " carbonothioate"


def _name_carbonodithioate(graph, pgrp, get_atom) -> str:
    """
    ジチオ炭酸エステル命名 (Phase 348): dialkyl carbonodithioate (RS)2C=S
    例: CSC(=S)SC → dimethyl carbonodithioate
    """
    from .constants import MULTIPLIER
    from .substituent import _name_carbon_substituent
    from .functional_group import get_bond_order
    from collections import Counter

    central_c = pgrp.atom_indices[0]
    s_idxs = [ai for ai in pgrp.atom_indices[1:]
              if get_atom(graph, ai).symbol == "S"
              and any(get_atom(graph, nb).symbol == "C"
                      for nb in graph.adjacency[ai] if nb != central_c)]
    alkyl_names = []
    for s_idx in s_idxs:
        c_nbrs = [nb for nb in graph.adjacency[s_idx]
                  if nb != central_c and get_atom(graph, nb).symbol == "C"]
        if c_nbrs:
            alkyl_names.append(_name_carbon_substituent(graph, c_nbrs[0], {s_idx, central_c}))
    if not alkyl_names:
        return "carbonodithioate"
    alkyl_names.sort()
    counts = Counter(alkyl_names)
    parts = []
    for name in sorted(counts):
        cnt = counts[name]
        if cnt == 1:
            parts.append(name)
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            parts.append(f"{mult}{name}")
    return " ".join(parts) + " carbonodithioate"


def _name_phosphate_ester(graph, pgrp, get_atom) -> str:
    """ホスフェートエステル: triR phosphate / diR hydrogen phosphate  (Phase 145)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from .molecule_analyzer import get_bond_order
    from collections import Counter
    p_idx = pgrp.atom_indices[0]
    # O-C 結合の O インデックスを収集
    o_ester_idxs = [nb for nb in graph.adjacency[p_idx]
                    if get_atom(graph, nb).symbol == "O"
                    and get_bond_order(graph, p_idx, nb) == 1.0]
    alkyl_names = []
    for o_idx in o_ester_idxs:
        c_neighbors = [nb for nb in graph.adjacency[o_idx]
                       if nb != p_idx and get_atom(graph, nb).symbol == "C"]
        if c_neighbors:
            alkyl_names.append(_name_carbon_substituent(graph, c_neighbors[0], {o_idx}))
        else:
            alkyl_names.append("hydrogen")
    alkyl_names.sort()
    counts = Counter(alkyl_names)
    # 有機置換基をアルファベット順に並べ、hydrogen を最後に
    organic_subs = sorted(s for s in counts if s != "hydrogen")
    parts = []
    for sub in organic_subs:
        n = counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    if "hydrogen" in counts:
        n = counts["hydrogen"]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}hydrogen")
    return " ".join(parts) + " phosphate"


def _name_phosphonate_halfester(graph, pgrp, get_atom) -> str:
    """{ester_alkyl} hydrogen {p_alkyl}phosphonate  (Phase 253)"""
    from .substituent import _name_carbon_substituent
    from .molecule_analyzer import get_bond_order
    p_idx = pgrp.atom_indices[0]
    pc_neighbors = [nb for nb in graph.adjacency[p_idx] if get_atom(graph, nb).symbol == "C"]
    o_single = [nb for nb in graph.adjacency[p_idx]
                if get_atom(graph, nb).symbol == "O" and get_bond_order(graph, p_idx, nb) == 1.0]
    ester_names = []
    for o_idx in o_single:
        c_nbrs = [nb for nb in graph.adjacency[o_idx]
                  if nb != p_idx and get_atom(graph, nb).symbol == "C"]
        if c_nbrs:
            ester_names.append(_name_carbon_substituent(graph, c_nbrs[0], {o_idx}))
    ester_names.sort()
    p_alkyl = _name_carbon_substituent(graph, pc_neighbors[0], {p_idx}) if pc_neighbors else ""
    return " ".join(ester_names) + " hydrogen " + p_alkyl + "phosphonate"


def _name_phosphonate_ester(graph, pgrp, get_atom) -> str:
    """ホスホネートエステル: diR' Rphosphonate  (Phase 145)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from .molecule_analyzer import get_bond_order
    from collections import Counter
    p_idx = pgrp.atom_indices[0]
    # P-C 直接結合 (P-alkyl)
    pc_neighbors = [nb for nb in graph.adjacency[p_idx] if get_atom(graph, nb).symbol == "C"]
    # O-アルキル結合
    o_ester_idxs = [nb for nb in graph.adjacency[p_idx]
                    if get_atom(graph, nb).symbol == "O"
                    and get_bond_order(graph, p_idx, nb) == 1.0]
    ester_names = []
    for o_idx in o_ester_idxs:
        c_neighbors = [nb for nb in graph.adjacency[o_idx]
                       if nb != p_idx and get_atom(graph, nb).symbol == "C"]
        if c_neighbors:
            ester_names.append(_name_carbon_substituent(graph, c_neighbors[0], {o_idx}))
        else:
            ester_names.append("hydrogen")
    ester_names.sort()
    counts_ester = Counter(ester_names)
    parts = []
    for sub in sorted(counts_ester):
        n = counts_ester[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    ester_part = " ".join(parts)
    # P 直接結合アルキル
    if pc_neighbors:
        p_alkyl = _name_carbon_substituent(graph, pc_neighbors[0], {p_idx})
    else:
        p_alkyl = "phosphonate"
    return f"{ester_part} {p_alkyl}phosphonate"


def _name_phosphinate_ester(graph, pgrp, get_atom) -> str:
    """ホスフィネートエステル: R' diRphosphinate  (Phase 145)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from .molecule_analyzer import get_bond_order
    from collections import Counter
    p_idx = pgrp.atom_indices[0]
    pc_neighbors = [nb for nb in graph.adjacency[p_idx] if get_atom(graph, nb).symbol == "C"]
    o_ester_idxs = [nb for nb in graph.adjacency[p_idx]
                    if get_atom(graph, nb).symbol == "O"
                    and get_bond_order(graph, p_idx, nb) == 1.0]
    ester_names = []
    for o_idx in o_ester_idxs:
        c_neighbors = [nb for nb in graph.adjacency[o_idx]
                       if nb != p_idx and get_atom(graph, nb).symbol == "C"]
        if c_neighbors:
            ester_names.append(_name_carbon_substituent(graph, c_neighbors[0], {o_idx}))
    if ester_names:
        ester_part = ester_names[0]
    else:
        ester_part = "hydrogen"
    pc_names = sorted(_name_carbon_substituent(graph, c, {p_idx}) for c in pc_neighbors)
    counts_pc = Counter(pc_names)
    pc_parts = []
    for sub in sorted(counts_pc):
        n = counts_pc[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        pc_parts.append(f"{mult}{sub}")
    p_part = "".join(pc_parts)
    return f"{ester_part} {p_part}phosphinate"


def _name_phosphonothioate_ester(graph, pgrp, get_atom) -> str:
    """ホスホノチオアートエステル: diR' Rphosphonothioate (Phase 372)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from .molecule_analyzer import get_bond_order
    from collections import Counter
    p_idx = pgrp.atom_indices[0]
    pc_neighbors = [nb for nb in graph.adjacency[p_idx] if get_atom(graph, nb).symbol == "C"]
    o_ester_idxs = [nb for nb in graph.adjacency[p_idx]
                    if get_atom(graph, nb).symbol == "O"
                    and get_bond_order(graph, p_idx, nb) == 1.0]
    ester_names = []
    for o_idx in o_ester_idxs:
        c_neighbors = [nb for nb in graph.adjacency[o_idx]
                       if nb != p_idx and get_atom(graph, nb).symbol == "C"]
        if c_neighbors:
            ester_names.append(_name_carbon_substituent(graph, c_neighbors[0], {o_idx}))
        else:
            ester_names.append("hydrogen")
    ester_names.sort()
    counts_ester = Counter(ester_names)
    parts = []
    for sub in sorted(counts_ester):
        n = counts_ester[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    ester_part = " ".join(parts)
    if pc_neighbors:
        p_alkyl = _name_carbon_substituent(graph, pc_neighbors[0], {p_idx})
    else:
        p_alkyl = ""
    return f"{ester_part} {p_alkyl}phosphonothioate"


def _name_phosphonic_acid(graph, pgrp, get_atom) -> str:
    """ホスホン酸: {alkyl}phosphonic acid  (Phase 143)"""
    from .substituent import _name_carbon_substituent
    p_idx = pgrp.atom_indices[0]
    c_neighbors = [nb for nb in graph.adjacency[p_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_neighbors:
        return "phosphonic acid"
    sub = _name_carbon_substituent(graph, c_neighbors[0], {p_idx})
    return f"{sub}phosphonic acid"


def _name_by_c_substituents(graph, pgrp, get_atom, suffix: str) -> str:
    """共通ヘルパー: 中心原子の C 隣接から置換基名を組み立て + suffix を返す。"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter
    central = pgrp.atom_indices[0]
    c_neighbors = [nb for nb in graph.adjacency[central] if get_atom(graph, nb).symbol == "C"]
    if not c_neighbors:
        return suffix
    names = sorted(_name_carbon_substituent(graph, c, {central}) for c in c_neighbors)
    counts = Counter(names)
    parts = []
    for sub in sorted(counts):
        n = counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    return "".join(parts) + suffix


def _name_phosphinic_acid(graph, pgrp, get_atom) -> str:
    """ホスフィン酸: di{alkyl}phosphinic acid  (Phase 143)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "phosphinic acid")


def _name_phosphonous_acid(graph, pgrp, get_atom) -> str:
    """ホスホナス酸: {alkyl}phosphonous acid  (Phase 241)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "phosphonous acid")


def _name_arsonic_acid(graph, pgrp, get_atom) -> str:
    """ヒ素酸: {alkyl}arsonic acid  (Phase 242)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "arsonic acid")


def _name_arsinic_acid(graph, pgrp, get_atom) -> str:
    """アルシン酸: di{alkyl}arsinic acid  (Phase 242)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "arsinic acid")


def _name_arsonous_acid(graph, pgrp, get_atom) -> str:
    """亜ヒ酸: {alkyl}arsonous acid  (Phase 242)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "arsonous acid")


def _name_arsinous_acid(graph, pgrp, get_atom) -> str:
    """亜アルシン酸: di{alkyl}arsinous acid  (Phase 242)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "arsinous acid")


def _name_phosphinous_acid(graph, pgrp, get_atom) -> str:
    """ホスフィナス酸: {alkyl}phosphinous acid  (Phase 200)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "phosphinous acid")


def _name_phosphane(graph, pgrp, get_atom) -> str:
    """ホスファン: {alkyl(s)}phosphane  (Phase 143)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "phosphane")


def _name_phosphine_oxide(graph, pgrp, get_atom) -> str:
    """ホスファンオキシド: {alkyl(s)}phosphane oxide (IUPAC 2013 P-68.3.2)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "phosphane oxide")


def _name_phosphite_ester(graph, pgrp, get_atom) -> str:
    """亜リン酸エステル: {alkyl(s)} phosphite (Phase 187)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter
    p_idx = pgrp.atom_indices[0]
    o_neighbors = [nb for nb in graph.adjacency[p_idx] if get_atom(graph, nb).symbol == "O"]
    alkyl_names: list[str] = []
    for o_idx in o_neighbors:
        c_on_o = [nb for nb in graph.adjacency[o_idx] if nb != p_idx
                  and get_atom(graph, nb).symbol == "C"]
        if c_on_o:
            alkyl_names.append(_name_carbon_substituent(graph, c_on_o[0], {o_idx}))
    if not alkyl_names:
        return "phosphite"
    counts = Counter(alkyl_names)
    parts = []
    for sub in sorted(counts):
        n = counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    return " ".join(parts) + " phosphite"


def _name_borate_ester(graph, pgrp, get_atom) -> str:
    """トリアルコキシボラン: {alkoxy}borane (Phase 227, IUPAC 2013 P-68.1.3)"""
    from .substituent import _name_carbon_substituent, _make_oxy_name
    from .constants import MULTIPLIER
    from collections import Counter
    b_idx = pgrp.atom_indices[0]
    o_neighbors = [nb for nb in graph.adjacency[b_idx] if get_atom(graph, nb).symbol == "O"]
    alkoxy_names: list[str] = []
    for o_idx in o_neighbors:
        c_on_o = [nb for nb in graph.adjacency[o_idx] if nb != b_idx
                  and get_atom(graph, nb).symbol == "C"]
        if c_on_o:
            alkyl = _name_carbon_substituent(graph, c_on_o[0], {o_idx})
            alkoxy_names.append(_make_oxy_name(alkyl))
    if not alkoxy_names:
        return "borane"
    counts = Counter(alkoxy_names)
    all_same = len(counts) == 1
    parts = []
    for sub in sorted(counts):
        n = counts[sub]
        mult = MULTIPLIER.get(n, f"{n}") if n > 1 else ""
        if all_same:
            # all identical: "trimethoxy" (no extra parens needed)
            parts.append(f"{mult}{sub}")
        else:
            # mixed: each gets its own parentheses for clarity
            inner = sub[1:-1] if sub.startswith("(") else sub
            if mult:
                parts.append(f"{mult}({inner})")
            else:
                parts.append(f"({inner})")
    return "".join(parts) + "borane"


def _name_boronate_ester(graph, pgrp, get_atom) -> str:
    """ボロン酸エステル: diR' Rboronate  (Phase 255)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter
    b_idx = pgrp.atom_indices[0]
    c_on_b = [nb for nb in graph.adjacency[b_idx] if get_atom(graph, nb).symbol == "C"]
    o_esters = [nb for nb in graph.adjacency[b_idx]
                if get_atom(graph, nb).symbol == "O"
                and any(get_atom(graph, onc).symbol == "C"
                        for onc in graph.adjacency[nb] if onc != b_idx)]
    # Alkyl on B (the boronate part)
    b_alkyl = _name_carbon_substituent(graph, c_on_b[0], {b_idx}) if c_on_b else ""
    # Ester alkyls
    ester_names = sorted(
        _name_carbon_substituent(graph, [nb for nb in graph.adjacency[o] if nb != b_idx and get_atom(graph, nb).symbol == "C"][0], {o})
        for o in o_esters
        if any(get_atom(graph, nb).symbol == "C" for nb in graph.adjacency[o] if nb != b_idx)
    )
    counts = Counter(ester_names)
    parts = []
    for sub in sorted(counts):
        n = counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    return " ".join(parts) + " " + b_alkyl + "boronate"


def _name_boronic_acid(graph, pgrp, get_atom) -> str:
    """ボロン酸: {alkyl/aryl}boronic acid  (Phase 143)"""
    from .substituent import _name_carbon_substituent
    b_idx = pgrp.atom_indices[0]
    c_neighbors = [nb for nb in graph.adjacency[b_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_neighbors:
        return "boronic acid"
    sub = _name_carbon_substituent(graph, c_neighbors[0], {b_idx})
    return f"{sub}boronic acid"


def _name_borinic_acid(graph, pgrp, get_atom) -> str:
    """ボリン酸: di{alkyl}borinic acid  (Phase 143)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "borinic acid")


def _name_organic_borane(graph, pgrp, get_atom) -> str:
    """有機ボラン: {alkyl(s)}borane  (Phase 143)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "borane")


def _name_isocyanide(graph, pgrp, get_atom) -> str:
    """
    イソシアニド命名: isocyano{chain}ane (IUPAC 2013 P-62.5.3.2 substitutive PIN)
    例: C[N+]#[C-]   → isocyanomethane
        CC[N+]#[C-]  → isocyanoethane
        c1ccccc1[N+]#[C-] → isocyanobenzene
    """
    from .constants import CHAIN_PREFIX
    n_idx = pgrp.atom_indices[0]
    c_cn_idxs = [nb for nb in graph.adjacency[n_idx]
                 if get_atom(graph, nb).symbol == "C"
                 and get_atom(graph, nb).formal_charge == -1]
    c_alkyl_idxs = [nb for nb in graph.adjacency[n_idx]
                    if nb not in c_cn_idxs and get_atom(graph, nb).symbol == "C"]
    if not c_alkyl_idxs:
        return "hydrogen isocyanide"
    alkyl_c = c_alkyl_idxs[0]
    c_atom = get_atom(graph, alkyl_c)
    if c_atom.in_ring and c_atom.is_aromatic:
        # aromatic ring: isocyanobenzene etc. — use ring namer via substituent prefix
        from .substituent import name_substituent
        ring_name = name_substituent(graph, alkyl_c, {n_idx})
        # ring_name is e.g. "phenyl" → parent "benzene"
        _RING_SUB_TO_PARENT = {"phenyl": "benzene", "naphthyl": "naphthalene"}
        parent = _RING_SUB_TO_PARENT.get(ring_name, ring_name.replace("yl", "ene"))
        return f"isocyano{parent}"
    # acyclic chain: use substituent name to capture E/Z and unsaturation
    from .substituent import _name_carbon_substituent
    import re as _re_ic
    yl_name_ic = _name_carbon_substituent(graph, alkyl_c, {n_idx})
    # Strip E/Z stereo prefix if present
    stereo_pfx_ic2 = ""
    yl_base_ic = yl_name_ic
    _m_stereo_ic = _re_ic.match(r"^(\([^)]+\)-)", yl_name_ic)
    if _m_stereo_ic:
        stereo_pfx_ic2 = _m_stereo_ic.group(1)
        yl_base_ic = yl_name_ic[len(stereo_pfx_ic2):]
    # General: "{chain_base}-{loc}-yl" → "{loc}-isocyano{chain_base}e"
    m_gen_ic = _re_ic.match(r"^(.*)-(\d+)-yl$", yl_base_ic)
    if m_gen_ic:
        chain_base_ic, loc_ic = m_gen_ic.group(1), m_gen_ic.group(2)
        return f"{stereo_pfx_ic2}{loc_ic}-isocyano{chain_base_ic}e"
    # Simple terminal: "methyl"/"ethyl" → isocyano{stem}ane
    stem_base_ic = yl_base_ic[:-2] if yl_base_ic.endswith("yl") else yl_base_ic
    return f"{stereo_pfx_ic2}isocyano{stem_base_ic}ane"


def _name_ammonium(graph, pgrp, get_atom) -> str:
    """アンモニウムイオン: tetraRazanium  (Phase 146)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "azanium")


def _name_phosphanium(graph, pgrp, get_atom) -> str:
    """ホスホニウム: tetraRphosphanium (IUPAC 2013 P-73.2, Phase 518)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "phosphanium")


def _name_sulfonium(graph, pgrp, get_atom) -> str:
    """スルホニウム: triRsulfonium (IUPAC 2013 P-73.4, Phase 518)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "sulfonium")


def _name_arsonium(graph, pgrp, get_atom) -> str:
    """アルソニウム: tetraRarsonium (IUPAC 2013 P-73.3, Phase 518)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "arsonium")


def _name_organic_silane(graph, pgrp, get_atom) -> str:
    """有機シラン: {alkyl(s)}silane or {halo}({alkyl})silane  (Phase 143/249)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter
    _halo_map = {"Cl": "chloro", "Br": "bromo", "F": "fluoro", "I": "iodo"}
    central = pgrp.atom_indices[0]
    c_nbrs = [nb for nb in graph.adjacency[central] if get_atom(graph, nb).symbol == "C"]
    hal_nbrs = [nb for nb in graph.adjacency[central] if get_atom(graph, nb).symbol in _halo_map]
    if not c_nbrs:
        return "silane"
    if not hal_nbrs:
        return _name_by_c_substituents(graph, pgrp, get_atom, "silane")
    # Mixed: halogen + alkyl substituents — sort alphabetically, alkyl in parens
    halo_names = sorted(_halo_map[get_atom(graph, h).symbol] for h in hal_nbrs)
    halo_counts = Counter(halo_names)
    alkyl_names = sorted(_name_carbon_substituent(graph, c, {central}) for c in c_nbrs)
    alkyl_counts = Counter(alkyl_names)
    parts: list[tuple[str, str]] = []
    for sub in sorted(halo_counts):
        n = halo_counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append((sub, f"{mult}{sub}"))
    for sub in sorted(alkyl_counts):
        n = alkyl_counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append((sub, f"({mult}{sub})"))
    parts.sort(key=lambda x: x[0])
    return "".join(p for _, p in parts) + "silane"


def _name_silyl_ether(graph, pgrp, get_atom) -> str:
    """シリルエーテル命名: {alkoxy(s)}{alkyl(s)}silane (Phase 376)"""
    from .substituent import _name_carbon_substituent, _make_oxy_name
    from .constants import MULTIPLIER
    from collections import Counter
    central = pgrp.atom_indices[0]
    c_nbrs = [nb for nb in graph.adjacency[central] if get_atom(graph, nb).symbol == "C"]
    o_nbrs = [nb for nb in graph.adjacency[central] if get_atom(graph, nb).symbol == "O"]
    alkoxy_groups: list[str] = []
    for o_idx in o_nbrs:
        c_on_o = next(
            (nb for nb in graph.adjacency[o_idx]
             if nb != central and get_atom(graph, nb).symbol == "C"),
            None,
        )
        if c_on_o is not None:
            alkyl = _name_carbon_substituent(graph, c_on_o, {o_idx})
            alkoxy_groups.append(_make_oxy_name(alkyl))
    alkyl_groups = [_name_carbon_substituent(graph, c, {central}) for c in c_nbrs]
    all_groups = sorted(alkoxy_groups) + sorted(alkyl_groups)
    cnt = Counter(all_groups)
    parts: list[tuple[str, str]] = []
    for sub in sorted(cnt):
        n = cnt[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append((sub, f"{mult}{sub}"))
    parts.sort(key=lambda x: x[0])
    return "".join(p for _, p in parts) + "silane"


def _name_disilane(graph, pgrp, get_atom) -> str:
    """ジシラン命名: {alkyl(s)}disilane (Phase 379)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter
    si1_idx = pgrp.atom_indices[0]
    si2_idx = pgrp.atom_indices[1]
    c1 = [nb for nb in graph.adjacency[si1_idx] if get_atom(graph, nb).symbol == "C"]
    c2 = [nb for nb in graph.adjacency[si2_idx] if get_atom(graph, nb).symbol == "C"]
    names = sorted(
        [_name_carbon_substituent(graph, c, {si1_idx}) for c in c1]
        + [_name_carbon_substituent(graph, c, {si2_idx}) for c in c2]
    )
    if not names:
        return "disilane"
    counts = Counter(names)
    parts = []
    for sub in sorted(counts):
        n = counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    return "".join(parts) + "disilane"


def _name_disiloxane(graph, pgrp, get_atom) -> str:
    """ジシロキサン命名: {alkyl(s)}disiloxane (Phase 378)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter
    si1_idx = pgrp.atom_indices[0]
    si2_idx = pgrp.atom_indices[1]
    c1 = [nb for nb in graph.adjacency[si1_idx] if get_atom(graph, nb).symbol == "C"]
    c2 = [nb for nb in graph.adjacency[si2_idx] if get_atom(graph, nb).symbol == "C"]
    names = sorted(
        [_name_carbon_substituent(graph, c, {si1_idx}) for c in c1]
        + [_name_carbon_substituent(graph, c, {si2_idx}) for c in c2]
    )
    if not names:
        return "disiloxane"
    counts = Counter(names)
    parts = []
    for sub in sorted(counts):
        n = counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    return "".join(parts) + "disiloxane"


def _name_disilazane(graph, pgrp, get_atom) -> str:
    """ジシラザン命名: {alkyl(s)}disilazane (Phase 378)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter
    si1_idx = pgrp.atom_indices[0]
    si2_idx = pgrp.atom_indices[1]
    c1 = [nb for nb in graph.adjacency[si1_idx] if get_atom(graph, nb).symbol == "C"]
    c2 = [nb for nb in graph.adjacency[si2_idx] if get_atom(graph, nb).symbol == "C"]
    names = sorted(
        [_name_carbon_substituent(graph, c, {si1_idx}) for c in c1]
        + [_name_carbon_substituent(graph, c, {si2_idx}) for c in c2]
    )
    if not names:
        return "disilazane"
    counts = Counter(names)
    parts = []
    for sub in sorted(counts):
        n = counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    return "".join(parts) + "disilazane"


def _name_organic_germane(graph, pgrp, get_atom) -> str:
    """有機ゲルマン: {alkyl(s)}germane  (Phase 243)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "germane")


def _name_organic_stannane(graph, pgrp, get_atom) -> str:
    """有機スタンナン: {alkyl(s)}stannane  (Phase 243)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "stannane")


def _name_arsane_org(graph, pgrp, get_atom) -> str:
    """有機ヒ化水素: {alkyl(s)}arsane  (Phase 245)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "arsane")


def _name_organomercury(graph, pgrp, get_atom) -> str:
    """有機水銀: di{alkyl}mercury or {halo}({alkyl})mercury  (Phase 245)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter
    hg_idx = pgrp.atom_indices[0]
    c_nbrs = [nb for nb in graph.adjacency[hg_idx] if get_atom(graph, nb).symbol == "C"]
    hal_nbrs = [nb for nb in graph.adjacency[hg_idx]
                if get_atom(graph, nb).symbol in ("Cl", "Br", "F", "I")]
    _halo_names = {"Cl": "chloro", "Br": "bromo", "F": "fluoro", "I": "iodo"}
    if not c_nbrs:
        return "mercury"
    alkyl_names = sorted(_name_carbon_substituent(graph, c, {hg_idx}) for c in c_nbrs)
    counts = Counter(alkyl_names)
    # Build alkyl part with multipliers
    alkyl_parts = []
    for sub in sorted(counts):
        n = counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        alkyl_parts.append(f"({mult}{sub})" if hal_nbrs else f"{mult}{sub}")
    if hal_nbrs:
        halo_parts = sorted(_halo_names.get(get_atom(graph, h).symbol, "halo") for h in hal_nbrs)
        # Sort all substituents alphabetically (strip parens for sorting)
        all_subs = [(p.strip("()"), p) for p in halo_parts + alkyl_parts]
        all_subs.sort(key=lambda x: x[0])
        return "".join(p for _, p in all_subs) + "mercury"
    return "".join(alkyl_parts) + "mercury"


def _name_organic_bismuthane(graph, pgrp, get_atom) -> str:
    """有機ビスマン: {alkyl(s)}bismuthane  (Phase 244)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "bismuthane")


def _name_organic_stibane(graph, pgrp, get_atom) -> str:
    """有機スチバン: {alkyl(s)}stibane  (Phase 244)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "stibane")


def _name_organic_plumbane(graph, pgrp, get_atom) -> str:
    """有機プルンバン: {alkyl(s)}plumbane  (Phase 244)"""
    return _name_by_c_substituents(graph, pgrp, get_atom, "plumbane")


def _name_chalcogen_oxyacid(graph, pgrp, get_atom, suffix: str) -> str:
    """Se/Te オキソ酸命名: {stem}ane{suffix}  (Phase 244)"""
    from .constants import CHAIN_PREFIX
    central = pgrp.atom_indices[0]
    c_on_x = [nb for nb in graph.adjacency[central] if get_atom(graph, nb).symbol == "C"]
    if not c_on_x:
        return suffix
    chain, locant = _chain_through_pivot(graph, c_on_x[0], {central}, get_atom)
    stem = CHAIN_PREFIX.get(len(chain), f"C{len(chain)}")
    _ene, _yne = _chain_multiple_bonds(graph, chain)
    if _ene or _yne:
        from .name_assembler import _format_multiple_bonds as _fmt
        mb_str = _fmt(_ene, _yne)
        if len(chain) >= 3:
            return f"{stem}{mb_str}e-{locant}-{suffix}"
        return f"{stem}{mb_str}e{suffix}"
    if len(chain) >= 3:
        return f"{stem}ane-{locant}-{suffix}"
    return f"{stem}ane{suffix}"


def _name_selenonic_acid(graph, pgrp, get_atom) -> str:
    """セレノン酸: {stem}aneselenonicacid  (Phase 244)"""
    return _name_chalcogen_oxyacid(graph, pgrp, get_atom, "selenonic acid")


def _name_seleninic_acid(graph, pgrp, get_atom) -> str:
    """セレニン酸: {stem}aneseleninic acid  (Phase 244)"""
    return _name_chalcogen_oxyacid(graph, pgrp, get_atom, "seleninic acid")


def _name_selenenic_acid(graph, pgrp, get_atom) -> str:
    """セレネン酸: {stem}aneselenenic acid  (Phase 244)"""
    return _name_chalcogen_oxyacid(graph, pgrp, get_atom, "selenenic acid")


def _name_telluronic_acid(graph, pgrp, get_atom) -> str:
    """テルロン酸: {stem}anetelluronic acid  (Phase 244)"""
    return _name_chalcogen_oxyacid(graph, pgrp, get_atom, "telluronic acid")


def _name_tellurinic_acid(graph, pgrp, get_atom) -> str:
    """テルリン酸: {stem}anetellurinic acid  (Phase 244)"""
    return _name_chalcogen_oxyacid(graph, pgrp, get_atom, "tellurinic acid")


def _name_tellurenic_acid(graph, pgrp, get_atom) -> str:
    """テルレン酸: {stem}anetellurenic acid  (Phase 244)"""
    return _name_chalcogen_oxyacid(graph, pgrp, get_atom, "tellurenic acid")


def _name_organic_silanol(graph, pgrp, get_atom) -> str:
    """シラノール/ジオール/トリオール: {alkyl(s)}silanol/silanediol/silanetriol (Phase 231/379)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter
    _OH_SUFFIX = {1: "silanol", 2: "silanediol", 3: "silanetriol", 4: "silanetetraol"}
    si_idx = pgrp.atom_indices[0]
    c_neighbors = [nb for nb in graph.adjacency[si_idx]
                   if get_atom(graph, nb).symbol == "C"]
    o_oh = [nb for nb in graph.adjacency[si_idx]
            if get_atom(graph, nb).symbol == "O"
            and any(get_atom(graph, onh).symbol == "H" for onh in graph.adjacency[nb])]
    suffix = _OH_SUFFIX.get(len(o_oh), "silanol")
    if not c_neighbors:
        return suffix
    alkyl_names = [_name_carbon_substituent(graph, c, {si_idx}) for c in c_neighbors]
    counts = Counter(alkyl_names)
    parts = []
    for sub in sorted(counts):
        n = counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    return "".join(parts) + suffix


def _name_carbamate(graph, pgrp, get_atom) -> str:
    """
    カルバメート命名: {alkyl} carbamate / {alkyl} N-{sub}carbamate (Phase 42)
    例: NC(=O)OC → methyl carbamate
        CNC(=O)OCC → ethyl N-methylcarbamate
        CN(C)C(=O)OC → methyl N,N-dimethylcarbamate
    """
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter

    carbonyl_c = pgrp.atom_indices[0]
    # O-R 側 (ester O): 単結合 O で C 隣接あり
    ester_o: int | None = None
    n_idx: int | None = None
    for nb_idx in graph.adjacency[carbonyl_c]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "O" and nb_idx in pgrp.atom_indices:
            from .molecule_analyzer import get_bond_order as _gbo
            if _gbo(graph, carbonyl_c, nb_idx) == 1.0:
                ester_o = nb_idx
        elif nb.symbol == "N" and nb_idx in pgrp.atom_indices:
            n_idx = nb_idx

    # アルキル基 (O 側)
    alkyl_name = "methyl"
    if ester_o is not None:
        c_on_o = [nb for nb in graph.adjacency[ester_o]
                  if nb != carbonyl_c and get_atom(graph, nb).symbol == "C"]
        if c_on_o:
            alkyl_name = _name_carbon_substituent(graph, c_on_o[0], {ester_o, carbonyl_c})

    base = f"{alkyl_name} carbamate"

    # N 置換基
    if n_idx is None:
        return base

    # Phase 390: N が環内にある場合 → {alkyl} {ring_name}-N-carboxylate (Phase 390)
    n_atom = get_atom(graph, n_idx)
    if n_atom.in_ring:
        ring_set = next(
            (rt for rt in (graph.ring_atom_sets or []) if n_idx in rt), None
        )
        if ring_set is not None:
            from .heterocycle_handler import (
                name_heterocycle as _nhc,
                _find_best_start as _fbs_cb,
                _build_locant_map as _blm_cb,
            )
            from .ring_handler import (
                _assign_ring_locants, collect_ring_substituents, assemble_ring_name,
            )
            ring_list = list(ring_set)
            is_all_c = all(get_atom(graph, a).symbol == "C" for a in ring_set)
            if is_all_c:
                ring_chain = _assign_ring_locants(graph, ring_list, False, "alkane", [carbonyl_c])
                ring_subs = collect_ring_substituents(graph, ring_chain, [carbonyl_c])
                ring_base = assemble_ring_name(ring_chain, ring_subs, "alkane", None, [])
                n_locant = next(
                    (i + 1 for i, a in enumerate(ring_chain) if a == n_idx), None
                )
                if n_locant is not None:
                    return f"{alkyl_name} {ring_base}-{n_locant}-carboxylate"
                return f"{alkyl_name} {ring_base}carboxylate"
            else:
                rot = _fbs_cb(ring_list, graph)
                lmap = _blm_cb(rot)
                n_locant = lmap.get(n_idx)
                # Build heterocycle base name from retained names
                from .heterocycle_handler import (
                    _is_aromatic_ring as _iar_cb,
                    _canonical_sig as _csig_cb,
                    _RETAINED_NAMES as _RN_cb,
                )
                is_arom = _iar_cb(ring_list, graph)
                csig = _csig_cb(rot, graph)
                retained = _RN_cb.get((is_arom, csig))
                if retained is not None:
                    rbase, is_nh = retained
                    if n_locant is not None:
                        return f"{alkyl_name} {rbase}-{n_locant}-carboxylate"
                    return f"{alkyl_name} {rbase}carboxylate"
                else:
                    from .heterocycle_handler import (
                        _match_hantzsch_widman as _mhw_cb,
                    )
                    hw_name = _mhw_cb(ring_list, graph)
                    if hw_name is not None:
                        if n_locant is not None:
                            return f"{alkyl_name} {hw_name}-{n_locant}-carboxylate"
                        return f"{alkyl_name} {hw_name}carboxylate"

    # Phase 258: N-N パターン (ヒドラジノカルボキシレート)
    n_n_nbrs = [nb for nb in graph.adjacency[n_idx]
                if nb != carbonyl_c and get_atom(graph, nb).symbol == "N"]
    if n_n_nbrs:
        return f"{alkyl_name} hydrazinecarboxylate"

    n_c_nbrs = [nb for nb in graph.adjacency[n_idx]
                if nb != carbonyl_c and get_atom(graph, nb).symbol == "C"]
    if not n_c_nbrs:
        return base

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in n_c_nbrs]
    sub_counts = Counter(n_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if sub.startswith("(") else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub_str}")
    n_prefix = "-".join(prefix_parts)
    return f"{alkyl_name} {n_prefix}carbamate"


def _name_o_thiocarbamate(graph, pgrp, get_atom) -> str:
    """
    O-アルキルチオカルバメート命名 (Phase 350): O-alkyl carbamothioate
    例: CCOC(=S)N  → O-ethyl carbamothioate
        CCOC(=S)NC → O-ethyl N-methylcarbamothioate
    """
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from .functional_group import get_bond_order
    from collections import Counter

    central_c = pgrp.atom_indices[0]
    # Find the O-ester atom
    o_idx: int | None = None
    n_idx: int | None = None
    for ai in pgrp.atom_indices[1:]:
        atom = get_atom(graph, ai)
        if atom.symbol == "O":
            o_idx = ai
        elif atom.symbol == "N":
            n_idx = ai

    # Name the O-alkyl substituent
    alkyl_name = "methyl"
    if o_idx is not None:
        c_on_o = [nb for nb in graph.adjacency[o_idx]
                  if nb != central_c and get_atom(graph, nb).symbol == "C"]
        if c_on_o:
            alkyl_name = _name_carbon_substituent(graph, c_on_o[0], {o_idx, central_c})

    base = f"O-{alkyl_name} carbamothioate"

    # N-substituents
    if n_idx is None:
        return base
    n_c_nbrs = [nb for nb in graph.adjacency[n_idx]
                if nb != central_c and get_atom(graph, nb).symbol == "C"]
    if not n_c_nbrs:
        return base

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in n_c_nbrs]
    sub_counts = Counter(n_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if sub.startswith("(") else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub_str}")
    n_prefix = "-".join(prefix_parts)
    return f"O-{alkyl_name} {n_prefix}carbamothioate"


def _name_s_carbamothioate(graph, pgrp, get_atom) -> str:
    """
    S-アルキルカルバモチオアート命名 (Phase 351): S-alkyl carbamothioate
    例: CCSC(=O)N  → S-ethyl carbamothioate
        CCSC(=O)NC → S-ethyl N-methylcarbamothioate
    """
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from .functional_group import get_bond_order
    from collections import Counter

    central_c = pgrp.atom_indices[0]
    s_idx: int | None = None
    n_idx: int | None = None
    for ai in pgrp.atom_indices[1:]:
        atom = get_atom(graph, ai)
        if atom.symbol == "S":
            s_idx = ai
        elif atom.symbol == "N":
            n_idx = ai

    alkyl_name = "methyl"
    if s_idx is not None:
        c_on_s = [nb for nb in graph.adjacency[s_idx]
                  if nb != central_c and get_atom(graph, nb).symbol == "C"]
        if c_on_s:
            alkyl_name = _name_carbon_substituent(graph, c_on_s[0], {s_idx, central_c})

    base = f"S-{alkyl_name} carbamothioate"
    if n_idx is None:
        return base
    n_c_nbrs = [nb for nb in graph.adjacency[n_idx]
                if nb != central_c and get_atom(graph, nb).symbol == "C"]
    if not n_c_nbrs:
        return base

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in n_c_nbrs]
    sub_counts = Counter(n_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if sub.startswith("(") else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub_str}")
    n_prefix = "-".join(prefix_parts)
    return f"S-{alkyl_name} {n_prefix}carbamothioate"


def _name_s_carbamodithioate(graph, pgrp, get_atom) -> str:
    """
    S-アルキルカルバモジチオアート命名 (Phase 351): S-alkyl carbamodithioate
    例: CSC(=S)N  → S-methyl carbamodithioate
        CSC(=S)NC → S-methyl N-methylcarbamodithioate
    """
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from .functional_group import get_bond_order
    from collections import Counter

    central_c = pgrp.atom_indices[0]
    s_ester: int | None = None
    n_idx: int | None = None
    for ai in pgrp.atom_indices[2:]:  # [0]=C, [1]=S_double, rest=N,S_ester
        atom = get_atom(graph, ai)
        if atom.symbol == "S":
            s_ester = ai
        elif atom.symbol == "N":
            n_idx = ai

    alkyl_name = "methyl"
    if s_ester is not None:
        c_on_s = [nb for nb in graph.adjacency[s_ester]
                  if nb != central_c and get_atom(graph, nb).symbol == "C"]
        if c_on_s:
            alkyl_name = _name_carbon_substituent(graph, c_on_s[0], {s_ester, central_c})

    base = f"S-{alkyl_name} carbamodithioate"
    if n_idx is None:
        return base
    n_c_nbrs = [nb for nb in graph.adjacency[n_idx]
                if nb != central_c and get_atom(graph, nb).symbol == "C"]
    if not n_c_nbrs:
        return base

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in n_c_nbrs]
    sub_counts = Counter(n_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if sub.startswith("(") else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub_str}")
    n_prefix = "-".join(prefix_parts)
    return f"S-{alkyl_name} {n_prefix}carbamodithioate"


def _name_peroxide(graph, pgrp, get_atom) -> str | None:
    """有機ペルオキシド命名: {dialkyl} peroxide (Phase 57)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    c1_idx = pgrp.atom_indices[0]
    o1_idx = pgrp.atom_indices[1]
    o2_idx = pgrp.atom_indices[2]
    c2_idx = pgrp.atom_indices[3]

    # Phase 260: O-O が環内にある場合 (1,2-dioxane 等) → 環命名へ委譲
    if get_atom(graph, o1_idx).in_ring and get_atom(graph, o2_idx).in_ring:
        return None

    alkyl1 = _name_carbon_substituent(graph, c1_idx, {o1_idx})
    alkyl2 = _name_carbon_substituent(graph, c2_idx, {o2_idx})

    if alkyl1 == alkyl2:
        mult = MULTIPLIER.get(2, "di")
        return f"{mult}{alkyl1} peroxide"
    return f"{min(alkyl1, alkyl2)} {max(alkyl1, alkyl2)} peroxide"


def _name_thioester(graph, pgrp, get_atom) -> str:
    """チオエステル命名: S-{alkyl} {acid}thioate (Phase 55)"""
    from .substituent import _name_carbon_substituent
    from .constants import CHAIN_PREFIX
    s_idx = pgrp.atom_indices[0]
    carbonyl_c = pgrp.atom_indices[1]
    alkyl_cs = pgrp.atom_indices[2:]

    # S-アルキル基
    s_alkyl = "methyl"
    if alkyl_cs:
        s_alkyl = _name_carbon_substituent(graph, alkyl_cs[0], {s_idx})

    # 芳香環直結 (benzene or heteroaromatic)
    for _rn_ts in graph.adjacency[carbonyl_c]:
        if _rn_ts == s_idx:
            continue
        _rna_ts = get_atom(graph, _rn_ts)
        if not (_rna_ts.symbol == "C" and _rna_ts.in_ring and _rna_ts.is_aromatic):
            continue
        _apfx_ts = _aryl_sulfonyl_prefix(graph, _rn_ts, carbonyl_c, get_atom)
        if _apfx_ts is not None:
            return f"S-{s_alkyl} {_apfx_ts}carbothioate"
        break

    # 酸鎖 (カルボニル C から S を除外した方向)
    acid_chain = _collect_acid_chain(graph, carbonyl_c, {s_idx}, get_atom)
    n_acid = len(acid_chain)
    stem = CHAIN_PREFIX.get(n_acid, f"C{n_acid}")
    _ene_ts, _yne_ts = _chain_multiple_bonds(graph, acid_chain)
    if _ene_ts or _yne_ts:
        from .name_assembler import _format_multiple_bonds as _fmt_ts
        acid_name = f"{stem}{_fmt_ts(_ene_ts, _yne_ts)}ethioate"
    else:
        acid_name = f"{stem}anethioate"

    # E/Z 立体化学
    stereo_prefix = ""
    if _ene_ts:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_ts = PrincipalChain(atom_indices=acid_chain,
                                locant_map={c: i + 1 for i, c in enumerate(acid_chain)})
        _stereo_ts = assign_stereochemistry(graph, _pc_ts)
        if _stereo_ts:
            _comb_ts = ",".join(d.strip("()") for d in _stereo_ts)
            stereo_prefix = f"({_comb_ts})-"

    return f"{stereo_prefix}S-{s_alkyl} {acid_name}"


def _name_o_thioester(graph, pgrp, get_atom) -> str:
    """O-チオエステル: O-{alkyl} {chain}thioate  (Phase 251)"""
    from .substituent import _name_carbon_substituent
    from .constants import CHAIN_PREFIX
    o_idx = pgrp.atom_indices[0]
    carbonyl_c = pgrp.atom_indices[1]
    alkyl_cs = pgrp.atom_indices[2:]
    o_alkyl = _name_carbon_substituent(graph, alkyl_cs[0], {o_idx}) if alkyl_cs else "methyl"
    for _rn_ot in graph.adjacency[carbonyl_c]:
        if _rn_ot == o_idx:
            continue
        _rna_ot = get_atom(graph, _rn_ot)
        if not (_rna_ot.symbol == "C" and _rna_ot.in_ring and _rna_ot.is_aromatic):
            continue
        _apfx_ot = _aryl_sulfonyl_prefix(graph, _rn_ot, carbonyl_c, get_atom)
        if _apfx_ot is not None:
            return f"O-{o_alkyl} {_apfx_ot}carbothioate"
        break
    acid_chain = _collect_acid_chain(graph, carbonyl_c, {o_idx}, get_atom)
    n_acid = len(acid_chain)
    stem = CHAIN_PREFIX.get(n_acid, f"C{n_acid}")
    return f"O-{o_alkyl} {stem}anethioate"


def _name_s_dithioate_ester(graph, pgrp, get_atom) -> str:
    """S-ジチオエステル: S-{alkyl} {chain}dithioate  (Phase 251)"""
    from .substituent import _name_carbon_substituent
    from .constants import CHAIN_PREFIX
    s_ester_idx = pgrp.atom_indices[0]
    carbonyl_c = pgrp.atom_indices[1]
    alkyl_cs = pgrp.atom_indices[2:]
    s_alkyl = _name_carbon_substituent(graph, alkyl_cs[0], {s_ester_idx}) if alkyl_cs else "methyl"
    for _rn_sd in graph.adjacency[carbonyl_c]:
        if _rn_sd == s_ester_idx:
            continue
        _rna_sd = get_atom(graph, _rn_sd)
        if not (_rna_sd.symbol == "C" and _rna_sd.in_ring and _rna_sd.is_aromatic):
            continue
        _apfx_sd = _aryl_sulfonyl_prefix(graph, _rn_sd, carbonyl_c, get_atom)
        if _apfx_sd is not None:
            return f"S-{s_alkyl} {_apfx_sd}carbodithioate"
        break
    acid_chain = _collect_acid_chain(graph, carbonyl_c, {s_ester_idx}, get_atom)
    n_acid = len(acid_chain)
    stem = CHAIN_PREFIX.get(n_acid, f"C{n_acid}")
    return f"S-{s_alkyl} {stem}anedithioate"


def _name_disulfide(graph, pgrp, get_atom) -> str:
    """ジスルフィド命名: IUPAC 2013 P-63.7.1 dialkyl disulfide 形式"""
    import re as _re_ds
    from .substituent import _name_carbon_substituent, name_substituent

    s1_idx = pgrp.atom_indices[0]
    s2_idx = pgrp.atom_indices[1]
    c1_atoms = pgrp.atom_indices[2:]
    c2_atoms = []
    for nb in graph.adjacency[s2_idx]:
        if nb != s1_idx and get_atom(graph, nb).symbol == "C":
            c2_atoms.append(nb)
            break

    if not c1_atoms or not c2_atoms:
        return "disulfide"

    c1, c2 = c1_atoms[0], c2_atoms[0]

    def _group_name(c_idx: int, excl_s: int) -> str:
        atom = get_atom(graph, c_idx)
        if atom.in_ring and atom.is_aromatic:
            return name_substituent(graph, c_idx, {excl_s}) or "phenyl"
        return _name_carbon_substituent(graph, c_idx, {excl_s})

    name1 = _group_name(c1, s1_idx)
    name2 = _group_name(c2, s2_idx)

    def _needs_parens(nm: str) -> bool:
        return bool(_re_ds.search(r"[0-9]", nm)) or nm.startswith("(")

    def _alpha_base(nm: str) -> str:
        s = _re_ds.sub(r"^\(", "", nm)
        s = _re_ds.sub(r"^(di|tri|tetra|bis|tris)", "", s)
        return s.lower()

    names = sorted([name1, name2], key=_alpha_base)

    if names[0] == names[1]:
        nm = names[0]
        prefix = f"bis({nm})" if _needs_parens(nm) else f"di{nm}"
    else:
        parts = [f"({nm})" if _needs_parens(nm) else nm for nm in names]
        prefix = " ".join(parts)

    return f"{prefix} disulfide"


def _name_polysulfide(graph, pgrp, get_atom) -> str:
    """ポリスルフィド命名: dialkyl tri/tetrasulfide (Phase 226)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    atom_indices = pgrp.atom_indices
    # atom_indices = [s1, s2, ..., sn, c1, c2]
    s_count = sum(1 for i in atom_indices if get_atom(graph, i).symbol == "S")
    c_indices = [i for i in atom_indices if get_atom(graph, i).symbol == "C"]
    if len(c_indices) < 2:
        return pgrp.group_type  # fallback
    alkyl1 = _name_carbon_substituent(graph, c_indices[0], {atom_indices[0]})
    alkyl2 = _name_carbon_substituent(graph, c_indices[1], {atom_indices[s_count - 1]})
    _SUFFIX = {2: "disulfide", 3: "trisulfide", 4: "tetrasulfide"}
    suffix = _SUFFIX.get(s_count, f"{s_count}sulfide")
    if alkyl1 == alkyl2:
        mult = MULTIPLIER.get(2, "di")
        return f"{mult}{alkyl1} {suffix}"
    return f"{min(alkyl1, alkyl2)} {max(alkyl1, alkyl2)} {suffix}"


def _name_selenide_telluride(graph, pgrp, get_atom) -> str:
    """セレニド / テルリド: IUPAC 2013 dialkyl selenide/telluride 形式"""
    import re as _re_ch
    from .substituent import _name_carbon_substituent, name_substituent

    chalcogen_idx = pgrp.atom_indices[0]
    chalcogen = get_atom(graph, chalcogen_idx)
    type_word = "selenide" if chalcogen.symbol == "Se" else "telluride"

    c_neighbors = [nb for nb in graph.adjacency[chalcogen_idx]
                   if get_atom(graph, nb).symbol == "C"]
    if len(c_neighbors) < 2:
        return type_word

    c1, c2 = c_neighbors[0], c_neighbors[1]

    def _group_name(c_idx: int) -> str:
        atom = get_atom(graph, c_idx)
        if atom.in_ring and atom.is_aromatic:
            return name_substituent(graph, c_idx, {chalcogen_idx}) or "phenyl"
        return _name_carbon_substituent(graph, c_idx, {chalcogen_idx})

    name1 = _group_name(c1)
    name2 = _group_name(c2)

    def _needs_parens(nm: str) -> bool:
        return bool(_re_ch.search(r"[0-9]", nm)) or nm.startswith("(")

    def _alpha_base(nm: str) -> str:
        s = _re_ch.sub(r"^\(", "", nm)
        s = _re_ch.sub(r"^(di|tri|tetra|bis|tris)", "", s)
        return s.lower()

    names = sorted([name1, name2], key=_alpha_base)

    if names[0] == names[1]:
        nm = names[0]
        prefix = f"bis({nm})" if _needs_parens(nm) else f"di{nm}"
    else:
        parts = [f"({nm})" if _needs_parens(nm) else nm for nm in names]
        prefix = " ".join(parts)

    return f"{prefix} {type_word}"


def _name_diselenide_ditelluride(graph, pgrp, get_atom) -> str:
    """ジセレニド / ジテルリド: IUPAC 2013 dialkyl diselenide/ditelluride 形式"""
    import re as _re_dch
    from .substituent import _name_carbon_substituent, name_substituent

    se1_idx = pgrp.atom_indices[0]
    se2_idx = pgrp.atom_indices[1]
    chalcogen = get_atom(graph, se1_idx)
    type_word = "diselenide" if chalcogen.symbol == "Se" else "ditelluride"

    c1_atoms = pgrp.atom_indices[2:]
    c2_atoms = []
    for nb in graph.adjacency[se2_idx]:
        if nb != se1_idx and get_atom(graph, nb).symbol == "C":
            c2_atoms.append(nb)
            break
    if not c1_atoms or not c2_atoms:
        return type_word

    c1, c2 = c1_atoms[0], c2_atoms[0]

    def _group_name(c_idx: int, excl: int) -> str:
        atom = get_atom(graph, c_idx)
        if atom.in_ring and atom.is_aromatic:
            return name_substituent(graph, c_idx, {excl}) or "phenyl"
        return _name_carbon_substituent(graph, c_idx, {excl})

    name1 = _group_name(c1, se1_idx)
    name2 = _group_name(c2, se2_idx)

    def _needs_parens(nm: str) -> bool:
        return bool(_re_dch.search(r"[0-9]", nm)) or nm.startswith("(")

    def _alpha_base(nm: str) -> str:
        s = _re_dch.sub(r"^\(", "", nm)
        s = _re_dch.sub(r"^(di|tri|tetra|bis|tris)", "", s)
        return s.lower()

    names = sorted([name1, name2], key=_alpha_base)

    if names[0] == names[1]:
        nm = names[0]
        prefix = f"bis({nm})" if _needs_parens(nm) else f"di{nm}"
    else:
        parts = [f"({nm})" if _needs_parens(nm) else nm for nm in names]
        prefix = " ".join(parts)

    return f"{prefix} {type_word}"


def _name_nitroso(graph, pgrp, get_atom) -> str:
    """ニトロソ化合物: nitroso{alkane} (Phase 52)"""
    from .constants import CHAIN_PREFIX
    alkyl_c = pgrp.atom_indices[0]
    n_idx = pgrp.atom_indices[1]
    c_atom = get_atom(graph, alkyl_c)
    # 芳香環に直結している場合: nitroso + ring name
    if c_atom.is_aromatic and c_atom.in_ring:
        ring_atoms = set()
        for rt in (graph.ring_atom_sets or []):
            if alkyl_c in rt:
                ring_atoms.update(rt)
        if all(get_atom(graph, a).symbol == "C" for a in ring_atoms) and len(ring_atoms) == 6:
            return "nitrosobenzene"
        # Heteroaromatic ring: use _aryl_sulfonyl_prefix for locant-aware naming
        _apfx_ns = _aryl_sulfonyl_prefix(graph, alkyl_c, n_idx, get_atom)
        if _apfx_ns is not None:
            _parts_ns = _apfx_ns.rstrip("-").rsplit("-", 1)
            if len(_parts_ns) == 2:
                _rname_ns, _loc_ns = _parts_ns
                _sep_ns = "-" if _rname_ns and _rname_ns[0].isdigit() else ""
                return f"{_loc_ns}-nitroso{_sep_ns}{_rname_ns}"
        return "nitrosobenzene"
    # 非芳香族環に直結している場合: nitroso + cycloalkane name (Phase 387)
    if c_atom.in_ring and not c_atom.is_aromatic:
        ring_set = next(
            (rt for rt in (graph.ring_atom_sets or []) if alkyl_c in rt), None
        )
        if ring_set is not None and all(get_atom(graph, a).symbol == "C" for a in ring_set):
            from .ring_handler import (
                _assign_ring_locants,
                collect_ring_substituents,
                assemble_ring_name,
            )
            ring_list = list(ring_set)
            ring_chain = _assign_ring_locants(graph, ring_list, False, "alkane", [n_idx])
            ring_subs = collect_ring_substituents(graph, ring_chain, [n_idx])
            ring_base = assemble_ring_name(ring_chain, ring_subs, "alkane", None, [])
            return f"nitroso{ring_base}"
    # 全炭素鎖を pivot アルゴリズムで取得し、ロカントを決定
    chain, locant = _chain_through_pivot(graph, alkyl_c, {n_idx}, get_atom)
    n_c = len(chain)
    stem = CHAIN_PREFIX.get(n_c, f"C{n_c}")
    if n_c >= 3:
        return f"{locant}-nitroso{stem}ane"
    return f"nitroso{stem}ane"


def _name_azide(graph, pgrp, get_atom) -> str:
    """アジド: azido{alkane} PIN (IUPAC 2013 P-65.3.1 置換命名)"""
    alkyl_c = pgrp.atom_indices[0]
    n_idx = pgrp.atom_indices[1]
    return _name_isocyanate_substitutive(graph, alkyl_c, n_idx, "azido", get_atom)


def _name_isocyanate(graph, pgrp, get_atom) -> str:
    """イソシアネート: {R} isocyanate (IUPAC 2013 P-65.3.1)"""
    alkyl_c = pgrp.atom_indices[0]
    n_idx = pgrp.atom_indices[1]
    from .substituent import _name_carbon_substituent
    yl = _name_carbon_substituent(graph, alkyl_c, {n_idx})
    return f"{yl} isocyanate"


def _name_isothiocyanate(graph, pgrp, get_atom) -> str:
    """イソチオシアネート: {R} isothiocyanate (IUPAC 2013 P-65.5.1.2)"""
    alkyl_c = pgrp.atom_indices[0]
    n_idx = pgrp.atom_indices[1]
    from .substituent import _name_carbon_substituent
    yl = _name_carbon_substituent(graph, alkyl_c, {n_idx})
    return f"{yl} isothiocyanate"


def _name_isocyanate_substitutive(graph, alkyl_c: int, n_idx: int, prefix: str, get_atom) -> str:
    """芳香族・脂肪族のイソシアナト/イソチオシアナト命名 (置換命名)。"""
    from .constants import CHAIN_PREFIX
    c_atom = get_atom(graph, alkyl_c)
    if c_atom.in_ring and c_atom.is_aromatic:
        # Heteroaromatic ring: use _aryl_sulfonyl_prefix for locant-aware naming
        _ra_ic = next((rt for rt in (graph.ring_atom_sets or []) if alkyl_c in rt), None)
        if _ra_ic and any(get_atom(graph, a).symbol != "C" for a in _ra_ic):
            _apfx_ic = _aryl_sulfonyl_prefix(graph, alkyl_c, n_idx, get_atom)
            if _apfx_ic is not None:
                # "pyridine-2-" → ring_name="pyridine", loc="2"
                _parts_ic = _apfx_ic.rstrip("-").rsplit("-", 1)
                if len(_parts_ic) == 2:
                    _rname_ic, _loc_ic = _parts_ic
                    _sep_ic = "-" if _rname_ic and _rname_ic[0].isdigit() else ""
                    return f"{_loc_ic}-{prefix}{_sep_ic}{_rname_ic}"
                return f"{prefix}{_apfx_ic.rstrip('-')}"
        # 芳香環親化合物として命名 (isocyanatobenzene 等)
        from .ring_handler import (
            find_rings, find_principal_ring,
            collect_ring_substituents, assemble_ring_name,
        )
        rings = find_rings(graph)
        if not rings:
            rings = [list(rt) for rt in (graph.ring_atom_sets or []) if alkyl_c in rt]
        if rings:
            ring_chain = find_principal_ring(graph, rings, "alkane", [])
            # isocyanato 基は substituent として収集される (substituent.py で対応済み)
            subs = collect_ring_substituents(graph, ring_chain, [])
            return assemble_ring_name(ring_chain, subs, "alkane", None, [])
    # 脂肪族: isocyanato{stem}ane / {loc}-isocyanato{stem}ane
    from .substituent import _name_carbon_substituent
    yl_name = _name_carbon_substituent(graph, alkyl_c, {n_idx})
    import re as _re2
    # Strip E/Z stereo prefix if present: "(2E)-but-2-en-1-yl" → stereo_pfx="(2E)-", rest="but-2-en-1-yl"
    stereo_pfx_ic = ""
    yl_base = yl_name
    _m_stereo = _re2.match(r"^(\([^)]+\)-)", yl_name)
    if _m_stereo:
        stereo_pfx_ic = _m_stereo.group(1)
        yl_base = yl_name[len(stereo_pfx_ic):]
    # General: "{chain_base}-{loc}-yl" → "{loc}-{prefix}{chain_base}e"
    m_gen = _re2.match(r"^(.*)-(\d+)-yl$", yl_base)
    if m_gen:
        chain_base, loc = m_gen.group(1), m_gen.group(2)
        return f"{stereo_pfx_ic}{loc}-{prefix}{chain_base}e"
    # Simple terminal: "methyl"/"ethyl"/"propyl" → {prefix}{stem}ane
    stem_base = yl_base[:-2] if yl_base.endswith("yl") else yl_base
    return f"{stereo_pfx_ic}{prefix}{stem_base}ane"


def _name_di_xisocyanate(graph, pgrp, get_atom, prefix: str) -> str:
    """ジイソシアナト/ジイソチオシアナト命名: X,Y-di{prefix}{alkane} (Phase 300)"""
    from .constants import CHAIN_PREFIX
    # Collect the two alkyl C atoms (the C bonded to each isocyanate N)
    n_atoms = [ai for ai in pgrp.atom_indices
               if get_atom(graph, ai).symbol == "N"]
    if len(n_atoms) < 2:
        return f"di{prefix}methane"
    # For each N, find the alkyl C (C bonded to N but not the =C= of isocyanate)
    alkyl_cs = []
    for n_idx in n_atoms:
        for nb in graph.adjacency[n_idx]:
            nb_atom = get_atom(graph, nb)
            if nb_atom.symbol == "C" and not any(
                get_atom(graph, nb2).symbol in ("O", "S")
                and _get_bond_order(graph, nb, nb2) == 2.0
                for nb2 in graph.adjacency[nb]
                if nb2 != n_idx
            ):
                alkyl_cs.append((nb, n_idx))
                break
    if len(alkyl_cs) < 2:
        return f"di{prefix}methane"
    c1, n1 = alkyl_cs[0]
    c2, n2 = alkyl_cs[1]
    chain, _ = _chain_through_pivot(graph, c1, {n1}, get_atom)
    if c2 not in chain:
        chain, _ = _chain_through_pivot(graph, c2, {n2}, get_atom)
    lm_fwd = {c: i + 1 for i, c in enumerate(chain)}
    lm_rev = {c: len(chain) - i for i, c in enumerate(chain)}
    locs_fwd = sorted(lm_fwd.get(c, 0) for c, _ in alkyl_cs)
    locs_rev = sorted(lm_rev.get(c, 0) for c, _ in alkyl_cs)
    locs = locs_fwd if locs_fwd <= locs_rev else locs_rev
    stem = CHAIN_PREFIX.get(len(chain), f"C{len(chain)}")
    loc_str = ",".join(str(l) for l in locs)
    return f"{loc_str}-di{prefix}{stem}ane"


def _get_bond_order(graph, a: int, b: int) -> float:
    """Helper to get bond order between atoms a and b."""
    from .molecule_analyzer import get_bond_order as _gbo
    return _gbo(graph, a, b)


def _name_diisocyanate(graph, pgrp, get_atom) -> str:
    """ジイソシアネート命名: X,Y-diisocyanato{alkane} (Phase 300)"""
    return _name_di_xisocyanate(graph, pgrp, get_atom, "isocyanato")


def _name_diisothiocyanate(graph, pgrp, get_atom) -> str:
    """ジイソチオシアネート命名: X,Y-diisothiocyanato{alkane} (Phase 300)"""
    return _name_di_xisocyanate(graph, pgrp, get_atom, "isothiocyanato")


def _name_cyanate(graph, pgrp, get_atom) -> str:
    """シアン酸エステル: cyanato{alkane} PIN (IUPAC 2013 P-65.3.1); HOCN → cyanic acid"""
    cyano_c = pgrp.atom_indices[0]
    n_idx = pgrp.atom_indices[1] if len(pgrp.atom_indices) > 1 else None
    o_idx = pgrp.atom_indices[2] if len(pgrp.atom_indices) > 2 else None
    if o_idx is not None:
        for nb in graph.adjacency[o_idx]:
            if nb != cyano_c and get_atom(graph, nb).symbol == "C":
                return _name_isocyanate_substitutive(graph, nb, o_idx, "cyanato", get_atom)
        return "cyanic acid"
    return "cyanate"


def _name_thiocyanate(graph, pgrp, get_atom) -> str:
    """チオシアン酸エステル: thiocyanato{alkane} PIN (IUPAC 2013 P-65.3.1); HSC≡N → thiocyanic acid"""
    cyano_c = pgrp.atom_indices[0]
    n_idx = pgrp.atom_indices[1] if len(pgrp.atom_indices) > 1 else None
    s_idx = pgrp.atom_indices[2] if len(pgrp.atom_indices) > 2 else None
    if s_idx is not None:
        for nb in graph.adjacency[s_idx]:
            if nb != cyano_c and get_atom(graph, nb).symbol == "C":
                return _name_isocyanate_substitutive(graph, nb, s_idx, "thiocyanato", get_atom)
        return "thiocyanic acid"
    return "thiocyanate"


def _name_carbamic_acid(graph, pgrp, get_atom) -> str:
    """カルバミン酸: carbamic acid / N-alkylcarbamic acid (Phase 71)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER

    carbonyl_c = pgrp.atom_indices[0]
    # N 原子を見つける
    n_idx = next(
        (ai for ai in pgrp.atom_indices[1:] if get_atom(graph, ai).symbol == "N"),
        None,
    )
    if n_idx is None:
        return "carbamic acid"

    # N-N 結合がある場合: ヒドラジンカルボン酸 (Phase 304)
    n2_on_n = [nb for nb in graph.adjacency[n_idx]
               if nb != carbonyl_c and get_atom(graph, nb).symbol == "N"]
    if n2_on_n:
        n2_idx = n2_on_n[0]
        c_on_n2 = [nb for nb in graph.adjacency[n2_idx]
                   if get_atom(graph, nb).symbol == "C"]
        if not c_on_n2:
            return "hydrazinecarboxylic acid"
        from collections import Counter
        names2 = [_name_carbon_substituent(graph, c, {n2_idx}) for c in c_on_n2]
        cnt_map2 = Counter(names2)
        parts2 = []
        for sub in sorted(cnt_map2):
            cnt = cnt_map2[sub]
            if cnt == 1:
                parts2.append(f"2-{sub}")
            else:
                mult = MULTIPLIER.get(cnt, f"{cnt}")
                parts2.append(f"2,2-{mult}{sub}")
        return ",".join(parts2) + "hydrazinecarboxylic acid"

    # N の C 置換基を収集
    n_c_subs = [
        nb for nb in graph.adjacency[n_idx]
        if nb != carbonyl_c and get_atom(graph, nb).symbol == "C"
    ]
    if not n_c_subs:
        return "carbamic acid"

    from collections import Counter
    names = [_name_carbon_substituent(graph, c, {n_idx, carbonyl_c}) for c in n_c_subs]
    sub_counts = Counter(names)
    parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if sub.startswith("(") else sub
        if cnt == 1:
            parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            parts.append(f"N,N-{mult}{sub_str}")
    prefix = ",".join(parts)
    return f"{prefix}carbamic acid"


def _name_carbodiimide(graph, pgrp, get_atom) -> str:
    """
    カルボジイミド: N,N'-di{alkyl}carbodiimide (Phase 73 / Phase 116)
    例: CN=C=NC → N,N'-dimethylcarbodiimide
        CCN=C=NCC → N,N'-diethylcarbodiimide
    """
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER

    central_c = pgrp.atom_indices[0]
    n_idxs = pgrp.atom_indices[1:]  # 2 N atoms

    alkyl_names = []
    for n_idx in n_idxs:
        c_sgl = [nb for nb in graph.adjacency[n_idx]
                 if nb != central_c and get_atom(graph, nb).symbol == "C"]
        if c_sgl:
            name = _name_carbon_substituent(graph, c_sgl[0], {n_idx, central_c})
            alkyl_names.append(name)

    if not alkyl_names:
        return "carbodiimide"

    def _wrap_ci(s: str) -> str:
        return f"({s})" if s.startswith("(") else s

    if len(alkyl_names) == 2 and alkyl_names[0] == alkyl_names[1]:
        mult = MULTIPLIER.get(2, "di")
        w = _wrap_ci(alkyl_names[0])
        return f"N,N'-{mult}{w}carbodiimide"

    names_sorted = sorted(alkyl_names)
    if len(names_sorted) == 2:
        w0, w1 = _wrap_ci(names_sorted[0]), _wrap_ci(names_sorted[1])
        return f"N-{w0}-N'-{w1}carbodiimide"
    return f"N-{_wrap_ci(names_sorted[0])}carbodiimide"


def _name_hydroperoxide(graph, pgrp, get_atom) -> str:
    """
    ヒドロペルオキシド命名: {alkyl} hydroperoxide (Phase 44)
    例: COO → methyl hydroperoxide
        CCOO → ethyl hydroperoxide
    """
    from .substituent import _name_carbon_substituent
    c_idx = pgrp.atom_indices[0]
    o1_idx = pgrp.atom_indices[1] if len(pgrp.atom_indices) > 1 else None
    excluded = {o1_idx} if o1_idx is not None else set()
    alkyl = _name_carbon_substituent(graph, c_idx, excluded)
    return f"{alkyl} hydroperoxide"


def _name_thiourea_if_match(graph, get_atom) -> str | None:
    """
    チオ尿素命名 (Phase 87)。
    C(=S) を中心に N が 2 つ付き、他の C が付かない構造を検出。
    例: NC(=S)N → "thiourea", CNC(=S)N → "N-methylthiourea"
    """
    from .molecule_analyzer import get_bond_order
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER

    center_idx = None
    n_neighbors: list[int] = []
    for a in graph.atoms:
        if a.symbol != "C":
            continue
        if a.in_ring:  # thiolactam C is in ring; thiourea C is not
            continue
        dbl_s = False
        nc: list[int] = []
        other_c = False
        for nb_idx in graph.adjacency[a.idx]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "S" and get_bond_order(graph, a.idx, nb_idx) == 2.0:
                dbl_s = True
            elif nb.symbol == "N":
                nc.append(nb_idx)
            elif nb.symbol == "C":
                other_c = True
        if dbl_s and len(nc) == 2 and not other_c:
            # N-N パターン (チオセミカルバゾン等) は除外
            n_n_bond = any(
                get_atom(graph, x).symbol == "N"
                for n_nb in nc
                for x in graph.adjacency[n_nb]
                if x != a.idx
            )
            if n_n_bond:
                continue
            center_idx = a.idx
            n_neighbors = nc
            break

    if center_idx is None:
        return None

    def get_subs(n_idx: int) -> list[str]:
        return sorted(
            _name_carbon_substituent(graph, c, {n_idx})
            for c in graph.adjacency[n_idx]
            if get_atom(graph, c).symbol == "C" and c != center_idx
        )

    def build_prefix(subs: list[str], prime: str = "") -> list[str]:
        from collections import Counter
        parts: list[str] = []
        for sub, cnt in sorted(Counter(subs).items()):
            tag = f"N{prime}"
            needs_p = sub.startswith("(") or any(c.isdigit() for c in sub)
            sub_str = f"({sub})" if needs_p else sub
            if cnt == 1:
                parts.append(f"{tag}-{sub_str}")
            else:
                mult = MULTIPLIER.get(cnt, str(cnt))
                parts.append(f"{tag},{tag}-{mult}{sub_str}")
        return parts

    n1, n2 = n_neighbors
    subs1 = get_subs(n1)
    subs2 = get_subs(n2)

    if not subs1 and not subs2:
        return "thiourea"

    if subs1 and not subs2:
        return "-".join(build_prefix(subs1)) + "thiourea"
    if subs2 and not subs1:
        return "-".join(build_prefix(subs2)) + "thiourea"

    if subs1 == subs2 and len(subs1) == 1:
        return f"N,N'-di{subs1[0]}thiourea"

    first1 = subs1[0] if subs1 else ""
    first2 = subs2[0] if subs2 else ""
    if first1 <= first2:
        parts = build_prefix(subs1, "") + build_prefix(subs2, "'")
    else:
        parts = build_prefix(subs2, "") + build_prefix(subs1, "'")
    return "-".join(parts) + "thiourea"


def _name_substituted_urea_if_match(graph, get_atom) -> str | None:
    """
    N-置換尿素命名 (Phase 83)。
    C(=O) を中心に N が 2 つ付き、他の C が付かない構造を検出して命名する。
    例: CNC(=O)N → "N-methylurea", CNC(=O)NC → "N,N'-dimethylurea"
    """
    from .molecule_analyzer import get_bond_order
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER

    # Find urea center: C with double O, 2 N, 0 other C (must be acyclic)
    center_idx = None
    n_neighbors: list[int] = []
    for a in graph.atoms:
        if a.symbol != "C" or a.in_ring:
            continue
        dbl_o = False
        nc: list[int] = []
        other_c = False
        for nb_idx in graph.adjacency[a.idx]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "O" and get_bond_order(graph, a.idx, nb_idx) == 2.0:
                dbl_o = True
            elif nb.symbol == "N":
                nc.append(nb_idx)
            elif nb.symbol == "C":
                other_c = True
        if dbl_o and len(nc) == 2 and not other_c:
            center_idx = a.idx
            n_neighbors = nc
            break

    if center_idx is None:
        return None

    def get_subs(n_idx: int) -> list[str]:
        return sorted(
            _name_carbon_substituent(graph, c, {n_idx})
            for c in graph.adjacency[n_idx]
            if get_atom(graph, c).symbol == "C" and c != center_idx
        )

    def build_prefix(subs: list[str], prime: str = "") -> list[str]:
        from collections import Counter
        parts: list[str] = []
        for sub, cnt in sorted(Counter(subs).items()):
            tag = f"N{prime}"
            needs_p = sub.startswith("(") or any(c.isdigit() for c in sub)
            sub_str = f"({sub})" if needs_p else sub
            if cnt == 1:
                parts.append(f"{tag}-{sub_str}")
            else:
                mult = MULTIPLIER.get(cnt, str(cnt))
                parts.append(f"{tag},{tag}-{mult}{sub_str}")
        return parts

    n1, n2 = n_neighbors
    subs1 = get_subs(n1)
    subs2 = get_subs(n2)

    if not subs1 and not subs2:
        return None  # unsubstituted urea — handled by _is_urea

    if subs1 and not subs2:
        return "-".join(build_prefix(subs1)) + "urea"
    if subs2 and not subs1:
        return "-".join(build_prefix(subs2)) + "urea"

    # Both Ns have substituents
    if subs1 == subs2 and len(subs1) == 1:
        return f"N,N'-di{subs1[0]}urea"

    # Different substituents — alphabetically-first sub goes to N (no prime)
    first1 = subs1[0] if subs1 else ""
    first2 = subs2[0] if subs2 else ""
    if first1 <= first2:
        parts = build_prefix(subs1, "") + build_prefix(subs2, "'")
    else:
        parts = build_prefix(subs2, "") + build_prefix(subs1, "'")
    return "-".join(parts) + "urea"


def _is_urea(graph, get_atom) -> bool:
    """H₂N-C(=O)-NH₂ パターン（尿素本体）を検出する (Phase 49)。
    N-置換尿素 (N-methyl 等) は除外する。
    """
    from .molecule_analyzer import get_bond_order
    for a in graph.atoms:
        if a.symbol != "C":
            continue
        dbl_o = None
        n_count = 0
        c_count = 0
        n_neighbors: list[int] = []
        for nb_idx in graph.adjacency[a.idx]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "O" and get_bond_order(graph, a.idx, nb_idx) == 2.0:
                dbl_o = nb_idx
            elif nb.symbol == "N":
                n_count += 1
                n_neighbors.append(nb_idx)
            elif nb.symbol == "C":
                c_count += 1
        if dbl_o is not None and n_count == 2 and c_count == 0:
            # 両 N が C 置換基を持たない (純粋な NH₂) こと
            # ただし carbonyl C 自身は除外して判定
            all_nh2 = all(
                not any(
                    get_atom(graph, nb2).symbol == "C" and nb2 != a.idx
                    for nb2 in graph.adjacency[n_nb]
                )
                for n_nb in n_neighbors
            )
            if all_nh2:
                return True
    return False


def _is_carbonohydrazide(graph, get_atom) -> bool:
    """H2N-NH-C(=O)-NH-NH2 パターン（カルボノヒドラジド）を検出 (Phase 297)。"""
    from .molecule_analyzer import get_bond_order
    for a in graph.atoms:
        if a.symbol != "C":
            continue
        dbl_o = None
        n_neighbors: list[int] = []
        c_count = 0
        for nb_idx in graph.adjacency[a.idx]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "O" and get_bond_order(graph, a.idx, nb_idx) == 2.0:
                dbl_o = nb_idx
            elif nb.symbol == "N":
                n_neighbors.append(nb_idx)
            elif nb.symbol == "C":
                c_count += 1
        if dbl_o is None or len(n_neighbors) != 2 or c_count != 0:
            continue
        # 両 N が N-N 単結合を持つこと (hydrazide パターン)
        both_hydrazide = all(
            any(
                get_atom(graph, nb2).symbol == "N"
                and get_bond_order(graph, n_nb, nb2) == 1.0
                and nb2 != a.idx
                for nb2 in graph.adjacency[n_nb]
            )
            for n_nb in n_neighbors
        )
        if both_hydrazide:
            return True
    return False


def _name_semicarbazone(graph, pgrp, get_atom) -> str:
    """
    セミカルバゾン / チオセミカルバゾン命名 (Phase 80)
    例: CC(=NNC(N)=O)C → propan-2-one semicarbazone
        CC=NNC(N)=O → ethanal semicarbazone
        CC(=NNC(N)=S)C → propan-2-one thiosemicarbazone
    """
    from .constants import CHAIN_PREFIX
    from .chain_finder import find_principal_chain

    gtype = pgrp.group_type
    is_thio = "thio" in gtype
    is_ald = gtype.startswith("ald")

    suffix = "thiosemicarbazone" if is_thio else "semicarbazone"

    hydrazone_c = pgrp.atom_indices[0]
    chain = find_principal_chain(graph, pgrp)
    n = chain.length
    stem = CHAIN_PREFIX.get(n, f"C{n}")

    _ene_sc, _yne_sc = _chain_multiple_bonds(graph, chain.atom_indices)
    stereo_pfx_sc = ""
    if _ene_sc:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_sc = PrincipalChain(atom_indices=chain.atom_indices,
                                locant_map=chain.locant_map)
        _stereo_sc = assign_stereochemistry(graph, _pc_sc)
        if _stereo_sc:
            stereo_pfx_sc = "(" + ",".join(d.strip("()") for d in _stereo_sc) + ")-"

    # Phase 349: C=N 結合の E/Z
    cn_pfx_sc = ""
    _cn_ns_sc = [ai for ai in pgrp.atom_indices[1:]
                 if get_atom(graph, ai).symbol == "N"]
    if _cn_ns_sc:
        from .stereochemistry import _get_bond_stereo
        _cn_st_sc = _get_bond_stereo(graph, hydrazone_c, _cn_ns_sc[0])
        if _cn_st_sc is not None:
            cn_pfx_sc = f"({_cn_st_sc})-"

    full_pfx = cn_pfx_sc + stereo_pfx_sc if not stereo_pfx_sc else stereo_pfx_sc  # C=C takes precedence when both present
    # Combine: C=N stereo goes first, then C=C stereo from chain
    if cn_pfx_sc and stereo_pfx_sc:
        _all = [d for d in [cn_pfx_sc.strip("-").strip("()"), stereo_pfx_sc.strip("-").strip("()")] if d]
        full_pfx = "(" + ",".join(_all) + ")-"
    elif cn_pfx_sc:
        full_pfx = cn_pfx_sc
    else:
        full_pfx = stereo_pfx_sc

    if _ene_sc or _yne_sc:
        from .name_assembler import _format_multiple_bonds as _fmt_sc
        mb_sc = _fmt_sc(_ene_sc, _yne_sc)
        if is_ald:
            return f"{full_pfx}{stem}{mb_sc}al {suffix}"
        else:
            loc = chain.locant_map.get(hydrazone_c, 2)
            return f"{full_pfx}{stem}{mb_sc}-{loc}-one {suffix}"

    if is_ald:
        return f"{full_pfx}{stem}anal {suffix}"
    else:
        loc = chain.locant_map.get(hydrazone_c, 2)
        return f"{full_pfx}{stem}an-{loc}-one {suffix}"


def _name_substituted_hydrazone(graph, pgrp, get_atom) -> str | None:
    """
    N-置換ヒドラゾン命名 (Phase 79)
    例: CC=NNC → ethanal N-methylhydrazone
        CC(=NNC)C → propan-2-one N-methylhydrazone
    N-substituent がない場合は None を返し通常パスへ。
    """
    from .constants import CHAIN_PREFIX, MULTIPLIER
    from .substituent import _name_carbon_substituent
    from .chain_finder import find_principal_chain
    from collections import Counter

    hydrazone_c = pgrp.atom_indices[0]
    # Phase 508: C=N 両原子が環内にある場合は環命名系に委譲
    if get_atom(graph, hydrazone_c).in_ring:
        n1_check = next((ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "N"), None)
        if n1_check is not None and get_atom(graph, n1_check).in_ring:
            return None

    n1_idx = next((ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "N"), None)
    if n1_idx is None:
        return None

    # N2 (terminal N, single bond from N1)
    n2_idx = next((nb for nb in graph.adjacency[n1_idx]
                   if nb != hydrazone_c and get_atom(graph, nb).symbol == "N"), None)
    if n2_idx is None:
        return None

    # N2 の C 置換基
    n2_sub_cs = [nb for nb in graph.adjacency[n2_idx]
                 if nb != n1_idx and get_atom(graph, nb).symbol == "C"]
    if not n2_sub_cs:
        return None  # 未置換 → 通常パスへ

    # find_principal_chain でロカント正確に決定
    chain = find_principal_chain(graph, pgrp)
    n = chain.length
    stem = CHAIN_PREFIX.get(n, f"C{n}")

    _ene_hz, _yne_hz = _chain_multiple_bonds(graph, chain.atom_indices)
    stereo_pfx_hz = ""
    if _ene_hz:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain as _PC_hz
        _pc_hz = _PC_hz(atom_indices=chain.atom_indices, locant_map=chain.locant_map)
        _stereo_hz = assign_stereochemistry(graph, _pc_hz)
        if _stereo_hz:
            stereo_pfx_hz = "(" + ",".join(d.strip("()") for d in _stereo_hz) + ")-"

    if _ene_hz or _yne_hz:
        from .name_assembler import _format_multiple_bonds as _fmt_hz
        mb_hz = _fmt_hz(_ene_hz, _yne_hz)
        if pgrp.group_type == "aldhydrazone":
            parent_name = f"{stem}{mb_hz}al"
        else:
            loc = chain.locant_map.get(hydrazone_c, 2)
            parent_name = f"{stem}{mb_hz}-{loc}-one"
    elif pgrp.group_type == "aldhydrazone":
        parent_name = f"{stem}anal"
    else:
        loc = chain.locant_map.get(hydrazone_c, 2)
        parent_name = f"{stem}an-{loc}-one"

    # N2 置換基名
    n2_sub_names = sorted(_name_carbon_substituent(graph, c, {n2_idx}) for c in n2_sub_cs)
    sub_count = Counter(n2_sub_names)
    sub_parts = []
    for sub_name in sorted(sub_count.keys()):
        cnt = sub_count[sub_name]
        mult = MULTIPLIER.get(cnt, "")
        sub_str_hz = f"({sub_name})" if sub_name.startswith("(") else sub_name
        sub_parts.append(f"N-{mult}{sub_str_hz}")
    n_prefix = ",".join(sub_parts) if sub_parts else ""

    # Phase 349: C=N 結合の E/Z
    from .stereochemistry import _get_bond_stereo as _gbs_hz
    _cn_st_hz = _gbs_hz(graph, hydrazone_c, n1_idx)
    cn_pfx_hz = f"({_cn_st_hz})-" if _cn_st_hz is not None else ""

    if n_prefix:
        return f"{cn_pfx_hz}{stereo_pfx_hz}{parent_name} {n_prefix}hydrazone"
    return f"{cn_pfx_hz}{stereo_pfx_hz}{parent_name} hydrazone"


def _name_peroxyacid(graph, pgrp, get_atom) -> str:
    """
    ペルオキシ酸命名: {stem}aneperoxoic acid (Phase 77)
    例: CC(=O)OO → ethaneperoxoic acid
        CCC(=O)OO → propaneperoxoic acid
        c1ccccc1C(=O)OO → benzeneperoxoic acid
    """
    from .constants import CHAIN_PREFIX

    carbonyl_c = pgrp.atom_indices[0]
    # O1 (O-O-H の外側 O) と O2 を除外して炭素鎖を収集
    excluded: set[int] = {ai for ai in pgrp.atom_indices
                          if get_atom(graph, ai).symbol == "O"}

    # 芳香族環に直接結合したカルボニル → {arene}peroxoic acid
    for _nb_idx in graph.adjacency[carbonyl_c]:
        if _nb_idx in excluded:
            continue
        _nb = get_atom(graph, _nb_idx)
        if _nb.symbol == "C" and _nb.is_aromatic and _nb.in_ring:
            _ring_atoms = next(
                (set(rt) for rt in (graph.ring_atom_sets or []) if _nb_idx in rt), set()
            )
            if (len(_ring_atoms) == 6
                    and all(get_atom(graph, a).symbol == "C" for a in _ring_atoms)):
                return "benzeneperoxoic acid"
            # Heteroaromatic
            if (any(get_atom(graph, a).symbol != "C" for a in _ring_atoms)
                    and all(get_atom(graph, a).is_aromatic for a in _ring_atoms)):
                _apfx_pa = _aryl_sulfonyl_prefix(graph, _nb_idx, carbonyl_c, get_atom)
                if _apfx_pa is not None:
                    return f"{_apfx_pa}carboperoxoic acid"

    acid_chain = _collect_acid_chain(graph, carbonyl_c, excluded, get_atom)
    stem = CHAIN_PREFIX.get(len(acid_chain), f"C{len(acid_chain)}")
    _ene_pa, _yne_pa = _chain_multiple_bonds(graph, acid_chain)
    if _ene_pa or _yne_pa:
        from .name_assembler import _format_multiple_bonds as _fmt_pa
        parent_pa = f"{stem}{_fmt_pa(_ene_pa, _yne_pa)}eperoxoic acid"
    else:
        parent_pa = f"{stem}aneperoxoic acid"
    stereo_pfx_pa = ""
    if _ene_pa:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_pa = PrincipalChain(atom_indices=acid_chain,
                                locant_map={c: i + 1 for i, c in enumerate(acid_chain)})
        _stereo_pa = assign_stereochemistry(graph, _pc_pa)
        if _stereo_pa:
            stereo_pfx_pa = "(" + ",".join(d.strip("()") for d in _stereo_pa) + ")-"
    return f"{stereo_pfx_pa}{parent_pa}"


def _name_peroxyester(graph, pgrp, get_atom) -> str:
    """
    ペルオキシエステル命名: {alkyl} {stem}aneperoxoate (Phase 393)
    例: CC(=O)OOC  → methyl ethaneperoxoate
        CC(=O)OOCC → ethyl ethaneperoxoate
    """
    from .constants import CHAIN_PREFIX
    from .substituent import collect_substituents, _name_carbon_substituent
    from .name_assembler import _build_prefix
    from .molecule_analyzer import get_bond_order

    carbonyl_c = pgrp.atom_indices[0]

    # O-O 鎖を特定してアルキル側 C を取得
    o1_idx: int | None = None
    o2_idx: int | None = None
    alkyl_c: int | None = None
    for nb_idx in graph.adjacency[carbonyl_c]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol != "O" or get_bond_order(graph, carbonyl_c, nb_idx) != 1.0:
            continue
        for o2_nb in graph.adjacency[nb_idx]:
            if o2_nb == carbonyl_c:
                continue
            if get_atom(graph, o2_nb).symbol == "O":
                o1_idx = nb_idx
                o2_idx = o2_nb
                for alc in graph.adjacency[o2_nb]:
                    if alc != nb_idx and get_atom(graph, alc).symbol == "C":
                        alkyl_c = alc
                break
        if o2_idx is not None:
            break

    if alkyl_c is None or o2_idx is None:
        return "peroxyester"

    alkyl_name = _name_carbon_substituent(graph, alkyl_c, {o2_idx})

    excluded: set[int] = {ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "O"}
    for _rn_pe in graph.adjacency[carbonyl_c]:
        if _rn_pe in excluded:
            continue
        _rna_pe = get_atom(graph, _rn_pe)
        if not (_rna_pe.symbol == "C" and _rna_pe.in_ring and _rna_pe.is_aromatic):
            continue
        _apfx_pe = _aryl_sulfonyl_prefix(graph, _rn_pe, carbonyl_c, get_atom)
        if _apfx_pe is not None:
            return f"{alkyl_name} {_apfx_pe}carboperoxoate"
        break
    acid_carbons = _collect_acid_chain(graph, carbonyl_c, excluded, get_atom)
    n_acid = len(acid_carbons)
    stem = CHAIN_PREFIX.get(n_acid, f"C{n_acid}")
    _ene_pe, _yne_pe = _chain_multiple_bonds(graph, acid_carbons)
    if _ene_pe or _yne_pe:
        from .name_assembler import _format_multiple_bonds as _fmt_pe
        acid_base = f"{stem}{_fmt_pe(_ene_pe, _yne_pe)}eperoxoate"
    else:
        acid_base = f"{stem}aneperoxoate"

    locant_map = {c: i + 1 for i, c in enumerate(acid_carbons)}
    principal_atoms: set[int] = set(pgrp.atom_indices)
    subs = collect_substituents(graph, acid_carbons, locant_map, list(principal_atoms))
    if n_acid == 1:
        subs = [(None, name) for _, name in subs]
    if subs:
        prefix = _build_prefix(subs)
        acid_base = f"{prefix}{acid_base}"

    return f"{alkyl_name} {acid_base}"


def _name_acyl_azide(graph, pgrp, get_atom) -> str:
    """
    アシルアジド命名: {acyl} azide (functional class)
    例: CC(=O)N=[N+]=[N-] → acetyl azide
        CCC(=O)N=[N+]=[N-] → propanoyl azide
    """
    from .constants import CHAIN_PREFIX
    c_idx = pgrp.atom_indices[0]
    excl = {nb for nb in graph.adjacency[c_idx] if get_atom(graph, nb).symbol in ("O", "N")}
    acid_chain = _collect_acid_chain(graph, c_idx, excl, get_atom)
    n_acid = len(acid_chain)
    stem = CHAIN_PREFIX.get(n_acid, f"C{n_acid}")
    # 保留名: formyl (1C) / acetyl (2C)
    if n_acid == 1:
        return "formyl azide"
    if n_acid == 2:
        return "acetyl azide"
    _ene_aa, _yne_aa = _chain_multiple_bonds(graph, acid_chain)
    if _ene_aa or _yne_aa:
        from .name_assembler import _format_multiple_bonds as _fmt_aa
        return f"{stem}{_fmt_aa(_ene_aa, _yne_aa)}oyl azide"
    return f"{stem}anoyl azide"


def _name_hydrazide(graph, pgrp, get_atom) -> str:
    """
    ヒドラジド命名: {stem}anohydrazide (Phase 75/178/186/257)
    例: CC(=O)NN → ethanohydrazide
        O=C(NN)c1ccccc1 → benzohydrazide
        CC(=O)NNC → N'-methylethanohydrazide
    """
    from .constants import CHAIN_PREFIX
    from .substituent import _name_carbon_substituent

    carbonyl_c = pgrp.atom_indices[0]
    n_acyl: int | None = None  # N adjacent to carbonyl (acyl side, N prefix)
    n_terminal: int | None = None  # terminal N (N' prefix)
    for nb_idx in graph.adjacency[carbonyl_c]:
        if get_atom(graph, nb_idx).symbol == "N":
            n_acyl = nb_idx
            break
    if n_acyl is not None:
        for nb in graph.adjacency[n_acyl]:
            if nb != carbonyl_c and get_atom(graph, nb).symbol == "N":
                n_terminal = nb
                break

    excluded: set[int] = set()
    if n_acyl is not None:
        excluded.add(n_acyl)
    if n_terminal is not None:
        excluded.add(n_terminal)

    # Phase 178/527: 芳香環に直接結合したカルボニル → benzo/hetaryl hydrazide
    for nb_idx in graph.adjacency[carbonyl_c]:
        if nb_idx in excluded:
            continue
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "C" and nb.is_aromatic and nb.in_ring:
            ring_atoms = next(
                (set(rt) for rt in (graph.ring_atom_sets or []) if nb_idx in rt), set()
            )
            if (len(ring_atoms) == 6
                    and all(get_atom(graph, a).symbol == "C" for a in ring_atoms)):
                return "benzohydrazide"
            # ヘテロ芳香環 → ring-N-carbohydrazide
            has_het_hy = any(get_atom(graph, a).symbol != "C" for a in ring_atoms)
            if has_het_hy and all(get_atom(graph, a).is_aromatic for a in ring_atoms):
                _apfx_hy = _aryl_sulfonyl_prefix(graph, nb_idx, carbonyl_c, get_atom)
                if _apfx_hy is not None:
                    return f"{_apfx_hy}carbohydrazide"
            break

    acid_chain = _collect_acid_chain(graph, carbonyl_c, excluded, get_atom)
    stem = CHAIN_PREFIX.get(len(acid_chain), f"C{len(acid_chain)}")
    _ene_hy, _yne_hy = _chain_multiple_bonds(graph, acid_chain)
    if _ene_hy or _yne_hy:
        from .name_assembler import _format_multiple_bonds as _fmt_hy
        base = f"{stem}{_fmt_hy(_ene_hy, _yne_hy)}ohydrazide"
    else:
        base = f"{stem}anohydrazide"

    stereo_pfx_hy = ""
    if _ene_hy:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_hy = PrincipalChain(atom_indices=acid_chain,
                                locant_map={c: i + 1 for i, c in enumerate(acid_chain)})
        _stereo_hy = assign_stereochemistry(graph, _pc_hy)
        if _stereo_hy:
            stereo_pfx_hy = "(" + ",".join(d.strip("()") for d in _stereo_hy) + ")-"

    # N-置換基 (acyl N 側, N prefix)
    n_subs: list[str] = []
    if n_acyl is not None:
        c_on_nacyl = [nb for nb in graph.adjacency[n_acyl]
                      if nb != carbonyl_c and nb != n_terminal
                      and get_atom(graph, nb).symbol == "C"]
        n_subs = [_name_carbon_substituent(graph, c, {n_acyl}) for c in c_on_nacyl]

    # N'-置換基 (terminal N 側, N' prefix)
    np_subs: list[str] = []
    if n_terminal is not None:
        c_on_nterm = [nb for nb in graph.adjacency[n_terminal]
                      if nb != n_acyl and get_atom(graph, nb).symbol == "C"]
        np_subs = [_name_carbon_substituent(graph, c, {n_terminal}) for c in c_on_nterm]

    if not n_subs and not np_subs:
        return f"{stereo_pfx_hy}{base}"

    from .constants import MULTIPLIER
    from .name_assembler import _needs_bis_tris as _nbp
    from collections import Counter
    prefix_parts: list[str] = []
    for sub in sorted(Counter(n_subs)):
        cnt = Counter(n_subs)[sub]
        sub_str = f"({sub})" if _nbp(sub) else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            prefix_parts.append(f"N,N-{MULTIPLIER.get(cnt, str(cnt))}{sub_str}")
    for sub in sorted(Counter(np_subs)):
        cnt = Counter(np_subs)[sub]
        sub_str = f"({sub})" if _nbp(sub) else sub
        if cnt == 1:
            prefix_parts.append(f"N'-{sub_str}")
        else:
            prefix_parts.append(f"N',N'-{MULTIPLIER.get(cnt, str(cnt))}{sub_str}")
    return f"{stereo_pfx_hy}" + "-".join(prefix_parts) + base


def _name_thiohydrazide(graph, pgrp, get_atom) -> str:
    """チオヒドラジド命名 (Phase 374): RC(=S)-NHNH2 → {stem}anethiohydrazide"""
    from .constants import CHAIN_PREFIX
    from .substituent import _name_carbon_substituent

    carbonyl_c = pgrp.atom_indices[0]
    n_acyl: int | None = None
    n_terminal: int | None = None
    for nb_idx in graph.adjacency[carbonyl_c]:
        if get_atom(graph, nb_idx).symbol == "N":
            n_acyl = nb_idx
            break
    if n_acyl is not None:
        for nb in graph.adjacency[n_acyl]:
            if nb != carbonyl_c and get_atom(graph, nb).symbol == "N":
                n_terminal = nb
                break

    excluded: set[int] = set()
    if n_acyl is not None:
        excluded.add(n_acyl)
    if n_terminal is not None:
        excluded.add(n_terminal)

    for nb_idx in graph.adjacency[carbonyl_c]:
        if nb_idx in excluded:
            continue
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "C" and nb.is_aromatic and nb.in_ring:
            ring_atoms = next(
                (set(rt) for rt in (graph.ring_atom_sets or []) if nb_idx in rt), set()
            )
            if (len(ring_atoms) == 6
                    and all(get_atom(graph, a).symbol == "C" for a in ring_atoms)):
                return "benzothiohydrazide"
            break

    acid_chain = _collect_acid_chain(graph, carbonyl_c, excluded, get_atom)
    stem = CHAIN_PREFIX.get(len(acid_chain), f"C{len(acid_chain)}")
    _ene_thz, _yne_thz = _chain_multiple_bonds(graph, acid_chain)
    if _ene_thz or _yne_thz:
        from .name_assembler import _format_multiple_bonds as _fmt_thz
        base = f"{stem}{_fmt_thz(_ene_thz, _yne_thz)}thiohydrazide"
    else:
        base = f"{stem}anethiohydrazide"

    n_subs: list[str] = []
    if n_acyl is not None:
        c_on_nacyl = [nb for nb in graph.adjacency[n_acyl]
                      if nb != carbonyl_c and nb != n_terminal
                      and get_atom(graph, nb).symbol == "C"]
        n_subs = [_name_carbon_substituent(graph, c, {n_acyl}) for c in c_on_nacyl]

    np_subs: list[str] = []
    if n_terminal is not None:
        c_on_nterm = [nb for nb in graph.adjacency[n_terminal]
                      if nb != n_acyl and get_atom(graph, nb).symbol == "C"]
        np_subs = [_name_carbon_substituent(graph, c, {n_terminal}) for c in c_on_nterm]

    if not n_subs and not np_subs:
        return base

    from .constants import MULTIPLIER
    from .name_assembler import _needs_bis_tris as _nbp_thz
    from collections import Counter as _Cthz
    prefix_parts: list[str] = []
    for sub in sorted(_Cthz(n_subs)):
        cnt = _Cthz(n_subs)[sub]
        sub_str = f"({sub})" if _nbp_thz(sub) else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            prefix_parts.append(f"N,N-{MULTIPLIER.get(cnt, str(cnt))}{sub_str}")
    for sub in sorted(_Cthz(np_subs)):
        cnt = _Cthz(np_subs)[sub]
        sub_str = f"({sub})" if _nbp_thz(sub) else sub
        if cnt == 1:
            prefix_parts.append(f"N'-{sub_str}")
        else:
            prefix_parts.append(f"N',N'-{MULTIPLIER.get(cnt, str(cnt))}{sub_str}")
    return "-".join(prefix_parts) + base


def _name_selenohydrazide(graph, pgrp, get_atom) -> str:
    """セレノヒドラジド命名 (Phase 374): RC(=[Se])-NHNH2 → {stem}aneselenohydrazide"""
    from .constants import CHAIN_PREFIX
    from .substituent import _name_carbon_substituent

    carbonyl_c = pgrp.atom_indices[0]
    n_acyl: int | None = None
    n_terminal: int | None = None
    for nb_idx in graph.adjacency[carbonyl_c]:
        if get_atom(graph, nb_idx).symbol == "N":
            n_acyl = nb_idx
            break
    if n_acyl is not None:
        for nb in graph.adjacency[n_acyl]:
            if nb != carbonyl_c and get_atom(graph, nb).symbol == "N":
                n_terminal = nb
                break

    excluded: set[int] = set()
    if n_acyl is not None:
        excluded.add(n_acyl)
    if n_terminal is not None:
        excluded.add(n_terminal)

    for nb_idx in graph.adjacency[carbonyl_c]:
        if nb_idx in excluded:
            continue
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "C" and nb.is_aromatic and nb.in_ring:
            ring_atoms = next(
                (set(rt) for rt in (graph.ring_atom_sets or []) if nb_idx in rt), set()
            )
            if (len(ring_atoms) == 6
                    and all(get_atom(graph, a).symbol == "C" for a in ring_atoms)):
                return "benzoselenohydrazide"
            break

    acid_chain = _collect_acid_chain(graph, carbonyl_c, excluded, get_atom)
    stem = CHAIN_PREFIX.get(len(acid_chain), f"C{len(acid_chain)}")
    _ene_shz, _yne_shz = _chain_multiple_bonds(graph, acid_chain)
    if _ene_shz or _yne_shz:
        from .name_assembler import _format_multiple_bonds as _fmt_shz
        base = f"{stem}{_fmt_shz(_ene_shz, _yne_shz)}selenohydrazide"
    else:
        base = f"{stem}aneselenohydrazide"

    n_subs: list[str] = []
    if n_acyl is not None:
        c_on_nacyl = [nb for nb in graph.adjacency[n_acyl]
                      if nb != carbonyl_c and nb != n_terminal
                      and get_atom(graph, nb).symbol == "C"]
        n_subs = [_name_carbon_substituent(graph, c, {n_acyl}) for c in c_on_nacyl]

    np_subs: list[str] = []
    if n_terminal is not None:
        c_on_nterm = [nb for nb in graph.adjacency[n_terminal]
                      if nb != n_acyl and get_atom(graph, nb).symbol == "C"]
        np_subs = [_name_carbon_substituent(graph, c, {n_terminal}) for c in c_on_nterm]

    if not n_subs and not np_subs:
        return base

    from .constants import MULTIPLIER
    from .name_assembler import _needs_bis_tris as _nbp_shz
    from collections import Counter as _Cshz
    prefix_parts: list[str] = []
    for sub in sorted(_Cshz(n_subs)):
        cnt = _Cshz(n_subs)[sub]
        sub_str = f"({sub})" if _nbp_shz(sub) else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            prefix_parts.append(f"N,N-{MULTIPLIER.get(cnt, str(cnt))}{sub_str}")
    for sub in sorted(_Cshz(np_subs)):
        cnt = _Cshz(np_subs)[sub]
        sub_str = f"({sub})" if _nbp_shz(sub) else sub
        if cnt == 1:
            prefix_parts.append(f"N'-{sub_str}")
        else:
            prefix_parts.append(f"N',N'-{MULTIPLIER.get(cnt, str(cnt))}{sub_str}")
    return "-".join(prefix_parts) + base


def _name_hetero_n_oxide(graph, get_atom) -> str | None:
    """ヘテロ芳香環 N-オキシド (Phase 165): ring N→O → "{ring_name} {locant}-oxide"。
    例: c1cc[n+]([O-])cc1 → pyridine 1-oxide
    """
    from .functional_group import get_bond_order
    from .ring_handler import find_hetero_rings, has_hetero_ring
    from .heterocycle_handler import (
        name_heterocycle, _find_best_start, _build_locant_map,
    )

    if not has_hetero_ring(graph):
        return None

    # ring N with N→O (charged or not) but no C=N double bonds
    for n_idx in range(len(graph.atoms)):
        n_atom = get_atom(graph, n_idx)
        if n_atom.symbol != "N" or not n_atom.in_ring or not n_atom.is_aromatic:
            continue
        # Find N→O
        o_nbs = [nb for nb in graph.adjacency[n_idx]
                 if get_atom(graph, nb).symbol == "O" and not get_atom(graph, nb).in_ring]
        if not o_nbs:
            continue
        o_idx = o_nbs[0]
        # O has no other non-H neighbors
        if any(get_atom(graph, nb).symbol not in ("N", "H")
               for nb in graph.adjacency[o_idx] if nb != n_idx):
            continue

        # Get the ring containing n_idx
        hetero_rings = find_hetero_rings(graph)
        ring = next((r for r in hetero_rings if n_idx in r), None)
        if ring is None:
            continue

        # Get the N locant: use _find_best_start which puts N at locant 1 for pyridine etc.
        rotation = _find_best_start(ring, graph)
        locant_map = _build_locant_map(rotation)
        n_locant = locant_map.get(n_idx, 1)

        # Get base ring name (without substituents, without N-oxide)
        # Build a modified graph description using the ring atoms only
        # Simpler: name the heterocycle but intercept the ring name from heterocycle handler
        from .heterocycle_handler import _match_retained, _canonical_sig
        sig = _canonical_sig(rotation, graph)
        from .heterocycle_handler import _RETAINED_NAMES
        is_arom = all(get_atom(graph, idx).is_aromatic for idx in ring)
        entry = _RETAINED_NAMES.get((is_arom, sig))
        if entry is None:
            continue
        base_name, _ = entry

        # Collect other substituents on the ring (excluding N-oxide O)
        from .heterocycle_handler import _collect_hetero_substituents
        subs = [(loc, nm) for loc, nm in _collect_hetero_substituents(
            graph, rotation, locant_map, excluded_atoms={o_idx}
        ) if nm not in ("", )]

        if not subs:
            # Pure ring N-oxide
            return f"{base_name} {n_locant}-oxide"
        else:
            # Ring with other substituents + N-oxide
            from .heterocycle_handler import _format_substituents
            sub_prefix = _format_substituents(base_name, subs)
            return f"{sub_prefix} {n_locant}-oxide"

    return None


def _name_amine_n_oxide(graph, get_atom) -> str | None:
    """
    アミン N-オキシドの命名 (Phase 125): R₃N→O パターン。
    例: C[N+](C)(C)[O-] → trimethylamine N-oxide
        C[N+](CC)(CC)[O-] → triethylamine N-oxide
    """
    from .functional_group import get_bond_order

    for idx in range(len(graph.atoms)):
        atom = get_atom(graph, idx)
        if atom.symbol != "N" or atom.in_ring:
            continue

        neighbors = graph.adjacency[idx]
        c_nbrs = [nb for nb in neighbors if get_atom(graph, nb).symbol == "C"]
        o_nbrs = [nb for nb in neighbors if get_atom(graph, nb).symbol == "O"]

        # N-oxide: 3個のC単結合 + 1個のO (末端、Hなし、単結合)
        if len(c_nbrs) < 1 or len(o_nbrs) != 1:
            continue

        o_idx = o_nbrs[0]
        o_atom = get_atom(graph, o_idx)

        # N-oxide O は末端 (N のみに結合) かつ H なし
        o_non_n_nbs = [nb for nb in graph.adjacency[o_idx]
                       if nb != idx and get_atom(graph, nb).symbol != "H"]
        if o_non_n_nbs:
            continue

        # O-H なし確認
        o_h_nbs = [nb for nb in graph.adjacency[o_idx] if get_atom(graph, nb).symbol == "H"]
        if o_h_nbs:
            continue

        # N-O 結合次数: 1.0 (単結合、N→O) または 2.0 (ニトロソの場合は除外済み)
        bo = get_bond_order(graph, idx, o_idx)
        if bo != 1.0:
            continue

        # N-C 結合次数: 全て 1.0 (単結合アミン)
        all_single_c = all(get_bond_order(graph, idx, c) == 1.0 for c in c_nbrs)
        if not all_single_c:
            continue

        # H 隣接 (第一級・第二級 N-oxide)
        h_nbrs = [nb for nb in neighbors if get_atom(graph, nb).symbol == "H"]

        # アミン部分を命名
        amine_name = _name_secondary_tertiary_amine(graph, idx, c_nbrs, get_atom)
        if amine_name is None:
            continue

        return f"{amine_name} N-oxide"

    return None


def _name_n_substituted_hydroxylamine(graph, get_atom) -> str | None:
    """
    N-置換ヒドロキシルアミンの命名 (Phase 202):
    R-NH-OH (一級) → N-alkylhydroxylamine
    R₂N-OH (二級) → N,N-dialkylhydroxylamine
    例: CNO  → N-methylhydroxylamine
        CCNO → N-ethylhydroxylamine
        CN(O)C → N,N-dimethylhydroxylamine
    """
    from .functional_group import get_bond_order
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter

    for idx in range(len(graph.atoms)):
        atom = get_atom(graph, idx)
        if atom.symbol != "N" or atom.in_ring:
            continue

        neighbors = graph.adjacency[idx]
        c_nbrs = [nb for nb in neighbors if get_atom(graph, nb).symbol == "C"]
        o_nbrs = [nb for nb in neighbors if get_atom(graph, nb).symbol == "O"]
        h_nbrs = [nb for nb in neighbors if get_atom(graph, nb).symbol == "H"]

        # Pattern: N has >=1 C (single bond), exactly 1 O (single bond with OH), 0-1 H
        if len(c_nbrs) < 1 or len(o_nbrs) != 1:
            continue

        o_idx = o_nbrs[0]
        # The O must be -OH (has exactly 1 H neighbor and only N as heavy neighbor)
        o_h_nbs = [nb for nb in graph.adjacency[o_idx] if get_atom(graph, nb).symbol == "H"]
        o_heavy_nbs = [nb for nb in graph.adjacency[o_idx]
                       if get_atom(graph, nb).symbol != "H" and nb != idx]
        if len(o_h_nbs) != 1 or o_heavy_nbs:
            continue

        # N-O bond must be single
        bo_no = get_bond_order(graph, idx, o_idx)
        if bo_no != 1.0:
            continue

        # All N-C bonds must be single
        if not all(get_bond_order(graph, idx, c) == 1.0 for c in c_nbrs):
            continue

        # Exclude N-substituted amides (hydroxamic acids): C neighbor has C=O
        # e.g. CC(=O)NO → N-hydroxyacetamide (handled separately, not here)
        from .functional_group import _has_double_bonded_oxygen
        if any(_has_double_bonded_oxygen(graph, c) for c in c_nbrs):
            continue

        # Exclude amidoximes: C neighbor has C=N (amidine/imine C)
        # e.g. CC(=N)NO → N'-hydroxyethanimidamide (amidoxime)
        def _has_double_bonded_nitrogen(c_idx: int) -> bool:
            for nb in graph.adjacency[c_idx]:
                if get_atom(graph, nb).symbol == "N" and get_bond_order(graph, c_idx, nb) == 2.0:
                    return True
            return False
        if any(_has_double_bonded_nitrogen(c) for c in c_nbrs):
            continue

        # Must have at least one H on N (hydroxylamine parent) or be N,N-disubstituted
        # CNO: N has 1C, 1O(OH), 1H → N-methylhydroxylamine
        # CN(O)C: N has 2C, 1O(OH), 0H → N,N-dimethylhydroxylamine
        # (but not NO itself which is handled as retained name 'hydroxylamine')

        # Name the N-substituents
        n_subs = [_name_carbon_substituent(graph, c, {idx}) for c in c_nbrs]
        counts = Counter(n_subs)

        # Build N-prefix: e.g. N-methyl, N,N-dimethyl
        n_prefix_parts = []
        for sub_name in sorted(counts.keys()):
            cnt = counts[sub_name]
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            sub_str = f"({sub_name})" if sub_name.startswith("(") else sub_name
            if cnt == 1:
                n_prefix_parts.append(f"N-{sub_str}")
            else:
                n_prefix_parts.append(f"N,N-{mult}{sub_str}")
        n_prefix = ",".join(n_prefix_parts)

        return f"{n_prefix}hydroxylamine"

    return None


def _name_o_substituted_hydroxylamine(graph, get_atom) -> str | None:
    """
    O-置換ヒドロキシルアミンの命名 (Phase 346): NH2-O-R → O-alkylhydroxylamine
    例: NOC  → O-methylhydroxylamine
        NOCC → O-ethylhydroxylamine
        NOC/C=C/C → O-((2E)-but-2-en-1-yl)hydroxylamine
    """
    from .functional_group import get_bond_order
    from .substituent import _name_carbon_substituent

    for idx in range(len(graph.atoms)):
        n_atom = get_atom(graph, idx)
        if n_atom.symbol != "N" or n_atom.in_ring:
            continue
        # N must have no C neighbors (unsubstituted N)
        c_nbrs = [nb for nb in graph.adjacency[idx] if get_atom(graph, nb).symbol == "C"]
        if c_nbrs:
            continue
        # N must have exactly one O neighbor (single bond)
        o_nbrs = [nb for nb in graph.adjacency[idx]
                  if get_atom(graph, nb).symbol == "O"
                  and get_bond_order(graph, idx, nb) == 1.0]
        if len(o_nbrs) != 1:
            continue
        o_idx = o_nbrs[0]
        # O must have exactly one C neighbor (the substituent)
        o_c_nbrs = [nb for nb in graph.adjacency[o_idx]
                    if nb != idx and get_atom(graph, nb).symbol == "C"]
        if not o_c_nbrs:
            continue
        alkyl = _name_carbon_substituent(graph, o_c_nbrs[0], {o_idx})
        alkyl_str = f"({alkyl})" if alkyl.startswith("(") else alkyl
        return f"O-{alkyl_str}hydroxylamine"

    return None


def _name_o_substituted_oxime(graph, get_atom) -> str | None:
    """
    O-置換オキシム命名 (Phase 347): C=N-O-R → O-alkyl{parent}al/one oxime
    例: CC=NOC  → O-methylethanal oxime
        CC(=NOC)C → O-methylpropan-2-one oxime
    """
    from .functional_group import get_bond_order
    from .substituent import _name_carbon_substituent
    from .constants import CHAIN_PREFIX
    from .name_assembler import _format_multiple_bonds

    for idx in range(len(graph.atoms)):
        n_atom = get_atom(graph, idx)
        if n_atom.symbol != "N" or n_atom.in_ring:
            continue
        # Find the double-bond C (oxime C)
        oxime_c = None
        for nb in graph.adjacency[idx]:
            if get_atom(graph, nb).symbol == "C" and get_bond_order(graph, idx, nb) == 2.0:
                oxime_c = nb
                break
        if oxime_c is None:
            continue
        # N must have exactly one single-bond O neighbor
        o_nbrs = [nb for nb in graph.adjacency[idx]
                  if get_atom(graph, nb).symbol == "O"
                  and get_bond_order(graph, idx, nb) == 1.0]
        if len(o_nbrs) != 1:
            continue
        o_idx = o_nbrs[0]
        # O must not have an H neighbor (that would be a regular oxime)
        if any(get_atom(graph, nb).symbol == "H" for nb in graph.adjacency[o_idx]):
            continue
        # O must have exactly one C neighbor (the O-alkyl substituent)
        o_c_nbrs = [nb for nb in graph.adjacency[o_idx]
                    if nb != idx and get_atom(graph, nb).symbol == "C"]
        if len(o_c_nbrs) != 1:
            continue

        alkyl = _name_carbon_substituent(graph, o_c_nbrs[0], {o_idx})
        alkyl_str = f"({alkyl})" if alkyl.startswith("(") else alkyl

        # aldoxime-type: oxime C has < 2 C neighbors; ketoxime-type: >= 2
        c_nbrs_of_oxime_c = [nb for nb in graph.adjacency[oxime_c]
                              if get_atom(graph, nb).symbol == "C"]
        is_ketoxime = len(c_nbrs_of_oxime_c) >= 2

        if is_ketoxime:
            chain, loc = _chain_through_pivot(graph, oxime_c, {idx}, get_atom)
        else:
            chain = _collect_acid_chain(graph, oxime_c, {idx}, get_atom)
            loc = 1

        chain_length = len(chain)
        stem = CHAIN_PREFIX.get(chain_length)
        if stem is None:
            continue

        locant_map = {c: i + 1 for i, c in enumerate(chain)}
        ene_locs, yne_locs = _chain_multiple_bonds(graph, chain)
        has_mb = bool(ene_locs or yne_locs)

        if is_ketoxime:
            loc = locant_map.get(oxime_c, loc)
            if has_mb:
                mb = _format_multiple_bonds(ene_locs, yne_locs)
                base = f"{stem}{mb}-{loc}-one oxime"
            else:
                base = f"{stem}an-{loc}-one oxime"
        else:
            if has_mb:
                mb = _format_multiple_bonds(ene_locs, yne_locs)
                base = f"{stem}{mb}al oxime"
            else:
                base = f"{stem}anal oxime"

        # Phase 349: C=N 結合の E/Z
        from .stereochemistry import _get_bond_stereo as _gbs_oox
        _cn_st_oox = _gbs_oox(graph, oxime_c, idx)
        cn_pfx_oox = f"({_cn_st_oox})-" if _cn_st_oox is not None else ""

        return f"{cn_pfx_oox}O-{alkyl_str}{base}"

    return None


def _name_nitrone(graph, get_atom) -> str | None:
    """
    ニトロン命名 (Phase 204): C=N+(R)-O- パターン (imine N-oxide).
    例: C=[N+]([O-])C  → N-methylmethanimine N-oxide
        CC=[N+]([O-])C → N-methylethan-1-imine N-oxide
    """
    from .functional_group import get_bond_order
    from .substituent import _name_carbon_substituent
    from .constants import CHAIN_PREFIX, MULTIPLIER
    from collections import Counter

    for idx in range(len(graph.atoms)):
        atom = get_atom(graph, idx)
        if atom.symbol != "N" or atom.in_ring or atom.formal_charge != 1:
            continue

        neighbors = graph.adjacency[idx]
        # Nitrone: N has one C with double bond (imine), one O with single bond (N-oxide), one C with single bond
        c_dbl = [nb for nb in neighbors if get_atom(graph, nb).symbol == "C"
                 and get_bond_order(graph, idx, nb) == 2.0]
        c_sgl = [nb for nb in neighbors if get_atom(graph, nb).symbol == "C"
                 and get_bond_order(graph, idx, nb) == 1.0]
        o_nbrs = [nb for nb in neighbors if get_atom(graph, nb).symbol == "O"]

        if len(c_dbl) != 1 or len(c_sgl) < 1 or len(o_nbrs) != 1:
            continue

        o_idx = o_nbrs[0]
        o_atom = get_atom(graph, o_idx)
        # O must have charge -1 and no heavy neighbors other than N
        if o_atom.formal_charge != -1:
            continue
        o_non_n = [nb for nb in graph.adjacency[o_idx]
                   if nb != idx and get_atom(graph, nb).symbol != "H"]
        if o_non_n:
            continue

        imine_c = c_dbl[0]

        # Build imine parent name: find chain length from imine_c
        # Walk from imine_c (excluding N) to get the chain
        visited: set[int] = {idx, o_idx}
        chain_atoms = [imine_c]
        visited.add(imine_c)
        queue = [nb for nb in graph.adjacency[imine_c]
                 if nb not in visited and get_atom(graph, nb).symbol == "C"
                 and not get_atom(graph, nb).in_ring]
        while queue:
            c = queue.pop(0)
            if c not in visited:
                visited.add(c)
                chain_atoms.append(c)
                for nb in graph.adjacency[c]:
                    if nb not in visited and get_atom(graph, nb).symbol == "C":
                        queue.append(nb)

        chain_len = len(chain_atoms)
        stem = CHAIN_PREFIX.get(chain_len, f"C{chain_len}")

        if chain_len <= 2:
            imine_name = f"{stem}animine"
        else:
            imine_name = f"{stem}an-1-imine"

        # N-substituents from c_sgl
        n_subs = [_name_carbon_substituent(graph, c, {idx}) for c in c_sgl]
        counts = Counter(n_subs)
        prefix_parts = []
        for sub in sorted(counts):
            cnt = counts[sub]
            needs_p = sub.startswith("(") or any(c.isdigit() for c in sub)
            sub_str = f"({sub})" if needs_p else sub
            if cnt == 1:
                prefix_parts.append(f"N-{sub_str}")
            else:
                mult = MULTIPLIER.get(cnt, f"{cnt}")
                prefix_parts.append(f"N,N-{mult}{sub_str}")
        n_prefix = ",".join(prefix_parts)

        return f"{n_prefix}{imine_name} N-oxide"

    return None


def _name_diazo_compound(graph, get_atom) -> str | None:
    """
    ジアゾ化合物の命名 (Phase 123): C=N=N パターン (allenic N-N).
    例: C=[N+]=[N-] → diazomethane, CC=[N+]=[N-] → diazoethane
    """
    from .functional_group import get_bond_order
    from .constants import CHAIN_PREFIX

    for idx in range(len(graph.atoms)):
        atom = get_atom(graph, idx)
        if atom.symbol != "N" or atom.in_ring:
            continue
        # この N が C と N の両方に二重結合を持つか確認
        n_neighbors = graph.adjacency[idx]
        c_double = [nb for nb in n_neighbors
                    if get_atom(graph, nb).symbol == "C"
                    and get_bond_order(graph, idx, nb) == 2.0]
        n_double = [nb for nb in n_neighbors
                    if get_atom(graph, nb).symbol == "N"
                    and get_bond_order(graph, idx, nb) == 2.0
                    and not get_atom(graph, nb).in_ring]
        if not c_double or not n_double:
            continue
        # 末端 N (もう一方の N の隣接は内側 N のみ)
        terminal_n = n_double[0]
        term_non_n_nbs = [nb for nb in graph.adjacency[terminal_n]
                          if nb != idx and get_atom(graph, nb).symbol != "H"]
        if term_non_n_nbs:
            continue  # 末端 N に他の接続があれば diazo ではない

        diazo_c = c_double[0]
        # diazo_c から C 鎖を収集 (N は除外)
        excluded = {idx}
        excluded.update(graph.adjacency[terminal_n])  # terminal N を除外

        # BFS で炭素鎖を収集
        chain: list[int] = []
        visited = {diazo_c}
        queue = [diazo_c]
        while queue:
            cur = queue.pop(0)
            chain.append(cur)
            for nb in graph.adjacency[cur]:
                if nb not in visited and get_atom(graph, nb).symbol == "C":
                    visited.add(nb)
                    queue.append(nb)

        n_c = len(chain)
        stem = CHAIN_PREFIX.get(n_c)
        if stem is None:
            return None
        return f"diazo{stem}ane"

    return None


def _name_nitrosamine(graph, get_atom) -> str | None:
    """N-ニトロソアミン命名 (Phase 164): R₂N-N=O パターン。
    例: CN(N=O)C → N-methyl-N-nitrosomethanamine
    """
    from .molecule_analyzer import get_bond_order as _gbo
    from .substituent import _name_carbon_substituent
    from .constants import CHAIN_PREFIX, MULTIPLIER

    # amine N: C 隣接あり + N 隣接(N=O)あり
    for n_idx in range(len(graph.atoms)):
        n_atom = get_atom(graph, n_idx)
        if n_atom.symbol != "N" or n_atom.in_ring:
            continue
        # N の隣接原子を分類
        c_nbs = [nb for nb in graph.adjacency[n_idx]
                 if get_atom(graph, nb).symbol == "C" and _gbo(graph, n_idx, nb) == 1.0]
        nitroso_ns = []
        for nb in graph.adjacency[n_idx]:
            if get_atom(graph, nb).symbol != "N" or nb == n_idx:
                continue
            if _gbo(graph, n_idx, nb) != 1.0:
                continue
            # 隣 N が =O を持つか
            if any(get_atom(graph, nb2).symbol == "O" and _gbo(graph, nb, nb2) == 2.0
                   for nb2 in graph.adjacency[nb]):
                nitroso_ns.append(nb)
        if not c_nbs or not nitroso_ns:
            continue
        # Found: n_idx is the amine N with N=O substituent
        # Name the amine part
        c_subs = sorted([_name_carbon_substituent(graph, c, {n_idx}) for c in c_nbs])
        from .constants import CHAIN_PREFIX
        # Parent = shortest amine chain, others are N- prefixes
        # Sort alphabetically and use the last as parent (naming convention for amines)
        c_subs_sorted = sorted(c_subs)
        if len(c_subs_sorted) == 0:
            return None
        # parent amine = first alphabetically (methanamine for methyl, ethanamine for ethyl...)
        # per IUPAC: longest chain is parent
        def _chain_len(sub_name: str) -> int:
            # Strip E/Z stereo prefix before matching
            clean = sub_name
            if sub_name.startswith("(") and ")-" in sub_name:
                clean = sub_name[sub_name.index(")-") + 2:]
            for n_c, stem in CHAIN_PREFIX.items():
                if clean.startswith(stem):
                    return n_c
            return 0
        parent_sub = max(c_subs, key=_chain_len)
        other_subs = [s for s in c_subs]
        other_subs.remove(parent_sub)

        # Build parent amine name from yl name
        import re as _re_na
        m_loc = _re_na.match(r"^(.+)-(\d+)-yl$", parent_sub)
        if m_loc:
            # "(2E)-but-2-en-1-yl" → "(2E)-but-2-en-1-amine"
            parent_amine = f"{m_loc.group(1)}-{m_loc.group(2)}-amine"
        else:
            # "methyl" → "meth" + "anamine" = "methanamine"
            parent_amine = parent_sub[:-2] + "anamine"

        # Extract stereo prefix from parent_amine to place at the front of the name
        stereo_pfx_na = ""
        parent_amine_base = parent_amine
        if parent_amine.startswith("(") and ")-" in parent_amine:
            _end_na = parent_amine.index(")-") + 2
            stereo_pfx_na = parent_amine[:_end_na]
            parent_amine_base = parent_amine[_end_na:]

        # Build N-substituent prefixes (other C subs + nitroso); alphabetical
        n_prefixes_raw = sorted(other_subs + ["nitroso"])
        n_prefixes_wrapped = [f"({p})" if p.startswith("(") else p for p in n_prefixes_raw]
        # Join as "N-{sub1}-N-{sub2}{parentamine}" (no trailing dash)
        n_part = "N-" + "-N-".join(n_prefixes_wrapped)
        return f"{stereo_pfx_na}{n_part}{parent_amine_base}"

    return None


def _name_hydrazine_compound(graph, get_atom) -> str | None:
    """
    ヒドラジン化合物の命名 (Phase 113): N-N 単結合を持つ非環状化合物。
    例: NNC → methylhydrazine, NNCC → ethylhydrazine,
        c1ccc(NN)cc1 → phenylhydrazine, NN → hydrazine,
        CN(N)C → 1,1-dimethylhydrazine
    ヒドラジド (C(=O)-N-N) とは区別される。
    """
    from .functional_group import get_bond_order
    from .substituent import _name_carbon_substituent, _name_aryl_substituent
    from .constants import MULTIPLIER

    # N-N 単結合ペアを探す
    nn_pairs: list[tuple[int, int]] = []
    for idx in range(len(graph.atoms)):
        atom = get_atom(graph, idx)
        if atom.symbol != "N" or atom.in_ring:
            continue
        for nb in graph.adjacency[idx]:
            nb_atom = get_atom(graph, nb)
            if nb_atom.symbol == "N" and not nb_atom.in_ring and nb > idx:
                bo = get_bond_order(graph, idx, nb)
                if bo == 1.0:
                    # ヒドラゾン除外: どちらかの N が C=N 二重結合を持つ場合
                    has_cn_double = False
                    for n_check in (idx, nb):
                        for nb2 in graph.adjacency[n_check]:
                            if nb2 == idx or nb2 == nb:
                                continue
                            if (get_atom(graph, nb2).symbol == "C"
                                    and get_bond_order(graph, n_check, nb2) == 2.0):
                                has_cn_double = True
                                break
                        if has_cn_double:
                            break
                    # ニトロソアミン除外: 一方の N が N=O を持つ場合 (N-N=O)
                    has_no = any(
                        get_atom(graph, nb2).symbol == "O"
                        and get_bond_order(graph, n_check, nb2) == 2.0
                        for n_check in (idx, nb)
                        for nb2 in graph.adjacency[n_check]
                        if nb2 != idx and nb2 != nb
                    )
                    if not has_cn_double and not has_no:
                        nn_pairs.append((idx, nb))

    if not nn_pairs:
        return None

    # 1ペアのみ対応
    if len(nn_pairs) != 1:
        return None

    n1_idx, n2_idx = nn_pairs[0]

    # カルボニル C に隣接している場合はヒドラジドなので除外
    for n_idx in (n1_idx, n2_idx):
        for nb in graph.adjacency[n_idx]:
            if nb == n1_idx or nb == n2_idx:
                continue
            nb_atom = get_atom(graph, nb)
            if nb_atom.symbol == "C":
                from .functional_group import _get_double_bonded_oxygen
                if _get_double_bonded_oxygen(graph, nb) is not None:
                    return None

    # 各 N の C 隣接原子を収集
    def _c_neighbors(n_idx: int) -> list[int]:
        return [nb for nb in graph.adjacency[n_idx]
                if nb not in (n1_idx, n2_idx) and get_atom(graph, nb).symbol == "C"]

    n1_c = _c_neighbors(n1_idx)
    n2_c = _c_neighbors(n2_idx)

    # ヒドラジン本体: 置換基が少ない側を「terminal NH2 側」とする
    # terminal N (NH2) は C 隣接がない or 少ない側
    if len(n1_c) <= len(n2_c):
        terminal_n, sub_n, sub_c_list = n1_idx, n2_idx, n2_c
        terminal_c = n1_c
    else:
        terminal_n, sub_n, sub_c_list = n2_idx, n1_idx, n1_c
        terminal_c = n2_c

    # 両 N に C なし → hydrazine
    if not sub_c_list and not terminal_c:
        return "hydrazine"

    # terminal 側にも C がある場合: 非対称二置換ヒドラジン
    # 例: 1,1-dimethylhydrazine: CN(N)C → sub_n に 2C, terminal_n に 0C
    # 例: 1,2-dimethylhydrazine: CNC → 1,2-…
    # 対称二置換: (CH3)2N-N(CH3)2 はスコープ外

    # 置換基を命名
    def _name_substituent(c_idx: int) -> str:
        c_atom = get_atom(graph, c_idx)
        if c_atom.in_ring and c_atom.is_aromatic:
            # アリール基
            return _name_aryl_substituent(graph, c_idx, {sub_n, terminal_n})
        return _name_carbon_substituent(graph, c_idx, {sub_n, terminal_n})

    sub_names = [_name_substituent(c) for c in sub_c_list]
    term_names = [_name_substituent(c) for c in terminal_c]

    # 1,1-二置換: sub 側に2つ、terminal 側に0つ
    if len(sub_c_list) == 2 and not terminal_c:
        if sub_names[0] == sub_names[1]:
            mult = MULTIPLIER.get(2, "di")
            return f"1,1-{mult}{sub_names[0]}hydrazine"
        s1, s2 = sorted(sub_names)
        return f"1,1-{s1}-1-{s2}hydrazine"

    # 1,2-二置換: 各 N に1つ
    if len(sub_c_list) == 1 and len(terminal_c) == 1:
        s1, s2 = sorted([sub_names[0], term_names[0]])
        if s1 == s2:
            return f"1,2-di{s1}hydrazine"
        return f"1-{s1}-2-{s2}hydrazine"

    # 一置換: sub 側に1つ
    if len(sub_c_list) == 1 and not terminal_c:
        return f"{sub_names[0]}hydrazine"

    return None


def _name_azo_compound(graph, get_atom) -> str | None:
    """
    アゾ化合物命名 (Phase 115): R-N=N-R' → azo{base} 形式。
    例: Ph-N=N-Ph → azobenzene, Me-N=N-Me → azomethane
    """
    from .functional_group import get_bond_order
    from .constants import CHAIN_PREFIX

    # N=N 二重結合ペアを探す (両方とも環外)
    nn_pairs: list[tuple[int, int]] = []
    for idx in range(len(graph.atoms)):
        atom = get_atom(graph, idx)
        if atom.symbol != "N" or atom.in_ring:
            continue
        for nb in graph.adjacency[idx]:
            nb_atom = get_atom(graph, nb)
            if nb_atom.symbol == "N" and not nb_atom.in_ring and nb > idx:
                if get_bond_order(graph, idx, nb) == 2.0:
                    nn_pairs.append((idx, nb))

    if len(nn_pairs) != 1:
        return None

    n1_idx, n2_idx = nn_pairs[0]

    def _c_neighbor(n_idx: int) -> int | None:
        """N に直結する C インデックスを返す (N=N 側は除く)"""
        for nb in graph.adjacency[n_idx]:
            if nb not in (n1_idx, n2_idx) and get_atom(graph, nb).symbol == "C":
                return nb
        return None

    def _chain_len(c_start: int) -> int:
        """c_start から N1/N2 を除外したグラフを DFS し炭素数を返す"""
        visited: set[int] = {n1_idx, n2_idx, c_start}
        stack = [c_start]
        count = 1
        while stack:
            node = stack.pop()
            for nb in graph.adjacency[node]:
                if nb in visited:
                    continue
                if get_atom(graph, nb).symbol == "C":
                    visited.add(nb)
                    stack.append(nb)
                    count += 1
        return count

    def _is_plain_benzene(n_idx: int) -> bool:
        """N が無置換ベンゼン環の C に直結しているか確認"""
        c = _c_neighbor(n_idx)
        if c is None:
            return False
        c_atom = get_atom(graph, c)
        if not (c_atom.in_ring and c_atom.is_aromatic):
            return False
        ring = next((rt for rt in (graph.ring_atom_sets or []) if c in rt), None)
        if ring is None or len(ring) != 6:
            return False
        for ra in ring:
            if get_atom(graph, ra).symbol != "C":
                return False
            for nb in graph.adjacency[ra]:
                sym = get_atom(graph, nb).symbol
                if sym not in ("C", "H", "N"):
                    return False
                if sym == "N" and nb not in (n1_idx, n2_idx):
                    return False
        return True

    r1_benzene = _is_plain_benzene(n1_idx)
    r2_benzene = _is_plain_benzene(n2_idx)

    if r1_benzene and r2_benzene:
        return "azobenzene"

    if not r1_benzene and not r2_benzene:
        c1 = _c_neighbor(n1_idx)
        c2 = _c_neighbor(n2_idx)
        if c1 is None or c2 is None:
            return None
        len1 = _chain_len(c1)
        len2 = _chain_len(c2)
        if len1 != len2 or len1 == 0:
            return None
        stem = CHAIN_PREFIX.get(len1)
        if stem is None:
            return None
        # Phase 352: E/Z on N=N bond
        from .stereochemistry import _get_bond_stereo as _gbs_azo
        _nn_stereo = _gbs_azo(graph, n1_idx, n2_idx)
        _nn_pfx = f"({_nn_stereo})-" if _nn_stereo is not None else ""
        return f"{_nn_pfx}azo{stem}ane"

    return None


def _name_thioamide(graph, pgrp, get_atom) -> str:
    """
    チオアミド命名: {stem}anethioamide / N-{sub}{stem}anethioamide (Phase 41)
    例: CC(=S)N → ethanethioamide
        CC(=S)NC → N-methylethanethioamide
        CC(=S)N(C)C → N,N-dimethylethanethioamide
    """
    from .constants import CHAIN_PREFIX, MULTIPLIER
    from .substituent import _name_carbon_substituent
    from collections import Counter

    carbonyl_c = pgrp.atom_indices[0]
    n_idx = next((ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "N"), None)
    if n_idx is None:
        return "thioamide"

    # Phase 526: ヘテロ芳香環直結チオアミド → ring-N-carbothioamide
    if not get_atom(graph, carbonyl_c).in_ring:
        for _nb_ta in graph.adjacency[carbonyl_c]:
            _nba_ta = get_atom(graph, _nb_ta)
            if _nba_ta.symbol != "C" or not _nba_ta.in_ring or not _nba_ta.is_aromatic:
                continue
            _ra_ta = next(
                (rt for rt in (graph.ring_atom_sets or []) if _nb_ta in rt), None
            )
            if _ra_ta is None:
                continue
            if all(get_atom(graph, a).symbol == "C" for a in _ra_ta):
                break  # benzene: fall through to existing return None below
            if (any(get_atom(graph, a).symbol != "C" for a in _ra_ta)
                    and all(get_atom(graph, a).is_aromatic for a in _ra_ta)):
                _apfx = _aryl_sulfonyl_prefix(graph, _nb_ta, carbonyl_c, get_atom)
                if _apfx is not None:
                    _het_base = f"{_apfx}carbothioamide"
                    _cn_ta = [nb for nb in graph.adjacency[n_idx]
                               if nb != carbonyl_c and get_atom(graph, nb).symbol == "C"]
                    if not _cn_ta:
                        return _het_base
                    _ns_ta = [_name_carbon_substituent(graph, c, {n_idx}) for c in _cn_ta]
                    _sc_ta = Counter(_ns_ta)
                    _pp_ta = []
                    for _s in sorted(_sc_ta):
                        _c = _sc_ta[_s]
                        _ss = f"({_s})" if _s.startswith("(") else _s
                        if _c == 1:
                            _pp_ta.append(f"N-{_ss}")
                        else:
                            _pp_ta.append(f"N,N-{MULTIPLIER.get(_c, str(_c))}{_ss}")
                    return f"{'-'.join(_pp_ta)}{_het_base}"
            break

    # 環に隣接する場合は ring path に委譲 (cyclopentanecarbothioamide 等)
    # N が環内の場合は ring N として後続処理で対応するため除外する (Phase 395)
    if not get_atom(graph, carbonyl_c).in_ring and any(
        get_atom(graph, nb).in_ring and get_atom(graph, nb).symbol == "C"
        for nb in graph.adjacency[carbonyl_c]
    ):
        return None

    excluded = {n_idx}
    excluded.update(nb for nb in graph.adjacency[n_idx]
                    if nb != carbonyl_c and get_atom(graph, nb).symbol == "C")
    acid_chain = _collect_acid_chain(graph, carbonyl_c, excluded, get_atom)
    stem = CHAIN_PREFIX.get(len(acid_chain), f"C{len(acid_chain)}")
    _ene_ta, _yne_ta = _chain_multiple_bonds(graph, acid_chain)
    if _ene_ta or _yne_ta:
        from .name_assembler import _format_multiple_bonds as _fmt_ta
        parent_name = f"{stem}{_fmt_ta(_ene_ta, _yne_ta)}ethioamide"
    else:
        parent_name = f"{stem}anethioamide"

    stereo_pfx_ta = ""
    if _ene_ta:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_ta = PrincipalChain(atom_indices=acid_chain,
                                locant_map={c: i + 1 for i, c in enumerate(acid_chain)})
        _stereo_ta = assign_stereochemistry(graph, _pc_ta)
        if _stereo_ta:
            stereo_pfx_ta = "(" + ",".join(d.strip("()") for d in _stereo_ta) + ")-"

    # Phase 395: N が環内にある場合 → {N_loc}-({thioamide_stem}thioyl){ring_base}
    n_atom = get_atom(graph, n_idx)
    if n_atom.in_ring:
        ring_set = next(
            (rt for rt in (graph.ring_atom_sets or []) if n_idx in rt), None
        )
        if ring_set is not None:
            from .heterocycle_handler import (
                _find_best_start as _fbs_ta,
                _build_locant_map as _blm_ta,
                _is_aromatic_ring as _iar_ta,
                _canonical_sig as _csig_ta,
                _RETAINED_NAMES as _RN_ta,
                _match_hantzsch_widman as _mhw_ta,
            )
            ring_list = list(ring_set)
            rot = _fbs_ta(ring_list, graph)
            lmap = _blm_ta(rot)
            n_locant = lmap.get(n_idx)
            is_arom = _iar_ta(ring_list, graph)
            csig = _csig_ta(rot, graph)
            retained = _RN_ta.get((is_arom, csig))
            if retained is not None:
                rbase, _ = retained
            else:
                rbase = _mhw_ta(ring_list, graph)
            if rbase is not None:
                thioyl_stem = parent_name.replace("thioamide", "thioyl")
                loc_str = f"{n_locant}-" if n_locant is not None else ""
                return f"{loc_str}({stereo_pfx_ta}{thioyl_stem}){rbase}"

    c_nbrs = [nb for nb in graph.adjacency[n_idx]
               if nb != carbonyl_c and get_atom(graph, nb).symbol == "C"]
    if not c_nbrs:
        return f"{stereo_pfx_ta}{parent_name}"

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in c_nbrs]
    sub_counts = Counter(n_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if sub.startswith("(") else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub_str}")
    prefix = "-".join(prefix_parts)
    return f"{stereo_pfx_ta}{prefix}{parent_name}"


def _name_dithioamide(graph, pgrp, get_atom) -> str:
    """ジチオアミド命名: {stem}anedithioamide / {ring}-X,Y-dicarbothioamide (Phase 306/319)
    例: NC(=S)C(=S)N → ethanedithioamide
    例: NC(=S)C1CCC(C(=S)N)CC1 → cyclohexane-1,4-dicarbothioamide
    """
    from .constants import CHAIN_PREFIX
    from .molecule_analyzer import get_bond_order

    thio_cs = []
    for ai in pgrp.atom_indices:
        if get_atom(graph, ai).symbol != "C":
            continue
        for nb_idx in graph.adjacency[ai]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "S" and get_bond_order(graph, ai, nb_idx) == 2.0:
                thio_cs.append(ai)
                break

    if len(thio_cs) < 2:
        return "dithioamide"

    c1, c2 = thio_cs[0], thio_cs[1]

    ring1, rc1, arom1 = _ring_for_exo_c(graph, c1, get_atom)
    ring2, rc2, arom2 = _ring_for_exo_c(graph, c2, get_atom)

    if ring1 is not None and ring2 is not None and ring1 == ring2:
        from .ring_handler import _assign_ring_locants
        ring_list = sorted(ring1)
        ring_chain = _assign_ring_locants(graph, ring_list, arom1, "arene" if arom1 else "alkane",
                                          [rc1, rc2])
        loc1 = ring_chain.locant_map.get(rc1, 1)
        loc2 = ring_chain.locant_map.get(rc2, 2)
        locs = sorted([loc1, loc2])
        n = len(ring1)
        ring_base = "benzene" if (arom1 and n == 6) else f"cyclo{CHAIN_PREFIX.get(n, f'C{n}')}ane"
        return f"{ring_base}-{locs[0]},{locs[1]}-dicarbothioamide"

    excl: set[int] = set()
    for ai in pgrp.atom_indices:
        if get_atom(graph, ai).symbol in ("S", "N"):
            excl.add(ai)

    chain_fwd_dta = _collect_acid_chain(graph, c1, excl, get_atom)
    chain_rev_dta = _collect_acid_chain(graph, c2, excl, get_atom)
    ene_fwd_dta, yne_fwd_dta = _chain_multiple_bonds(graph, chain_fwd_dta)
    ene_rev_dta, yne_rev_dta = _chain_multiple_bonds(graph, chain_rev_dta)
    mb_fwd_dta = sorted(ene_fwd_dta + yne_fwd_dta)
    mb_rev_dta = sorted(ene_rev_dta + yne_rev_dta)
    if mb_rev_dta and (not mb_fwd_dta or mb_rev_dta < mb_fwd_dta):
        chain_dta = chain_rev_dta
        _ene_dta, _yne_dta = ene_rev_dta, yne_rev_dta
    else:
        chain_dta = chain_fwd_dta
        _ene_dta, _yne_dta = ene_fwd_dta, yne_fwd_dta

    stem = CHAIN_PREFIX.get(len(chain_dta), f"C{len(chain_dta)}")
    if _ene_dta or _yne_dta:
        from .name_assembler import _format_multiple_bonds as _fmt_dta
        parent_dta = f"{stem}{_fmt_dta(_ene_dta, _yne_dta)}edithioamide"
    else:
        parent_dta = f"{stem}anedithioamide"

    stereo_pfx_dta = ""
    if _ene_dta:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_dta = PrincipalChain(atom_indices=chain_dta,
                                 locant_map={c: i + 1 for i, c in enumerate(chain_dta)})
        _stereo_dta = assign_stereochemistry(graph, _pc_dta)
        if _stereo_dta:
            stereo_pfx_dta = "(" + ",".join(d.strip("()") for d in _stereo_dta) + ")-"
    return f"{stereo_pfx_dta}{parent_dta}"


def _name_diselenoamide(graph, pgrp, get_atom) -> str:
    """ジセレノアミド命名: {stem}anediselenoamide / {ring}-X,Y-dicarboselenoamide (Phase 308/320)"""
    from .constants import CHAIN_PREFIX
    from .molecule_analyzer import get_bond_order

    seleno_cs = []
    for ai in pgrp.atom_indices:
        if get_atom(graph, ai).symbol != "C":
            continue
        for nb_idx in graph.adjacency[ai]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "Se" and get_bond_order(graph, ai, nb_idx) == 2.0:
                seleno_cs.append(ai)
                break

    if len(seleno_cs) < 2:
        return "diselenoamide"

    c1, c2 = seleno_cs[0], seleno_cs[1]

    ring1, rc1, arom1 = _ring_for_exo_c(graph, c1, get_atom)
    ring2, rc2, arom2 = _ring_for_exo_c(graph, c2, get_atom)

    if ring1 is not None and ring2 is not None and ring1 == ring2:
        from .ring_handler import _assign_ring_locants
        ring_list = sorted(ring1)
        ring_chain = _assign_ring_locants(graph, ring_list, arom1, "arene" if arom1 else "alkane",
                                          [rc1, rc2])
        loc1 = ring_chain.locant_map.get(rc1, 1)
        loc2 = ring_chain.locant_map.get(rc2, 2)
        locs = sorted([loc1, loc2])
        n = len(ring1)
        ring_base = "benzene" if (arom1 and n == 6) else f"cyclo{CHAIN_PREFIX.get(n, f'C{n}')}ane"
        return f"{ring_base}-{locs[0]},{locs[1]}-dicarboselenoamide"

    excl: set[int] = set()
    for ai in pgrp.atom_indices:
        if get_atom(graph, ai).symbol in ("Se", "N"):
            excl.add(ai)

    chain_fwd_dsa = _collect_acid_chain(graph, c1, excl, get_atom)
    chain_rev_dsa = _collect_acid_chain(graph, c2, excl, get_atom)
    ene_fwd_dsa, yne_fwd_dsa = _chain_multiple_bonds(graph, chain_fwd_dsa)
    ene_rev_dsa, yne_rev_dsa = _chain_multiple_bonds(graph, chain_rev_dsa)
    mb_fwd_dsa = sorted(ene_fwd_dsa + yne_fwd_dsa)
    mb_rev_dsa = sorted(ene_rev_dsa + yne_rev_dsa)
    if mb_rev_dsa and (not mb_fwd_dsa or mb_rev_dsa < mb_fwd_dsa):
        chain_dsa = chain_rev_dsa
        _ene_dsa, _yne_dsa = ene_rev_dsa, yne_rev_dsa
    else:
        chain_dsa = chain_fwd_dsa
        _ene_dsa, _yne_dsa = ene_fwd_dsa, yne_fwd_dsa

    stem = CHAIN_PREFIX.get(len(chain_dsa), f"C{len(chain_dsa)}")
    if _ene_dsa or _yne_dsa:
        from .name_assembler import _format_multiple_bonds as _fmt_dsa
        parent_dsa = f"{stem}{_fmt_dsa(_ene_dsa, _yne_dsa)}ediselenoamide"
    else:
        parent_dsa = f"{stem}anediselenoamide"

    stereo_pfx_dsa = ""
    if _ene_dsa:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_dsa = PrincipalChain(atom_indices=chain_dsa,
                                 locant_map={c: i + 1 for i, c in enumerate(chain_dsa)})
        _stereo_dsa = assign_stereochemistry(graph, _pc_dsa)
        if _stereo_dsa:
            stereo_pfx_dsa = "(" + ",".join(d.strip("()") for d in _stereo_dsa) + ")-"
    return f"{stereo_pfx_dsa}{parent_dsa}"


def _ring_for_exo_c(graph, c_idx: int, get_atom):
    """exocyclic C に隣接する全炭素環を返す: (ring_set, ring_C_idx, is_aromatic)"""
    for nb_idx in graph.adjacency[c_idx]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "C" and nb.in_ring:
            ring_set = next(
                (set(rt) for rt in (graph.ring_atom_sets or []) if nb_idx in rt), None
            )
            if ring_set is not None and all(
                get_atom(graph, a).symbol == "C" for a in ring_set
            ):
                return ring_set, nb_idx, nb.is_aromatic
    return None, None, False


def _name_dinitrile(graph, pgrp, get_atom):
    """ベンゼン/シクロアルカン環上の2つのニトリル基 → {ring}-X,Y-dicarbonitrile (Phase 309/318)"""
    from .molecule_analyzer import get_bond_order
    from .constants import CHAIN_PREFIX

    nitrile_cs = []
    for ai in pgrp.atom_indices:
        if get_atom(graph, ai).symbol != "C":
            continue
        for nb_idx in graph.adjacency[ai]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "N" and get_bond_order(graph, ai, nb_idx) == 3.0:
                nitrile_cs.append(ai)
                break

    if len(nitrile_cs) < 2:
        return None

    c1, c2 = nitrile_cs[0], nitrile_cs[1]

    ring1, rc1, arom1 = _ring_for_exo_c(graph, c1, get_atom)
    ring2, rc2, arom2 = _ring_for_exo_c(graph, c2, get_atom)

    if ring1 is not None and ring2 is not None and ring1 == ring2:
        from .ring_handler import _assign_ring_locants
        ring_list = sorted(ring1)
        ring_chain = _assign_ring_locants(graph, ring_list, arom1, "arene" if arom1 else "alkane",
                                          [rc1, rc2])
        loc1 = ring_chain.locant_map.get(rc1, 1)
        loc2 = ring_chain.locant_map.get(rc2, 2)
        locs = sorted([loc1, loc2])
        n = len(ring1)
        ring_base = "benzene" if (arom1 and n == 6) else f"cyclo{CHAIN_PREFIX.get(n, f'C{n}')}ane"
        return f"{ring_base}-{locs[0]},{locs[1]}-dicarbonitrile"

    return None


def _name_diamide(graph, pgrp, get_atom):
    """ベンゼン/シクロアルカン環上の2つのアミド基 → {ring}-X,Y-dicarboxamide (Phase 307/318)"""
    from .molecule_analyzer import get_bond_order
    from .constants import CHAIN_PREFIX

    amide_cs = []
    for ai in pgrp.atom_indices:
        if get_atom(graph, ai).symbol != "C":
            continue
        for nb_idx in graph.adjacency[ai]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "O" and get_bond_order(graph, ai, nb_idx) == 2.0:
                amide_cs.append(ai)
                break

    if len(amide_cs) < 2:
        return None

    c1, c2 = amide_cs[0], amide_cs[1]

    ring1, rc1, arom1 = _ring_for_exo_c(graph, c1, get_atom)
    ring2, rc2, arom2 = _ring_for_exo_c(graph, c2, get_atom)

    if ring1 is not None and ring2 is not None and ring1 == ring2:
        from .ring_handler import _assign_ring_locants
        ring_list = sorted(ring1)
        ring_chain = _assign_ring_locants(graph, ring_list, arom1, "arene" if arom1 else "alkane",
                                          [rc1, rc2])
        loc1 = ring_chain.locant_map.get(rc1, 1)
        loc2 = ring_chain.locant_map.get(rc2, 2)
        locs = sorted([loc1, loc2])
        n = len(ring1)
        ring_base = "benzene" if (arom1 and n == 6) else f"cyclo{CHAIN_PREFIX.get(n, f'C{n}')}ane"
        return f"{ring_base}-{locs[0]},{locs[1]}-dicarboxamide"

    return None


def _name_diamidine(graph, pgrp, get_atom):
    """ジアミジン命名: {ring}-X,Y-dicarboximidamide / {stem}anediimidamide (Phase 320)"""
    from .constants import CHAIN_PREFIX
    from .molecule_analyzer import get_bond_order

    amidine_cs = []
    for ai in pgrp.atom_indices:
        if get_atom(graph, ai).symbol != "C":
            continue
        n_nbs = [nb for nb in graph.adjacency[ai] if get_atom(graph, nb).symbol == "N"]
        if len(n_nbs) >= 2:
            amidine_cs.append(ai)

    if len(amidine_cs) < 2:
        return None

    c1, c2 = amidine_cs[0], amidine_cs[1]

    ring1, rc1, arom1 = _ring_for_exo_c(graph, c1, get_atom)
    ring2, rc2, arom2 = _ring_for_exo_c(graph, c2, get_atom)

    if ring1 is not None and ring2 is not None and ring1 == ring2:
        from .ring_handler import _assign_ring_locants
        ring_list = sorted(ring1)
        ring_chain = _assign_ring_locants(graph, ring_list, arom1, "arene" if arom1 else "alkane",
                                          [rc1, rc2])
        loc1 = ring_chain.locant_map.get(rc1, 1)
        loc2 = ring_chain.locant_map.get(rc2, 2)
        locs = sorted([loc1, loc2])
        n = len(ring1)
        ring_base = "benzene" if (arom1 and n == 6) else f"cyclo{CHAIN_PREFIX.get(n, f'C{n}')}ane"
        return f"{ring_base}-{locs[0]},{locs[1]}-dicarboximidamide"

    return None


def _name_selenoamide(graph, pgrp, get_atom) -> str:
    """セレノアミド命名: {stem}aneselenoamide / N-{sub}{stem}aneselenoamide (Phase 296)"""
    from .constants import CHAIN_PREFIX, MULTIPLIER
    from .substituent import _name_carbon_substituent
    from collections import Counter

    carbonyl_c = pgrp.atom_indices[0]
    n_idx = next((ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "N"), None)
    if n_idx is None:
        return "selenoamide"

    # Phase 527: ヘテロ芳香環直結セレノアミド → ring-N-carboselenoamide
    if not get_atom(graph, carbonyl_c).in_ring:
        for _nb_sea in graph.adjacency[carbonyl_c]:
            _nba_sea = get_atom(graph, _nb_sea)
            if _nba_sea.symbol != "C" or not _nba_sea.in_ring or not _nba_sea.is_aromatic:
                continue
            _ra_sea = next(
                (rt for rt in (graph.ring_atom_sets or []) if _nb_sea in rt), None
            )
            if _ra_sea is None:
                continue
            if all(get_atom(graph, a).symbol == "C" for a in _ra_sea):
                break  # benzene: fall through to existing return None below
            if (any(get_atom(graph, a).symbol != "C" for a in _ra_sea)
                    and all(get_atom(graph, a).is_aromatic for a in _ra_sea)):
                _apfx_sea = _aryl_sulfonyl_prefix(graph, _nb_sea, carbonyl_c, get_atom)
                if _apfx_sea is not None:
                    _het_base_sea = f"{_apfx_sea}carboselenoamide"
                    _cn_sea = [nb for nb in graph.adjacency[n_idx]
                                if nb != carbonyl_c and get_atom(graph, nb).symbol == "C"]
                    if not _cn_sea:
                        return _het_base_sea
                    _ns_sea = [_name_carbon_substituent(graph, c, {n_idx}) for c in _cn_sea]
                    _sc_sea = Counter(_ns_sea)
                    _pp_sea = []
                    for _s in sorted(_sc_sea):
                        _c = _sc_sea[_s]
                        _ss = f"({_s})" if _s.startswith("(") else _s
                        if _c == 1:
                            _pp_sea.append(f"N-{_ss}")
                        else:
                            _pp_sea.append(f"N,N-{MULTIPLIER.get(_c, str(_c))}{_ss}")
                    return f"{'-'.join(_pp_sea)}{_het_base_sea}"
            break

    if not get_atom(graph, carbonyl_c).in_ring and any(
        get_atom(graph, nb).in_ring for nb in graph.adjacency[carbonyl_c]
    ):
        return None

    excluded = {n_idx}
    excluded.update(nb for nb in graph.adjacency[n_idx]
                    if nb != carbonyl_c and get_atom(graph, nb).symbol == "C")
    acid_chain = _collect_acid_chain(graph, carbonyl_c, excluded, get_atom)
    stem = CHAIN_PREFIX.get(len(acid_chain), f"C{len(acid_chain)}")
    _ene_sa, _yne_sa = _chain_multiple_bonds(graph, acid_chain)
    if _ene_sa or _yne_sa:
        from .name_assembler import _format_multiple_bonds as _fmt_sa
        parent_name = f"{stem}{_fmt_sa(_ene_sa, _yne_sa)}eselenoamide"
    else:
        parent_name = f"{stem}aneselenoamide"

    stereo_pfx_sa = ""
    if _ene_sa:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_sa = PrincipalChain(atom_indices=acid_chain,
                                locant_map={c: i + 1 for i, c in enumerate(acid_chain)})
        _stereo_sa = assign_stereochemistry(graph, _pc_sa)
        if _stereo_sa:
            stereo_pfx_sa = "(" + ",".join(d.strip("()") for d in _stereo_sa) + ")-"

    c_nbrs = [nb for nb in graph.adjacency[n_idx]
               if nb != carbonyl_c and get_atom(graph, nb).symbol == "C"]
    if not c_nbrs:
        return f"{stereo_pfx_sa}{parent_name}"

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in c_nbrs]
    sub_counts = Counter(n_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if sub.startswith("(") else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub_str}")
    prefix = "-".join(prefix_parts)
    return f"{stereo_pfx_sa}{prefix}{parent_name}"


def _name_telluramide(graph, pgrp, get_atom) -> str:
    """テルラミド命名: {stem}aneteluramide / N-{sub}{stem}aneteluramide (Phase 298)"""
    from .constants import CHAIN_PREFIX, MULTIPLIER
    from .substituent import _name_carbon_substituent
    from collections import Counter

    carbonyl_c = pgrp.atom_indices[0]
    n_idx = next((ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "N"), None)
    if n_idx is None:
        return "teluramide"

    if not get_atom(graph, carbonyl_c).in_ring and any(
        get_atom(graph, nb).in_ring for nb in graph.adjacency[carbonyl_c]
    ):
        return None

    excluded = {n_idx}
    excluded.update(nb for nb in graph.adjacency[n_idx]
                    if nb != carbonyl_c and get_atom(graph, nb).symbol == "C")
    acid_chain = _collect_acid_chain(graph, carbonyl_c, excluded, get_atom)
    stem = CHAIN_PREFIX.get(len(acid_chain), f"C{len(acid_chain)}")
    _ene_te, _yne_te = _chain_multiple_bonds(graph, acid_chain)
    if _ene_te or _yne_te:
        from .name_assembler import _format_multiple_bonds as _fmt_te
        parent_name = f"{stem}{_fmt_te(_ene_te, _yne_te)}eteluramide"
    else:
        parent_name = f"{stem}aneteluramide"

    stereo_pfx_te = ""
    if _ene_te:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_te = PrincipalChain(atom_indices=acid_chain,
                                locant_map={c: i + 1 for i, c in enumerate(acid_chain)})
        _stereo_te = assign_stereochemistry(graph, _pc_te)
        if _stereo_te:
            stereo_pfx_te = "(" + ",".join(d.strip("()") for d in _stereo_te) + ")-"

    c_nbrs = [nb for nb in graph.adjacency[n_idx]
               if nb != carbonyl_c and get_atom(graph, nb).symbol == "C"]
    if not c_nbrs:
        return f"{stereo_pfx_te}{parent_name}"

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in c_nbrs]
    sub_counts = Counter(n_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if sub.startswith("(") else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub_str}")
    prefix = "-".join(prefix_parts)
    return f"{stereo_pfx_te}{prefix}{parent_name}"


def _name_anhydride(graph, pgrp, get_atom) -> str:
    """
    酸無水物命名: {acid1} {acid2} anhydride
    例: CC(=O)OC(=O)C → ethanoic anhydride
        CC(=O)OC(=O)CCC → butanoic ethanoic anhydride
    """
    from .constants import CHAIN_PREFIX

    c1 = pgrp.atom_indices[0]
    o_link = pgrp.atom_indices[1] if len(pgrp.atom_indices) > 1 else None
    c2 = pgrp.atom_indices[2] if len(pgrp.atom_indices) > 2 else None

    if c2 is None and o_link is not None:
        for nb in graph.adjacency[o_link]:
            if nb != c1 and get_atom(graph, nb).symbol == "C":
                c2 = nb
                break

    if c2 is None:
        return "anhydride"

    chain1 = _collect_acid_chain(graph, c1, set(), get_atom)
    chain2 = _collect_acid_chain(graph, c2, set(), get_atom)

    from .name_assembler import _format_multiple_bonds as _fmt_anh

    _ACID_RETAINED = {1: "formic", 2: "acetic"}

    def _aryl_acid_name(carbonyl_c: int) -> str | None:
        """Carbonyl C が芳香族環に直結している場合の酸名を返す。ベンゼン環 → 'benzoic'。"""
        for nb_idx in graph.adjacency[carbonyl_c]:
            nb = get_atom(graph, nb_idx)
            if nb.symbol == "C" and nb.is_aromatic and nb.in_ring:
                ring_atoms = next(
                    (set(rt) for rt in (graph.ring_atom_sets or []) if nb_idx in rt), set()
                )
                if (len(ring_atoms) == 6
                        and all(get_atom(graph, a).symbol == "C" for a in ring_atoms)):
                    return "benzoic"
                # Heteroaromatic
                if (any(get_atom(graph, a).symbol != "C" for a in ring_atoms)
                        and all(get_atom(graph, a).is_aromatic for a in ring_atoms)):
                    _apfx_anh = _aryl_sulfonyl_prefix(graph, nb_idx, carbonyl_c, get_atom)
                    if _apfx_anh is not None:
                        return f"{_apfx_anh}carboxylic"
        return None

    def _acid_stem_name(carbonyl_c: int, chain: list[int]) -> str:
        aryl = _aryl_acid_name(carbonyl_c)
        if aryl is not None:
            return aryl
        _ene, _yne = _chain_multiple_bonds(graph, chain)
        if not _ene and not _yne and len(chain) in _ACID_RETAINED:
            return _ACID_RETAINED[len(chain)]
        stem = CHAIN_PREFIX.get(len(chain), f"C{len(chain)}")
        if _ene or _yne:
            acid_name = f"{stem}{_fmt_anh(_ene, _yne)}oic"
            if _ene:
                from .stereochemistry import assign_stereochemistry
                from .chain_finder import PrincipalChain
                _pc = PrincipalChain(atom_indices=chain,
                                     locant_map={c: i + 1 for i, c in enumerate(chain)})
                _stereo = assign_stereochemistry(graph, _pc)
                if _stereo:
                    _combined = ",".join(d.strip("()") for d in _stereo)
                    return f"({_combined})-{acid_name}"
            return acid_name
        return f"{stem}anoic"

    acid1 = _acid_stem_name(c1, chain1)
    acid2 = _acid_stem_name(c2, chain2)

    acids = sorted([acid1, acid2])
    if acids[0] == acids[1]:
        return f"{acids[0]} anhydride"
    return f"{acids[0]} {acids[1]} anhydride"



def _alkylidene_name(graph, imine_c: int, n_idx: int, get_atom) -> str:
    """
    C=N の C 側 (imine_c) からアルキリデン名を返す (Phase 118)。
    例: imine_c + CH3 chain → "ethylidene"
        imine_c + Ph → "benzylidene"
    """
    from .constants import CHAIN_PREFIX

    c_nbs = [nb for nb in graph.adjacency[imine_c]
             if nb != n_idx and get_atom(graph, nb).symbol == "C"]

    ring_cs = [c for c in c_nbs if get_atom(graph, c).in_ring]
    chain_cs = [c for c in c_nbs if not get_atom(graph, c).in_ring]

    # アルキル鎖の長さを DFS で収集
    chain = [imine_c]
    visited: set[int] = {imine_c, n_idx}
    visited.update(ring_cs)
    q = list(chain_cs)
    for c in q:
        if c not in visited:
            visited.add(c)
            chain.append(c)
            for nb in graph.adjacency[c]:
                if (nb not in visited
                        and get_atom(graph, nb).symbol == "C"
                        and not get_atom(graph, nb).in_ring):
                    q.append(nb)

    stem = CHAIN_PREFIX.get(len(chain), f"C{len(chain)}")

    if ring_cs:
        ring_c = ring_cs[0]
        ring_set = next((rt for rt in (graph.ring_atom_sets or []) if ring_c in rt), None)
        if (ring_set and len(ring_set) == 6
                and all(get_atom(graph, a).symbol == "C" and get_atom(graph, a).is_aromatic
                        for a in ring_set)):
            if len(chain) == 1:
                return "benzylidene"
            return f"phenyl{stem}ylidene"

    return f"{stem}ylidene"


def _name_n_substituted_imine(
    graph, imine_c: int, n_idx: int, n_c_subs: list[int],
    find_principal_chain, get_multiple_bond_locants,
    collect_substituents, assemble_name, get_atom,
) -> str:
    """N-置換イミン: N-{sub}propan-2-imine 形式 (Phase 64)"""
    from .constants import MULTIPLIER
    from .substituent import _name_carbon_substituent
    from .functional_group import FunctionalGroup, FUNCTIONAL_GROUP_PRIORITY
    from collections import Counter

    # Phase 118: N が芳香環に直結 → N-alkylidene{arylamine} 形式 (Schiff base)
    if len(n_c_subs) == 1:
        n_c = n_c_subs[0]
        n_c_atom = get_atom(graph, n_c)
        if n_c_atom.in_ring and n_c_atom.is_aromatic:
            # アリールアミン親名を特定
            ring_set = next((rt for rt in (graph.ring_atom_sets or []) if n_c in rt), None)
            if ring_set is not None:
                ring_has_hetero = any(get_atom(graph, a).symbol != "C" for a in ring_set)
                ring_size = len(ring_set)
                if ring_size == 6 and not ring_has_hetero:
                    # ベンゼン環 → aniline / N-sub-aniline 系
                    # 環上の置換基を確認
                    from .ring_handler import _assign_ring_locants, collect_ring_substituents, assemble_ring_name
                    ring_chain = _assign_ring_locants(graph, list(ring_set), True, "alkane", [n_idx])
                    ring_subs = collect_ring_substituents(graph, ring_chain, [n_idx])
                    if ring_subs:
                        ring_base = assemble_ring_name(ring_chain, ring_subs, "alkane", None, [])
                        arylamine = f"{ring_base}aniline"
                    else:
                        arylamine = "aniline"
                    alk = _alkylidene_name(graph, imine_c, n_idx, get_atom)
                    return f"N-{alk}{arylamine}"

    pseudo_pgrp = FunctionalGroup(
        group_type="imine",
        atom_indices=[imine_c, n_idx],
        priority=FUNCTIONAL_GROUP_PRIORITY["imine"],
    )
    chain = find_principal_chain(graph, pseudo_pgrp)
    mb = get_multiple_bond_locants(graph, chain)
    suffix_locant = chain.locant_map.get(imine_c)
    substituents = collect_substituents(
        graph, chain.atom_indices, chain.locant_map, [imine_c, n_idx]
    )

    # Phase 349: C=N 結合の E/Z
    from .stereochemistry import _get_bond_stereo as _gbs_im
    _cn_st_im = _gbs_im(graph, imine_c, n_idx)
    _cn_pfx_im = f"({_cn_st_im})-" if _cn_st_im is not None else ""

    imine_base = assemble_name(
        chain_length=chain.length,
        principal_group_type="imine",
        multiple_bonds=mb,
        substituents=substituents,
        stereo_descriptors=[],
        suffix_locant=suffix_locant,
    )

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in n_c_subs]
    sub_counts = Counter(n_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if sub.startswith("(") else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub_str}")
    n_prefix = "-".join(prefix_parts)
    return f"{_cn_pfx_im}{n_prefix}{imine_base}"


def _name_secondary_tertiary_amide(graph, carbonyl_c: int, n_idx: int, get_atom) -> str:
    """
    二級・三級アミドの命名: N-alkyl-alkanamide 形式。
    例: CC(=O)NC → N-methylethanamide
        CC(=O)N(C)C → N,N-dimethylethanamide
        CC(C)C(=O)NC → 2-methyl-N-methylpropanamide
    """
    from .constants import CHAIN_PREFIX, MULTIPLIER
    from .substituent import _name_carbon_substituent
    from collections import Counter

    # 酸側炭素鎖: carbonyl_c から N 方向を除いた DFS
    excluded = {n_idx}
    excluded.update(nb for nb in graph.adjacency[n_idx]
                    if nb != carbonyl_c and get_atom(graph, nb).symbol == "C")
    acid_carbons = _collect_acid_chain(graph, carbonyl_c, excluded, get_atom)
    stem = CHAIN_PREFIX.get(len(acid_carbons), f"C{len(acid_carbons)}")

    # 鎖内の多重結合を検出して ene/yne 接尾辞を付加
    n_acid = len(acid_carbons)
    _ene_a, _yne_a = _chain_multiple_bonds(graph, acid_carbons)
    if _ene_a or _yne_a:
        from .name_assembler import _format_multiple_bonds as _fmt_a
        parent_name = f"{stem}{_fmt_a(_ene_a, _yne_a)}amide"
    elif n_acid == 1:
        parent_name = "formamide"
    elif n_acid == 2:
        parent_name = "acetamide"
    else:
        parent_name = f"{stem}anamide"

    # E/Z stereo prefix
    stereo_prefix_sta = ""
    if _ene_a:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_sta = PrincipalChain(atom_indices=acid_carbons,
                                 locant_map={c: i + 1 for i, c in enumerate(acid_carbons)})
        _stereo_sta = assign_stereochemistry(graph, _pc_sta)
        if _stereo_sta:
            _comb_sta = ",".join(d.strip("()") for d in _stereo_sta)
            stereo_prefix_sta = f"({_comb_sta})-"

    # 鎖上の置換基 (e.g. 2-methyl branch on acid chain)
    from .substituent import collect_substituents as _cs_amid
    _o_atoms_amid = [nb for nb in graph.adjacency[carbonyl_c]
                     if get_atom(graph, nb).symbol == "O"]
    _chain_lmap_amid = {c: i + 1 for i, c in enumerate(acid_carbons)}
    _chain_subs_amid = _cs_amid(graph, acid_carbons, _chain_lmap_amid,
                                 [n_idx] + _o_atoms_amid)
    chain_sub_parts: list[str] = []
    for _loc_am, _snm_am in sorted(_chain_subs_amid):
        chain_sub_parts.append(f"{_loc_am}-{_snm_am}")

    # N 置換基: N の隣接原子（carbonyl_c を除く）
    c_nbrs = [nb for nb in graph.adjacency[n_idx]
               if nb != carbonyl_c and get_atom(graph, nb).symbol == "C"]
    # N-hydroxy 置換基 (Phase 61: hydroxamic acid)
    o_nbrs = [nb for nb in graph.adjacency[n_idx]
              if get_atom(graph, nb).symbol == "O"]
    n_hydroxy_subs = []
    for o_nb in o_nbrs:
        has_h = any(get_atom(graph, nb2).symbol == "H"
                    for nb2 in graph.adjacency[o_nb])
        if has_h:
            n_hydroxy_subs.append("hydroxy")

    if not c_nbrs and not n_hydroxy_subs and not chain_sub_parts:
        return f"{stereo_prefix_sta}{parent_name}"

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in c_nbrs]
    n_subs += n_hydroxy_subs
    sub_counts = Counter(n_subs)
    n_prefix_parts: list[str] = []
    from .name_assembler import _needs_bis_tris as _nbp
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if _nbp(sub) else sub
        if cnt == 1:
            n_prefix_parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            n_prefix_parts.append(f"N,N-{mult}{sub_str}")

    if not n_prefix_parts and not chain_sub_parts:
        return f"{stereo_prefix_sta}{parent_name}"

    def _alpha_key_amide(part: str) -> tuple:
        import re as _re
        is_letter = bool(_re.match(r'^N[,N]*-', part))
        s = _re.sub(r'^[\d,]+-', '', part)
        s = _re.sub(r'^N[,N]*-', '', s)
        s = _re.sub(r'^\(', '', s)
        s = _re.sub(r'^(di|tri|tetra|penta|hexa|hepta|octa|nona|deca|bis|tris)', '', s)
        return (s.lower(), is_letter)

    all_parts = sorted(chain_sub_parts + n_prefix_parts, key=_alpha_key_amide)
    prefix = "-".join(all_parts)
    return f"{stereo_prefix_sta}{prefix}{parent_name}"


def _name_substituted_amidine(graph, amidine_c: int,
                               n_imine: int | None, n_amine: int | None,
                               get_atom) -> str:
    """
    アミジン命名 (N-置換含む): N-alkyl-alkan-imidamide 形式 (Phase 185)
    例: CC(=N)NC  → N-methylethanimidamide  (N = amine N, single bond)
        CC(=NC)N  → N'-methylethanimidamide (N' = imine N, double bond)
    """
    from .constants import CHAIN_PREFIX, MULTIPLIER
    from .substituent import _name_carbon_substituent
    from collections import Counter
    from .name_assembler import _needs_bis_tris as _nbp

    # グアニジン: アミジン炭素に N が 3 つ付いている場合 (NC(=N)N)
    n_all = [nb for nb in graph.adjacency[amidine_c] if get_atom(graph, nb).symbol == "N"]
    if len(n_all) >= 3:
        # N-置換グアニジン: N 上のアルキル置換基を収集する
        # グアニジンの命名規則:
        #   imine N (=N, double bond) → N'' (not used in IUPAC, typically unprimed for C=NH)
        #   amine N 1 (single bond) → N
        #   amine N 2 (single bond) → N' (the second amine N)
        # Note: IUPAC convention uses N and N' for the two amine N's in guanidine
        from .functional_group import get_bond_order
        from .constants import MULTIPLIER

        imine_n = [n for n in n_all if get_bond_order(graph, amidine_c, n) == 2.0]
        amine_ns = [n for n in n_all if get_bond_order(graph, amidine_c, n) == 1.0]

        def _get_n_subs(n_idx: int) -> list[str]:
            c_nbrs = [nb for nb in graph.adjacency[n_idx]
                      if nb != amidine_c and get_atom(graph, nb).symbol == "C"]
            return [_name_carbon_substituent(graph, c, {n_idx}) for c in c_nbrs]

        sub_parts: list[str] = []

        # imine N substituents (=N, labeled as the C=N nitrogen)
        for n_idx in imine_n:
            subs = _get_n_subs(n_idx)
            for alky in sorted(set(subs)):
                cnt = subs.count(alky)
                if cnt == 1:
                    sub_parts.append(("N''", alky))
                else:
                    mult = MULTIPLIER.get(cnt, f"{cnt}")
                    sub_parts.append((",".join(["N''"] * cnt), f"{mult}{alky}"))

        # amine N substituents: more-substituted N first → "N", other → "N'"
        # Sort so the N with more substituents (or heavier by sorted name) comes first,
        # making naming independent of SMILES parsing order.
        amine_ns_sorted = sorted(
            amine_ns,
            key=lambda n: (-len(_get_n_subs(n)), sorted(_get_n_subs(n))),
        )
        amine_labels = ["N", "N'"]
        for i, n_idx in enumerate(amine_ns_sorted):
            label = amine_labels[i] if i < len(amine_labels) else f"N{i+1}"
            subs = _get_n_subs(n_idx)
            for alky in sorted(set(subs)):
                cnt = subs.count(alky)
                if cnt == 1:
                    sub_parts.append((label, alky))
                else:
                    mult = MULTIPLIER.get(cnt, f"{cnt}")
                    labels_str = ",".join([label] * cnt)
                    sub_parts.append((labels_str, f"{mult}{alky}"))

        if not sub_parts:
            return "guanidine"

        # Merge identical alkyl groups on different N/N' into N,N'-di... form
        # Group by alkyl name, collect N-labels
        from collections import defaultdict
        by_alkyl: dict[str, list[str]] = defaultdict(list)
        for label, alky in sub_parts:
            by_alkyl[alky].append(label)

        result_parts: list[str] = []
        for alky in sorted(by_alkyl):
            labels = sorted(by_alkyl[alky])
            alky_str = f"({alky})" if alky.startswith("(") else alky
            if len(labels) == 1:
                result_parts.append(f"{labels[0]}-{alky_str}")
            else:
                labels_str = ",".join(labels)
                mult = MULTIPLIER.get(len(labels), f"{len(labels)}")
                result_parts.append(f"{labels_str}-{mult}{alky_str}")

        prefix = "-".join(result_parts)
        return f"{prefix}guanidine"

    excl: set[int] = set()
    if n_imine is not None:
        excl.add(n_imine)
    if n_amine is not None:
        excl.add(n_amine)
    acid_carbons = _collect_acid_chain(graph, amidine_c, excl, get_atom)
    n_acid = len(acid_carbons)
    stem = CHAIN_PREFIX.get(n_acid, f"C{n_acid}")
    _ene_ai, _yne_ai = _chain_multiple_bonds(graph, acid_carbons)
    if _ene_ai or _yne_ai:
        from .name_assembler import _format_multiple_bonds as _fmt_ai
        parent_name = f"{stem}{_fmt_ai(_ene_ai, _yne_ai)}imidamide"
    else:
        parent_name = f"{stem}animidamide"

    stereo_pfx_ai = ""
    if _ene_ai:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_ai = PrincipalChain(atom_indices=acid_carbons,
                                locant_map={c: i + 1 for i, c in enumerate(acid_carbons)})
        _stereo_ai = assign_stereochemistry(graph, _pc_ai)
        if _stereo_ai:
            stereo_pfx_ai = "(" + ",".join(d.strip("()") for d in _stereo_ai) + ")-"

    # N-置換基 (amine N, single-bond): "N-" prefix
    n_subs: list[str] = []
    n_hydroxy = False  # N-OH (amidoxime pattern)
    if n_amine is not None:
        c_nbrs = [nb for nb in graph.adjacency[n_amine]
                  if nb != amidine_c and get_atom(graph, nb).symbol == "C"]
        n_subs = [_name_carbon_substituent(graph, c, {n_amine}) for c in c_nbrs]
        # アミドキシム: amine N に O-H が付いている場合 (N-hydroxy)
        from .functional_group import get_bond_order as _gbo
        for o_nb in graph.adjacency[n_amine]:
            if get_atom(graph, o_nb).symbol == "O" and _gbo(graph, n_amine, o_nb) == 1.0:
                if any(get_atom(graph, hh).symbol == "H" for hh in graph.adjacency[o_nb]):
                    n_hydroxy = True
                    break

    # N'-置換基 (imine N, double-bond): "N'-" prefix
    np_subs: list[str] = []
    np_hydroxy = False  # N'-OH pattern
    if n_imine is not None:
        c_nbrs_i = [nb for nb in graph.adjacency[n_imine]
                    if nb != amidine_c and get_atom(graph, nb).symbol == "C"]
        np_subs = [_name_carbon_substituent(graph, c, {n_imine}) for c in c_nbrs_i]
        # N'-OH on imine N
        from .functional_group import get_bond_order as _gbo2
        for o_nb in graph.adjacency[n_imine]:
            if get_atom(graph, o_nb).symbol == "O" and _gbo2(graph, n_imine, o_nb) == 1.0:
                if any(get_atom(graph, hh).symbol == "H" for hh in graph.adjacency[o_nb]):
                    np_hydroxy = True
                    break

    if not n_subs and not np_subs and not n_hydroxy and not np_hydroxy:
        return f"{stereo_pfx_ai}{parent_name}"

    prefix_parts: list[str] = []

    # N-hydroxy (amidoxime): add "N-hydroxy" prefix
    if n_hydroxy:
        prefix_parts.append("N-hydroxy")

    sub_counts = Counter(n_subs)
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if _nbp(sub) else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub_str}")

    # N'-hydroxy
    if np_hydroxy:
        prefix_parts.append("N'-hydroxy")

    np_counts = Counter(np_subs)
    for sub in sorted(np_counts):
        cnt = np_counts[sub]
        sub_str = f"({sub})" if _nbp(sub) else sub
        if cnt == 1:
            prefix_parts.append(f"N'-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N',N'-{mult}{sub_str}")

    prefix = "-".join(prefix_parts)
    return f"{stereo_pfx_ai}{prefix}{parent_name}"


def _name_secondary_tertiary_amine(graph, n_idx: int, c_neighbors: list[int], get_atom) -> str:
    """
    二級・三級アミンの命名: N-alkyl-alkanamine 形式。
    例: CNC → N-methylmethanamine, CN(C)C → N,N-dimethylmethanamine
        CNc1ccccc1 → N-methylaniline (アリールアミン)
        CNC1CCCCC1 → N-methylcyclohexanamine (シクロアルキルアミン)
    """
    from .constants import CHAIN_PREFIX, MULTIPLIER
    from .substituent import _name_carbon_substituent
    from collections import Counter

    # Phase 36: 環炭素を隣接する場合は環を親とする
    ring_c_nbrs = [c for c in c_neighbors if get_atom(graph, c).in_ring]
    non_ring_c_nbrs = [c for c in c_neighbors if not get_atom(graph, c).in_ring]

    # Phase 112: ジフェニルアミン / トリフェニルアミン保留名 (IUPAC P-62.2.3)
    if ring_c_nbrs and not non_ring_c_nbrs:
        def _is_plain_benzene(c_idx: int) -> bool:
            ra = next((rt for rt in (graph.ring_atom_sets or []) if c_idx in rt), None)
            return (ra is not None and len(ra) == 6
                    and all(get_atom(graph, a).symbol == "C" and get_atom(graph, a).is_aromatic
                            for a in ra))
        if all(_is_plain_benzene(c) for c in ring_c_nbrs):
            _n_ph = len(ring_c_nbrs)
            if _n_ph == 2:
                return "diphenylamine"
            if _n_ph == 3:
                return "triphenylamine"

    if ring_c_nbrs:
        ring_c = ring_c_nbrs[0]
        rc_atom = get_atom(graph, ring_c)

        # 環サイズを特定
        ring_size = 0
        for ring_set in (graph.ring_atom_sets or []):
            if ring_c in ring_set:
                ring_size = len(ring_set)
                break

        # ヘテロ芳香族チェック
        ring_has_hetero = any(
            get_atom(graph, a).symbol != "C"
            for ring_set in (graph.ring_atom_sets or [])
            if ring_c in ring_set
            for a in ring_set
        )
        if rc_atom.is_aromatic and ring_size == 6 and not ring_has_hetero:
            parent_name = "aniline"
        elif rc_atom.is_aromatic and ring_size == 6 and ring_has_hetero:
            # ヘテロ芳香族アミン: "{base}-{locant}-amine" (e.g. pyridin-4-amine)
            from .substituent import _name_aryl_substituent
            from .heterocycle_handler import _match_retained
            from .ring_handler import _order_ring
            ring_tuple = next(
                (rs for rs in (graph.ring_atom_sets or []) if ring_c in rs), None
            )
            if ring_tuple is not None:
                ring_ordered = _order_ring(list(ring_tuple), graph)
                match = _match_retained(ring_ordered, graph)
                if match is not None:
                    base_nm, _is_nh, rotation = match
                    locant = rotation.index(ring_c) + 1 if ring_c in rotation else 1
                    stem_nm = base_nm[:-1] if base_nm.endswith("e") else base_nm
                    parent_name = f"{stem_nm}-{locant}-amine"
                else:
                    parent_name = "amine"  # fallback
            else:
                parent_name = "amine"
        elif rc_atom.is_aromatic:
            # ナフタレン等: フォールバック (benzen-amine 系)
            stem = CHAIN_PREFIX.get(ring_size, f"C{ring_size}")
            parent_name = f"{stem}anamine"
        else:
            stem = CHAIN_PREFIX.get(ring_size, f"C{ring_size}")
            parent_name = f"cyclo{stem}anamine"

        # 非環炭素がN置換基
        n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in non_ring_c_nbrs]
        # 余分な環隣接炭素（複数環がある場合）
        for rc2 in ring_c_nbrs[1:]:
            rc2_atom = get_atom(graph, rc2)
            n_subs.append("phenyl" if rc2_atom.is_aromatic else
                          _name_carbon_substituent(graph, rc2, {n_idx}))
        # N-ハロゲン置換基 (N-chloro 等)
        _HAL_NAMES_N2 = {"F": "fluoro", "Cl": "chloro", "Br": "bromo", "I": "iodo"}
        for _nb_idx2 in graph.adjacency[n_idx]:
            _nb_sym2 = get_atom(graph, _nb_idx2).symbol
            if _nb_sym2 in _HAL_NAMES_N2:
                n_subs.append(_HAL_NAMES_N2[_nb_sym2])

        if not n_subs:
            return parent_name

        sub_counts = Counter(n_subs)
        prefix_parts = []
        for sub in sorted(sub_counts):
            cnt = sub_counts[sub]
            needs_p = sub.startswith("(") or any(c.isdigit() for c in sub)
            sub_str = f"({sub})" if needs_p else sub
            if cnt == 1:
                prefix_parts.append(f"N-{sub_str}")
            else:
                mult = MULTIPLIER.get(cnt, f"{cnt}")
                prefix_parts.append(f"N,N-{mult}{sub_str}")
        prefix = "-".join(prefix_parts)
        return f"{prefix}{parent_name}"

    # 非環アミン: chain_through_pivot で N を内点として含む最長鎖を選ぶ
    from .chain_finder import chain_through_pivot as _ctp_am

    def chain_through_c(start_c: int) -> tuple[list[int], int]:
        if get_atom(graph, start_c).in_ring:
            return [], 1
        return _ctp_am(graph, start_c, {n_idx}, get_atom)

    chain_locants = [chain_through_c(c) for c in c_neighbors]
    chains = [cl[0] for cl in chain_locants]
    locants_n = [cl[1] for cl in chain_locants]
    lengths = [len(ch) for ch in chains]
    max_len = max(lengths) if lengths else 1

    # 最長鎖を親鎖とする（最初に見つかったものを選択）
    parent_pos = next((i for i, l in enumerate(lengths) if l == max_len), 0)
    parent_chain = chains[parent_pos]
    n_locant = locants_n[parent_pos]
    stem = CHAIN_PREFIX.get(max_len, f"C{max_len}")

    # 多重結合の検出: ene/yne がある場合は親鎖の locant を明示
    _ene_am, _yne_am = _chain_multiple_bonds(graph, parent_chain)
    if _ene_am or _yne_am:
        from .name_assembler import _format_multiple_bonds as _fmt_am
        if max_len == 2 and _ene_am == [1] and not _yne_am:
            parent_name = f"{stem}enamine"  # ethenamine: locants unambiguous for 2-carbon
        else:
            parent_name = f"{stem}{_fmt_am(_ene_am, _yne_am)}-{n_locant}-amine"
    elif n_locant > 1:
        parent_name = f"{stem}an-{n_locant}-amine"
    elif max_len >= 3:
        parent_name = f"{stem}an-1-amine"
    else:
        parent_name = f"{stem}anamine"

    # E/Z stereo for parent chain double bonds
    stereo_pfx_am = ""
    if _ene_am:
        from .stereochemistry import assign_stereochemistry
        from .chain_finder import PrincipalChain
        _pc_am = PrincipalChain(atom_indices=parent_chain,
                                locant_map={c: i + 1 for i, c in enumerate(parent_chain)})
        _stereo_am = assign_stereochemistry(graph, _pc_am)
        if _stereo_am:
            stereo_pfx_am = "(" + ",".join(d.strip("()") for d in _stereo_am) + ")-"

    # 親鎖上の置換基（分岐鎖 etc.）
    from .substituent import collect_substituents as _cs_am
    _chain_locant_map = {c: i + 1 for i, c in enumerate(parent_chain)}
    _pgrp_excluded = {n_idx}
    _chain_subs = _cs_am(graph, parent_chain, _chain_locant_map, list(_pgrp_excluded))
    chain_sub_parts: list[str] = []
    for _loc, _sname in sorted(_chain_subs):
        _needs_p2 = _sname.startswith("(") or any(ch.isdigit() for ch in _sname)
        _sstr2 = f"({_sname})" if _needs_p2 else _sname
        chain_sub_parts.append(f"{_loc}-{_sstr2}")

    # N 上の置換基（親鎖以外）
    n_subs = [
        _name_carbon_substituent(graph, c, {n_idx})
        for i, c in enumerate(c_neighbors) if i != parent_pos
    ]

    # Phase 199: N-ハロゲン置換基 (N-chloro 等)
    _HAL_NAMES_N = {"F": "fluoro", "Cl": "chloro", "Br": "bromo", "I": "iodo"}
    for nb_idx in graph.adjacency[n_idx]:
        nb_sym = get_atom(graph, nb_idx).symbol
        if nb_sym in _HAL_NAMES_N:
            n_subs.append(_HAL_NAMES_N[nb_sym])

    if not n_subs and not chain_sub_parts:
        return f"{stereo_pfx_am}{parent_name}"

    sub_counts = Counter(n_subs)
    n_prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        needs_p = sub.startswith("(") or any(c.isdigit() for c in sub)
        sub_str = f"({sub})" if needs_p else sub
        if cnt == 1:
            n_prefix_parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            n_prefix_parts.append(f"N,N-{mult}{sub_str}")

    def _alpha_key(part: str) -> tuple:
        import re
        is_letter = bool(re.match(r'^N[,N]*-', part))
        s = re.sub(r'^[\d,]+-', '', part)
        s = re.sub(r'^N[,N]*-', '', s)
        s = re.sub(r'^\(', '', s)
        s = re.sub(r'^(di|tri|tetra|penta|hexa|hepta|octa|nona|deca|bis|tris)', '', s)
        return (s.lower(), is_letter)

    all_parts = sorted(chain_sub_parts + n_prefix_parts, key=_alpha_key)
    prefix = "-".join(all_parts)
    return f"{stereo_pfx_am}{prefix}{parent_name}"


# ─── Ring-conditional dispatch wrappers ───────────────────────────────────────

def _dispatch_diester(graph, pgrp, get_atom):
    """Returns None when ring-cyclic path should be taken."""
    d_cs = [ai for ai in pgrp.atom_indices
            if get_atom(graph, ai).symbol == "C" and get_atom(graph, ai).in_ring]
    d_ring_os = [ai for ai in pgrp.atom_indices
                 if get_atom(graph, ai).symbol == "O" and get_atom(graph, ai).in_ring]
    if not (d_cs and d_ring_os):
        return _name_diester(graph, pgrp, get_atom)
    return None


def _dispatch_ester(graph, pgrp, get_atom):
    ester_c = pgrp.atom_indices[0]
    ring_o_in_ester = [ai for ai in pgrp.atom_indices[1:] if get_atom(graph, ai).in_ring]
    if not (get_atom(graph, ester_c).in_ring and ring_o_in_ester):
        return _name_ester(graph, pgrp, get_atom)
    return None


def _dispatch_carbamate(graph, pgrp, get_atom):
    if not get_atom(graph, pgrp.atom_indices[0]).in_ring:
        return _name_carbamate(graph, pgrp, get_atom)
    return None


def _name_amidine_pgrp(graph, pgrp, get_atom) -> str | None:
    """PGRP_DISPATCH entry for amidine: handle heteroaromatic ring case only.

    For benzene/aliphatic amidines, return None so the ring_handler
    (spec.benzene_name) or _name_acyclic path handles it correctly.
    """
    from .functional_group import _get_amidine_nitrogens
    from .constants import MULTIPLIER
    from .substituent import _name_carbon_substituent
    from collections import Counter

    amidine_c = pgrp.atom_indices[0]
    if get_atom(graph, amidine_c).in_ring:
        return None

    for _nb_am in graph.adjacency[amidine_c]:
        _nba_am = get_atom(graph, _nb_am)
        if _nba_am.symbol != "C" or not _nba_am.in_ring or not _nba_am.is_aromatic:
            continue
        _ra_am = next(
            (rt for rt in (graph.ring_atom_sets or []) if _nb_am in rt), None
        )
        if _ra_am is None:
            continue
        if all(get_atom(graph, a).symbol == "C" for a in _ra_am):
            return None  # benzene → let ring_handler use spec.benzene_name
        if (any(get_atom(graph, a).symbol != "C" for a in _ra_am)
                and all(get_atom(graph, a).is_aromatic for a in _ra_am)):
            _apfx_am = _aryl_sulfonyl_prefix(graph, _nb_am, amidine_c, get_atom)
            if _apfx_am is None:
                return None
            _n_imine, _n_amine = _get_amidine_nitrogens(graph, amidine_c)
            _het_base_am = f"{_apfx_am}carboximidamide"
            # N-amine substituents
            _pfx_parts: list[str] = []
            if _n_amine is not None:
                _cn_am = [nb for nb in graph.adjacency[_n_amine]
                           if nb != amidine_c and get_atom(graph, nb).symbol == "C"]
                if _cn_am:
                    _ns = [_name_carbon_substituent(graph, c, {_n_amine}) for c in _cn_am]
                    _sc = Counter(_ns)
                    for _s in sorted(_sc):
                        _c = _sc[_s]
                        _ss = f"({_s})" if _s.startswith("(") else _s
                        if _c == 1:
                            _pfx_parts.append(f"N-{_ss}")
                        else:
                            _pfx_parts.append(f"N,N-{MULTIPLIER.get(_c, str(_c))}{_ss}")
            # N'-imine substituents
            if _n_imine is not None:
                _cn_im = [nb for nb in graph.adjacency[_n_imine]
                           if nb != amidine_c and get_atom(graph, nb).symbol == "C"]
                if _cn_im:
                    _ns2 = [_name_carbon_substituent(graph, c, {_n_imine}) for c in _cn_im]
                    _sc2 = Counter(_ns2)
                    for _s in sorted(_sc2):
                        _c = _sc2[_s]
                        _ss = f"({_s})" if _s.startswith("(") else _s
                        if _c == 1:
                            _pfx_parts.append(f"N'-{_ss}")
                        else:
                            _pfx_parts.append(f"N',N'-{MULTIPLIER.get(_c, str(_c))}{_ss}")
            pfx = "-".join(_pfx_parts)
            return f"{pfx}{_het_base_am}"
    return None


PGRP_DISPATCH: dict = {
    "imidic_acid": _name_imidic_acid,
    "imidate_ester": _name_imidate_ester,
    "diester": _dispatch_diester,
    "ester": _dispatch_ester,
    "diacid_halide": _name_diacid_halide,
    "acid_halide": _name_acid_halide,
    "anhydride": _name_anhydride,
    "carbonate": _name_carbonate,
    "carbonothioate": _name_carbonothioate,
    "carbonodithioate": _name_carbonodithioate,
    "o_thiocarbamate": _name_o_thiocarbamate,
    "s_carbamothioate": _name_s_carbamothioate,
    "s_carbamodithioate": _name_s_carbamodithioate,
    "sulfonyl_chloride": _name_sulfonyl_chloride,
    "sulfinyl_chloride": _name_sulfinyl_chloride,
    "chloroformate": _name_chloroformate,
    "sulfonic_acid": _name_sulfonic_acid,
    "sulfonate_anion": _name_sulfonate_anion,
    "disulfonic_acid": _name_disulfonic_acid,
    "sulfinic_acid": _name_sulfinic_acid,
    "sulfenic_acid": _name_sulfenic_acid,
    "sulfenyl_halide": _name_sulfenyl_halide,
    "sulfenate_ester": _name_sulfenate_ester,
    "sulfonate_ester": _name_sulfonate_sulfinate_ester,
    "sulfinate_ester": _name_sulfonate_sulfinate_ester,
    "sulfoxide": _name_sulfoxide_sulfone,
    "sulfone": _name_sulfoxide_sulfone,
    "selenoxide": _name_selenoxide_selenone,
    "selenone": _name_selenoxide_selenone,
    "sulfonyl_azide": _name_sulfonyl_azide,
    "sulfonohydrazide": _name_sulfonohydrazide,
    "sulfinylhydrazide": _name_sulfinylhydrazide,
    "sulfonamide": _name_sulfonamide,
    "sulfamic_acid": _name_sulfamic_acid,
    "disulfonamide": _name_disulfonamide,
    "sulfinamide": _name_sulfinamide,
    "sulfide": _name_sulfide,
    "phosphate_ester": _name_phosphate_ester,
    "phosphonate_halfester": _name_phosphonate_halfester,
    "phosphonate_ester": _name_phosphonate_ester,
    "phosphonothioate_ester": _name_phosphonothioate_ester,
    "phosphinate_ester": _name_phosphinate_ester,
    "phosphonic_acid": _name_phosphonic_acid,
    "phosphinic_acid": _name_phosphinic_acid,
    "phosphonous_acid": _name_phosphonous_acid,
    "phosphinous_acid": _name_phosphinous_acid,
    "arsonic_acid": _name_arsonic_acid,
    "arsinic_acid": _name_arsinic_acid,
    "arsonous_acid": _name_arsonous_acid,
    "arsinous_acid": _name_arsinous_acid,
    "phosphane": _name_phosphane,
    "phosphine_oxide": _name_phosphine_oxide,
    "phosphite_ester": _name_phosphite_ester,
    "boronic_acid": _name_boronic_acid,
    "boronate_ester": _name_boronate_ester,
    "borate_ester": _name_borate_ester,
    "borinic_acid": _name_borinic_acid,
    "borane_org": _name_organic_borane,
    "silane_org": _name_organic_silane,
    "silyl_ether_org": _name_silyl_ether,
    "disilane_org": _name_disilane,
    "disiloxane_org": _name_disiloxane,
    "disilazane_org": _name_disilazane,
    "germane_org": _name_organic_germane,
    "stannane_org": _name_organic_stannane,
    "arsane_org": _name_arsane_org,
    "organomercury": _name_organomercury,
    "bismuthane_org": _name_organic_bismuthane,
    "stibane_org": _name_organic_stibane,
    "plumbane_org": _name_organic_plumbane,
    "selenonic_acid": _name_selenonic_acid,
    "seleninic_acid": _name_seleninic_acid,
    "selenenic_acid": _name_selenenic_acid,
    "telluronic_acid": _name_telluronic_acid,
    "tellurinic_acid": _name_tellurinic_acid,
    "tellurenic_acid": _name_tellurenic_acid,
    "silanol_org": _name_organic_silanol,
    "isocyanide": _name_isocyanide,
    "ammonium": _name_ammonium,
    "phosphanium": _name_phosphanium,
    "sulfonium": _name_sulfonium,
    "arsonium": _name_arsonium,
    "carbamate": _dispatch_carbamate,
    "hydroperoxide": _name_hydroperoxide,
    "thioamide": _name_thioamide,
    "dithioamide": _name_dithioamide,
    "diselenoamide": _name_diselenoamide,
    "diamide": _name_diamide,
    "amidine": _name_amidine_pgrp,
    "diamidine": _name_diamidine,
    "dinitrile": _name_dinitrile,
    "selenoamide": _name_selenoamide,
    "telluramide": _name_telluramide,
    "diisocyanate": _name_diisocyanate,
    "diisothiocyanate": _name_diisothiocyanate,
    "diimine": _name_diimine,
    "dial": _name_dial,
    "peroxide": _name_peroxide,
    "thioester": _name_thioester,
    "o_thioester": _name_o_thioester,
    "s_dithioate_ester": _name_s_dithioate_ester,
    "disulfide": _name_disulfide,
    "trisulfide": _name_polysulfide,
    "tetrasulfide": _name_polysulfide,
    "selenide": _name_selenide_telluride,
    "telluride": _name_selenide_telluride,
    "diselenide": _name_diselenide_ditelluride,
    "ditelluride": _name_diselenide_ditelluride,
    "nitroso": _name_nitroso,
    "azide": _name_azide,
    "isocyanate": _name_isocyanate,
    "isothiocyanate": _name_isothiocyanate,
    "cyanate": _name_cyanate,
    "thiocyanate": _name_thiocyanate,
    "carboxylate": _name_carboxylate,
    "dicarboxylate": _name_dicarboxylate,
    "dioic_acid": _name_dioic_acid,
    "thioic_s_acid": _name_thioic_acid,
    "thioic_o_acid": _name_thioic_acid,
    "dithioic_acid": _name_thioic_acid,
    "nitrate_ester": _name_nitrate_ester,
    "nitrite_ester": _name_nitrite_ester,
    "carbamic_acid": _name_carbamic_acid,
    "carbodiimide": _name_carbodiimide,
    "acyl_azide": _name_acyl_azide,
    "hydrazide": _name_hydrazide,
    "thiohydrazide": _name_thiohydrazide,
    "selenohydrazide": _name_selenohydrazide,
    "peroxyacid": _name_peroxyacid,
    "peroxy_ester": _name_peroxyester,
    "semicarbazone": _name_semicarbazone,
    "aldsemicarbazone": _name_semicarbazone,
    "thiosemicarbazone": _name_semicarbazone,
    "aldthiosemicarbazone": _name_semicarbazone,
    "aldhydrazone": _name_substituted_hydrazone,
    "kethydrazone": _name_substituted_hydrazone,
}


def _name_diazonium(graph, get_atom) -> str | None:
    """ジアゾニウム塩の命名 (Phase 215, IUPAC 2013 P-68.3.2.4):

    R-[N+]#N → R-diazonium
    例:
        C[N+]#N            → methanediazonium
        [N+](#N)c1ccccc1   → benzenediazonium
    """
    from .functional_group import get_bond_order
    from .constants import CHAIN_PREFIX

    for idx in range(len(graph.atoms)):
        atom = get_atom(graph, idx)
        if atom.symbol != "N" or atom.formal_charge != 1:
            continue

        # N+ must have one triple bond to N (terminal)
        n_triple = [nb for nb in graph.adjacency[idx]
                    if get_atom(graph, nb).symbol == "N"
                    and get_bond_order(graph, idx, nb) == 3.0]
        if len(n_triple) != 1:
            continue

        # N+ must have one C neighbor (single bond)
        c_nbrs = [nb for nb in graph.adjacency[idx]
                  if get_atom(graph, nb).symbol == "C"
                  and get_bond_order(graph, idx, nb) == 1.0]
        if len(c_nbrs) != 1:
            continue

        c_idx = c_nbrs[0]
        c_atom = get_atom(graph, c_idx)

        # Benzene ring case: benzenediazonium or substituted arene-diazonium
        if c_atom.in_ring and c_atom.is_aromatic:
            # Collect all aromatic ring atoms (BFS)
            ring_members: list[int] = []
            visited: set[int] = set()
            stack = [c_idx]
            while stack:
                cur = stack.pop()
                if cur in visited:
                    continue
                visited.add(cur)
                cur_atom = get_atom(graph, cur)
                if not cur_atom.in_ring or cur_atom.symbol not in ("C", "N", "O", "S"):
                    continue
                ring_members.append(cur)
                for nb in graph.adjacency[cur]:
                    if nb not in visited and get_atom(graph, nb).in_ring:
                        stack.append(nb)

            # Only handle simple 6-membered carbocyclic ring (benzene)
            if len(ring_members) == 6 and all(get_atom(graph, a).symbol == "C" for a in ring_members):
                # Check for no other substituents on ring (simple benzenediazonium)
                ring_set = set(ring_members)
                non_ring_subs = [(rc, nb) for rc in ring_members
                                 for nb in graph.adjacency[rc]
                                 if nb not in ring_set
                                 and get_atom(graph, nb).symbol not in ("H",)
                                 and nb != idx]
                if not non_ring_subs:
                    return "benzenediazonium"
            # Fallthrough for complex aryl cases (not handled here)
            continue

        # Chain case: aliphatic diazonium
        chain = [c_idx]
        visited_c = {c_idx, idx}
        # Extend chain to find the alkyl chain
        current = c_idx
        while True:
            c_nbrs_chain = [nb for nb in graph.adjacency[current]
                           if get_atom(graph, nb).symbol == "C"
                           and nb not in visited_c
                           and not get_atom(graph, nb).in_ring]
            if not c_nbrs_chain:
                break
            # Pick longest chain (simplified: just follow first)
            nxt = c_nbrs_chain[0]
            chain.append(nxt)
            visited_c.add(nxt)
            current = nxt

        n_chain = len(chain)
        stem = CHAIN_PREFIX.get(n_chain)
        if stem is None:
            continue
        return f"{stem}anediazonium"

    return None


def _name_sulfate_ester(graph, get_atom) -> str | None:
    """硫酸エステル・スルファミン酸エステルの命名 (Phase 211, IUPAC 2013 P-67.2):

    COS(=O)(=O)O   → methyl hydrogen sulfate  (モノエステル)
    COS(=O)(=O)OC  → dimethyl sulfate          (ジエステル)
    NS(=O)(=O)OC   → methyl sulfamate          (スルファミン酸エステル)

    硫酸エステル: S が 2 つの =O と 2 つの -O- を持ち、かつ S に直結 C なし。
    """
    from .functional_group import get_bond_order
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER

    for idx in range(len(graph.atoms)):
        atom = get_atom(graph, idx)
        if atom.symbol != "S" or atom.in_ring:
            continue

        neighbors = graph.adjacency[idx]
        o_double = [nb for nb in neighbors
                    if get_atom(graph, nb).symbol == "O"
                    and get_bond_order(graph, idx, nb) == 2.0]
        o_single = [nb for nb in neighbors
                    if get_atom(graph, nb).symbol == "O"
                    and get_bond_order(graph, idx, nb) == 1.0]
        n_single = [nb for nb in neighbors
                    if get_atom(graph, nb).symbol == "N"
                    and get_bond_order(graph, idx, nb) == 1.0]
        c_direct = [nb for nb in neighbors if get_atom(graph, nb).symbol == "C"]

        # 硫酸エステル: 2 C=O (sulfate), no direct C on S
        if len(o_double) != 2 or c_direct:
            continue

        # スルファミン酸エステル: N-S(=O)₂-O-C
        if len(n_single) == 1 and len(o_single) == 1:
            # NS(=O)(=O)OC pattern: N must be NH2-like
            n_idx = n_single[0]
            n_h = [nb for nb in graph.adjacency[n_idx] if get_atom(graph, nb).symbol == "H"]
            n_c = [nb for nb in graph.adjacency[n_idx]
                   if nb != idx and get_atom(graph, nb).symbol == "C"]
            if not n_h and not n_c:
                pass  # no H, no C on N - unusual
            o_idx = o_single[0]
            oc_nbrs = [nb for nb in graph.adjacency[o_idx]
                       if nb != idx and get_atom(graph, nb).symbol == "C"]
            if oc_nbrs:
                alkyl = _name_carbon_substituent(graph, oc_nbrs[0], {o_idx})
                return f"{alkyl} sulfamate"
            continue

        # 硫酸エステル: 2-O single bonds, no direct C, no N
        if len(o_single) != 2 or n_single:
            continue

        # 各 O に付いている C を収集
        ester_cs: list[int] = []
        oh_count = 0
        for o_idx in o_single:
            c_nbrs = [nb for nb in graph.adjacency[o_idx]
                      if nb != idx and get_atom(graph, nb).symbol == "C"]
            h_nbrs = [nb for nb in graph.adjacency[o_idx]
                      if get_atom(graph, nb).symbol == "H"]
            if c_nbrs:
                ester_cs.append((o_idx, c_nbrs[0]))
            elif h_nbrs:
                oh_count += 1

        if len(ester_cs) == 0:
            continue  # neither ester nor monoester found

        alkyl_names = sorted(
            _name_carbon_substituent(graph, c, {o}) for o, c in ester_cs
        )

        if len(ester_cs) == 1:
            # モノエステル: alkyl hydrogen sulfate
            return f"{alkyl_names[0]} hydrogen sulfate"
        elif len(ester_cs) == 2:
            # ジエステル: dialkyl sulfate (or alkyl1 alkyl2 sulfate)
            if alkyl_names[0] == alkyl_names[1]:
                mult = MULTIPLIER.get(2, "di")
                return f"{mult}{alkyl_names[0]} sulfate"
            return f"{alkyl_names[0]} {alkyl_names[1]} sulfate"

    return None


def _name_sulfite_ester(graph, get_atom) -> str | None:
    """亜硫酸エステルの命名 (Phase 225, IUPAC 2013 P-67.2):

    COS(=O)OC  → dimethyl sulfite  (ジエステル)
    COS(=O)O   → methyl hydrogen sulfite  (モノエステル)

    亜硫酸エステル: S が 1 つの =O と 2 つの -O- を持ち、かつ S に直結 C なし。
    """
    from .functional_group import get_bond_order
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER

    for idx in range(len(graph.atoms)):
        atom = get_atom(graph, idx)
        if atom.symbol != "S" or atom.in_ring:
            continue
        neighbors = graph.adjacency[idx]
        o_double = [nb for nb in neighbors
                    if get_atom(graph, nb).symbol == "O"
                    and get_bond_order(graph, idx, nb) == 2.0]
        o_single = [nb for nb in neighbors
                    if get_atom(graph, nb).symbol == "O"
                    and get_bond_order(graph, idx, nb) == 1.0]
        c_direct = [nb for nb in neighbors if get_atom(graph, nb).symbol == "C"]
        n_direct = [nb for nb in neighbors if get_atom(graph, nb).symbol == "N"]
        halogen = [nb for nb in neighbors
                   if get_atom(graph, nb).symbol in ("Cl", "F", "Br", "I")]
        if len(o_double) != 1 or c_direct or n_direct or halogen:
            continue
        if len(o_single) != 2:
            continue
        ester_cs: list[tuple[int, int]] = []
        oh_count = 0
        for o_idx in o_single:
            c_nbrs = [nb for nb in graph.adjacency[o_idx]
                      if nb != idx and get_atom(graph, nb).symbol == "C"]
            h_nbrs = [nb for nb in graph.adjacency[o_idx]
                      if get_atom(graph, nb).symbol == "H"]
            if c_nbrs:
                ester_cs.append((o_idx, c_nbrs[0]))
            elif h_nbrs:
                oh_count += 1
        if not ester_cs:
            continue
        alkyl_names = sorted(
            _name_carbon_substituent(graph, c, {o}) for o, c in ester_cs
        )
        if len(ester_cs) == 1:
            return f"{alkyl_names[0]} hydrogen sulfite"
        if alkyl_names[0] == alkyl_names[1]:
            mult = MULTIPLIER.get(2, "di")
            return f"{mult}{alkyl_names[0]} sulfite"
        return f"{alkyl_names[0]} {alkyl_names[1]} sulfite"

    return None


def _name_acyl_peroxide(graph, get_atom) -> str | None:
    """アシルペルオキシドの命名 (Phase 209, IUPAC 2013 P-65.1.5.3):

    RC(=O)-O-O-C(=O)R' → di{acyl} peroxide (対称) / {acyl1} {acyl2} peroxide (非対称)

    例:
        CC(=O)OOC(=O)C  → diethanoyl peroxide
        CC(=O)OOC(=O)CC → ethanoyl propanoyl peroxide
    """
    from .functional_group import get_bond_order, _get_double_bonded_oxygen
    from .constants import CHAIN_PREFIX
    from .name_assembler import _format_multiple_bonds as _fmt_mb
    from .chain_finder import chain_multiple_bonds as _chain_mb

    # O-O 結合を探す
    for idx in range(len(graph.atoms)):
        atom = get_atom(graph, idx)
        if atom.symbol != "O":
            continue
        o1_idx = idx
        for o2_idx in graph.adjacency[o1_idx]:
            if get_atom(graph, o2_idx).symbol != "O":
                continue
            # O1-O2 が単結合でなければスキップ
            if get_bond_order(graph, o1_idx, o2_idx) != 1.0:
                continue
            # O1, O2 それぞれに隣接する C を探す
            c1_candidates = [nb for nb in graph.adjacency[o1_idx]
                             if nb != o2_idx and get_atom(graph, nb).symbol == "C"]
            c2_candidates = [nb for nb in graph.adjacency[o2_idx]
                             if nb != o1_idx and get_atom(graph, nb).symbol == "C"]
            if not c1_candidates or not c2_candidates:
                continue
            c1_idx = c1_candidates[0]
            c2_idx = c2_candidates[0]
            # 両 C に C=O があること (acyl groups)
            if (_get_double_bonded_oxygen(graph, c1_idx) is None
                    or _get_double_bonded_oxygen(graph, c2_idx) is None):
                continue
            # O1, O2 に H がないこと (peroxyacid は O2 に H あり → 除外)
            if (any(get_atom(graph, nb).symbol == "H" for nb in graph.adjacency[o1_idx])
                    or any(get_atom(graph, nb).symbol == "H" for nb in graph.adjacency[o2_idx])):
                continue

            # アシル鎖を収集して名前を生成
            chain1 = _collect_acid_chain(graph, c1_idx, {o1_idx}, get_atom)
            chain2 = _collect_acid_chain(graph, c2_idx, {o2_idx}, get_atom)

            def _acyl_name(chain: list[int]) -> str:
                n = len(chain)
                stem = CHAIN_PREFIX.get(n, f"C{n}")
                ene_locs, yne_locs = _chain_mb(graph, chain)
                if ene_locs or yne_locs:
                    return f"{stem}{_fmt_mb(ene_locs, yne_locs)}oyl"
                if n == 1:
                    return "formyl"
                return f"{stem}anoyl"

            acyl1 = _acyl_name(chain1)
            acyl2 = _acyl_name(chain2)

            if acyl1 == acyl2:
                return f"di{acyl1} peroxide"
            acyls = sorted([acyl1, acyl2])
            return f"{acyls[0]} {acyls[1]} peroxide"

    return None


def _name_cyanamide(graph, get_atom) -> str | None:
    """N-置換シアナミドの命名 (Phase 208, IUPAC 2013 P-66.4.1.1):

    R-NH-C≡N → N-alkylcyanamide
    R₂N-C≡N → N,N-dialkylcyanamide

    非置換 H₂N-C≡N は _RETAINED_NAMES テーブルで 'cyanamide' として処理する。
    """
    from .functional_group import get_bond_order
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter

    for idx in range(len(graph.atoms)):
        atom = get_atom(graph, idx)
        if atom.symbol != "N" or atom.in_ring:
            continue

        neighbors = graph.adjacency[idx]
        c_nbrs = [nb for nb in neighbors if get_atom(graph, nb).symbol == "C"]
        h_nbrs = [nb for nb in neighbors if get_atom(graph, nb).symbol == "H"]

        # Need >=1 C neighbor (for substitution)
        if len(c_nbrs) < 1:
            continue

        # Find the nitrile C (C≡N): must have exactly one triple-bonded N neighbor
        cn_c: int | None = None
        for c in c_nbrs:
            for c_nb in graph.adjacency[c]:
                c_nb_atom = get_atom(graph, c_nb)
                if c_nb_atom.symbol == "N" and get_bond_order(graph, c, c_nb) == 3.0:
                    cn_c = c
                    break
            if cn_c is not None:
                break

        if cn_c is None:
            continue

        # N-CN bond must be single
        if get_bond_order(graph, idx, cn_c) != 1.0:
            continue

        # Alkyl C neighbors = all C neighbors except the nitrile C
        alkyl_cs = [c for c in c_nbrs if c != cn_c]

        # Must have 1 or 2 alkyl substituents (otherwise it's unsubstituted cyanamide)
        if len(alkyl_cs) == 0:
            continue

        # Collect alkyl names
        alkyl_names = []
        for c in alkyl_cs:
            name = _name_carbon_substituent(graph, c, {idx})
            alkyl_names.append(name)

        counts = Counter(alkyl_names)
        prefix_parts = []
        for nm in sorted(set(alkyl_names)):
            cnt = counts[nm]
            mult = MULTIPLIER.get(cnt, f"{cnt}") if cnt > 1 else ""
            prefix_parts.append(f"{mult}{nm}")
        n_prefix = ",".join(["N"] * len(alkyl_cs))
        return f"{n_prefix}-{''.join(prefix_parts)}cyanamide"

    return None
