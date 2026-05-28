"""
smiles2iupac: SMILES から IUPAC 名を生成するライブラリ。

IUPACの命名ロジックはスクラッチ実装。
SMILESのパース・分子グラフ構築は RDKit に委譲。

対応範囲（Phase 1-3）:
  Phase 1-2:
  - 直鎖・分岐アルカン、アルコール、ハロゲン化物
  - カルボン酸、アルデヒド、ケトン
  - アルケン・アルキン
  - 立体化学 (R/S, E/Z)
  Phase 3 (complete):
  - シクロアルカン (cyclopropane ~ cyclodecane)
  - 置換シクロアルカン (methylcyclohexane 等)
  - シクロアルカンの官能基 (cyclohexanol, cyclohexanone 等)
  - ベンゼン
  - 置換ベンゼン (chlorobenzene, 1,2-dimethylbenzene 等)
"""

from __future__ import annotations

__all__ = ["smiles_to_iupac"]


# Phase 135/136: 保留名テーブル (IUPAC 2013 P-12.1, P-65.1.1.4, P-65.1.3)
# キー = RDKit canonical SMILES (stereo あり/なし両方)
_RETAINED_NAMES: dict[str, str] = {
    "NCC(=O)O":                          "glycine",
    # alanine
    "C[C@H](N)C(=O)O":                  "L-alanine",
    "C[C@@H](N)C(=O)O":                 "D-alanine",
    "CC(N)C(=O)O":                       "alanine",
    # valine
    "CC(C)[C@H](N)C(=O)O":              "L-valine",
    "CC(C)[C@@H](N)C(=O)O":             "D-valine",
    "CC(C)C(N)C(=O)O":                   "valine",
    # leucine
    "CC(C)C[C@H](N)C(=O)O":             "L-leucine",
    "CC(C)C[C@@H](N)C(=O)O":            "D-leucine",
    "CC(C)CC(N)C(=O)O":                  "leucine",
    # isoleucine
    "CC[C@H](C)[C@H](N)C(=O)O":         "L-isoleucine",
    # proline
    "O=C(O)[C@@H]1CCCN1":               "L-proline",
    "O=C(O)[C@H]1CCCN1":                "D-proline",
    "OC(=O)C1CCCN1":                     "proline",
    # phenylalanine
    "N[C@@H](Cc1ccccc1)C(=O)O":          "L-phenylalanine",
    "N[C@H](Cc1ccccc1)C(=O)O":          "D-phenylalanine",
    "NC(Cc1ccccc1)C(=O)O":               "phenylalanine",
    # tryptophan
    "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O": "L-tryptophan",
    "N[C@H](Cc1c[nH]c2ccccc12)C(=O)O":  "D-tryptophan",
    # methionine
    "CSCC[C@H](N)C(=O)O":               "L-methionine",
    "CSCC[C@@H](N)C(=O)O":              "D-methionine",
    "CSCCC(N)C(=O)O":                    "methionine",
    # serine
    "N[C@@H](CO)C(=O)O":                 "L-serine",
    "N[C@H](CO)C(=O)O":                 "D-serine",
    "NC(CO)C(=O)O":                      "serine",
    # threonine
    "C[C@H](O)[C@H](N)C(=O)O":          "L-threonine",
    "C[C@@H](O)[C@@H](N)C(=O)O":        "D-threonine",
    # cysteine
    "N[C@@H](CS)C(=O)O":                 "L-cysteine",
    "N[C@H](CS)C(=O)O":                 "D-cysteine",
    "NC(CS)C(=O)O":                      "cysteine",
    # tyrosine
    "N[C@@H](Cc1ccc(O)cc1)C(=O)O":       "L-tyrosine",
    "N[C@H](Cc1ccc(O)cc1)C(=O)O":       "D-tyrosine",
    # asparagine
    "NC(=O)C[C@H](N)C(=O)O":            "L-asparagine",
    "NC(=O)C[C@@H](N)C(=O)O":           "D-asparagine",
    "NC(=O)CC(N)C(=O)O":                 "asparagine",
    # glutamine
    "NC(=O)CC[C@H](N)C(=O)O":           "L-glutamine",
    "NC(=O)CC[C@@H](N)C(=O)O":          "D-glutamine",
    "NC(=O)CCC(N)C(=O)O":                "glutamine",
    # lysine
    "NCCCC[C@H](N)C(=O)O":              "L-lysine",
    "NCCCC[C@@H](N)C(=O)O":             "D-lysine",
    "NCCCCC(N)C(=O)O":                   "lysine",
    # arginine
    "N=C(N)NCCC[C@H](N)C(=O)O":         "L-arginine",
    "N=C(N)NCCC[C@@H](N)C(=O)O":        "D-arginine",
    # histidine
    "N[C@@H](Cc1cnc[nH]1)C(=O)O":        "L-histidine",
    "N[C@H](Cc1cnc[nH]1)C(=O)O":        "D-histidine",
    "NC(Cc1cnc[nH]1)C(=O)O":             "histidine",
    # aspartic acid
    "N[C@@H](CC(=O)O)C(=O)O":            "L-aspartic acid",
    "N[C@H](CC(=O)O)C(=O)O":            "D-aspartic acid",
    "NC(CC(=O)O)C(=O)O":                 "aspartic acid",
    # glutamic acid
    "N[C@@H](CCC(=O)O)C(=O)O":           "L-glutamic acid",
    "N[C@H](CCC(=O)O)C(=O)O":           "D-glutamic acid",
    "NC(CCC(=O)O)C(=O)O":                "glutamic acid",

    # ── Phase 136: 二塩基酸・ヒドロキシ酸 保留名 (P-65.1.1.4, P-65.1.3) ──
    # aliphatic diacids
    "O=C(O)C(=O)O":                      "oxalic acid",
    "O=C(O)CC(=O)O":                     "malonic acid",
    "O=C(O)CCC(=O)O":                    "succinic acid",
    "O=C(O)CCCC(=O)O":                   "glutaric acid",
    "O=C(O)CCCCC(=O)O":                  "adipic acid",
    "O=C(O)CCCCCC(=O)O":                 "pimelic acid",
    "O=C(O)CCCCCCC(=O)O":               "suberic acid",
    "O=C(O)CCCCCCCC(=O)O":              "azelaic acid",
    "O=C(O)CCCCCCCCC(=O)O":             "sebacic acid",
    # unsaturated diacids
    "O=C(O)/C=C\\C(=O)O":               "maleic acid",
    "O=C(O)/C=C/C(=O)O":                "fumaric acid",
    # cinnamic acid
    "O=C(O)/C=C/c1ccccc1":              "(E)-cinnamic acid",
    "O=C(O)/C=C\\c1ccccc1":             "(Z)-cinnamic acid",
    # aromatic diacids (phthalic, isophthalic, terephthalic) P-65.1.1.4
    "O=C(O)c1ccccc1C(=O)O":             "phthalic acid",
    "O=C(O)c1cccc(C(=O)O)c1":           "isophthalic acid",
    "O=C(O)c1ccc(C(=O)O)cc1":           "terephthalic acid",
    # lactic acid P-65.1.1.4
    "C[C@@H](O)C(=O)O":                 "L-lactic acid",
    "C[C@H](O)C(=O)O":                  "D-lactic acid",
    "CC(O)C(=O)O":                       "lactic acid",
    # malic acid
    "O=C(O)C[C@H](O)C(=O)O":           "L-malic acid",
    "O=C(O)C[C@@H](O)C(=O)O":          "D-malic acid",
    "O=C(O)CC(O)C(=O)O":               "malic acid",
    # tartaric acid
    "O=C(O)[C@@H](O)[C@H](O)C(=O)O":   "L-tartaric acid",
    "O=C(O)[C@H](O)[C@@H](O)C(=O)O":   "D-tartaric acid",
    "O=C(O)[C@@H](O)[C@@H](O)C(=O)O":  "meso-tartaric acid",
    "O=C(O)C(O)C(O)C(=O)O":            "tartaric acid",
    # citric acid
    "O=C(O)CC(O)(CC(=O)O)C(=O)O":      "citric acid",

    # ── Phase 137: 核酸塩基 保留名 (IUPAC 2013 P-14.5) ──
    "Nc1ncnc2[nH]cnc12":               "adenine",
    "Nc1nc2[nH]cnc2c(=O)[nH]1":        "guanine",
    "Nc1ccnc(=O)[nH]1":                "cytosine",
    "O=c1cc[nH]c(=O)[nH]1":           "uracil",
    "Cc1c[nH]c(=O)[nH]c1=O":          "thymine",
    "O=c1[nH]cnc2[nH]cnc12":          "hypoxanthine",
    "O=c1[nH]c(=O)c2[nH]cnc2[nH]1":  "xanthine",
    "c1ncc2[nH]cnc2n1":               "9H-purine",

    # ── Phase 146: 無機オキソ酸・オキソ酸アニオン 保留名 ──
    "O=C([O-])O":                      "bicarbonate",
    "O=C([O-])[O-]":                   "carbonate",
    "O=S(=O)(O)O":                     "sulfuric acid",
    "O=S(O)O":                         "sulfurous acid",
    "O=[N+]([O-])O":                   "nitric acid",
    "O=NO":                            "nitrous acid",
    "Cl":                              "hydrochloric acid",
    "OB(O)O":                          "boric acid",
    "O=C(O)O":                         "carbonic acid",
    "O=[N+]([O-])[O-]":               "nitrate",
    "O=S(=O)([O-])[O-]":              "sulfate",

    # ── Phase 144: 無機イオン・ヒドロキシルアミン 保留名 ──
    "NO":                              "hydroxylamine",
    "[NH4+]":                          "azanium",
    "[OH-]":                           "hydroxide",
    "[F-]":                            "fluoride",
    "[Cl-]":                           "chloride",
    "[Br-]":                           "bromide",
    "[I-]":                            "iodide",
    "[Na+]":                           "sodium",
    "[K+]":                            "potassium",
    "[Li+]":                           "lithium",
    "[Ca+2]":                          "calcium",
    "[Mg+2]":                          "magnesium",
    "[Zn+2]":                          "zinc",
    "[Fe+2]":                          "iron(2+)",
    "[Fe+3]":                          "iron(3+)",
    "[Cu+]":                           "copper(1+)",
    "[Cu+2]":                          "copper(2+)",
    "[Al+3]":                          "aluminium",

    # ── Phase 143: 無機リン・ホウ素・ケイ素 保留名 ──
    "P":                               "phosphane",
    "PP":                              "diphosphane",
    "O=P(O)(O)O":                      "phosphoric acid",
    "OP(O)O":                          "phosphorous acid",
    "O=[PH](O)O":                      "hypophosphorous acid",
    "B":                               "borane",
    "BB":                              "diborane",
    "[SiH4]":                          "silane",
    "[SiH3][SiH3]":                    "disilane",
}


def _try_retained_name(graph) -> str | None:
    """canonical SMILES で保留名テーブルを照合する。"""
    if graph.rdkit_mol is None:
        return None
    from rdkit.Chem import MolToSmiles
    canon = MolToSmiles(graph.rdkit_mol)
    return _RETAINED_NAMES.get(canon)


def _name_multicomponent(smiles: str) -> str:
    """複数成分 SMILES (ドット区切り) の命名 — 塩形式に対応  (Phase 147)"""
    from rdkit import Chem
    from collections import Counter
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    frags = Chem.GetMolFrags(mol, asMols=True)
    # 各成分を名前付きリストに
    cat_names: list[str] = []
    ani_names: list[str] = []
    neu_names: list[str] = []
    for frag in frags:
        frag_smi = Chem.MolToSmiles(frag)
        name = smiles_to_iupac(frag_smi)
        charge = sum(a.GetFormalCharge() for a in frag.GetAtoms())
        if charge > 0:
            cat_names.append(name)
        elif charge < 0:
            ani_names.append(name)
        else:
            neu_names.append(name)

    def _compact(names: list[str]) -> str:
        counts = Counter(names)
        from .constants import MULTIPLIER
        parts = []
        for nm in sorted(set(names), key=lambda n: names.index(n)):
            cnt = counts[nm]
            if cnt > 1:
                mult = MULTIPLIER.get(cnt, f"{cnt}")
                parts.append(f"{mult}{nm}")
            else:
                parts.append(nm)
        return " ".join(parts)

    parts = []
    if cat_names:
        parts.append(_compact(cat_names))
    if ani_names:
        parts.append(_compact(ani_names))
    if neu_names:
        parts.append(_compact(neu_names))
    return " ".join(parts)


def _is_isocyanic_acid(graph) -> bool:
    """H-N=C=O パターン (isocyanic acid) を検出する。"""
    from .molecule_analyzer import get_atom, get_bond_order
    for c_idx in range(len(graph.atoms)):
        if get_atom(graph, c_idx).symbol != "C":
            continue
        has_dbl_n_h = False
        has_dbl_o = False
        for nb_idx in graph.adjacency[c_idx]:
            nb = get_atom(graph, nb_idx)
            bo = get_bond_order(graph, c_idx, nb_idx)
            if nb.symbol == "N" and bo == 2.0:
                n_non_ch = [n for n in graph.adjacency[nb_idx]
                            if n != c_idx and get_atom(graph, n).symbol not in ("H",)]
                n_h = sum(1 for n in graph.adjacency[nb_idx]
                          if get_atom(graph, n).symbol == "H")
                if not n_non_ch and n_h >= 1:
                    has_dbl_n_h = True
            elif nb.symbol == "O" and bo == 2.0:
                o_non_ch = [n for n in graph.adjacency[nb_idx]
                            if n != c_idx and get_atom(graph, n).symbol not in ("H",)]
                if not o_non_ch:
                    has_dbl_o = True
        if has_dbl_n_h and has_dbl_o:
            heavy_other = [nb for nb in graph.adjacency[c_idx]
                           if get_atom(graph, nb).symbol not in ("H", "N", "O")]
            if not heavy_other:
                return True
    return False


def smiles_to_iupac(smiles: str) -> str:
    """
    SMILES 文字列から IUPAC 系統名を生成する。

    Args:
        smiles: 有効な SMILES 文字列

    Returns:
        IUPAC 系統名文字列

    Raises:
        ValueError: 無効な SMILES
        NotImplementedError: 未対応構造（縮合多環系など）
    """
    # 遅延インポート（rdkit が不要なモジュールを単独テスト可能にする）
    from .molecule_analyzer import build_molecule_graph, get_atom
    from .functional_group import detect_groups, principal_group, FunctionalGroup
    from .chain_finder import find_principal_chain, get_multiple_bond_locants
    from .substituent import collect_substituents
    from .stereochemistry import assign_stereochemistry
    from .name_assembler import assemble_name
    from .ring_handler import (
        has_ring, find_rings, find_principal_ring,
        collect_ring_substituents, assemble_ring_name,
        has_hetero_ring,
    )
    from .heterocycle_handler import name_heterocycle, _try_fused_hetero_retained
    from .polycyclic_handler import name_polycyclic

    # ─── Phase 147: 複数成分 SMILES (塩形式) ─────────────────────────────
    if "." in smiles:
        return _name_multicomponent(smiles)

    # ─── 1. 分子グラフ構築 ───────────────────────────────────────────
    graph = build_molecule_graph(smiles)

    # Phase 135/136: 保留名チェック (早期チェック)
    _aa = _try_retained_name(graph)
    if _aa is not None:
        return _aa

    # ─── 2. 官能基を先行検出（ester/acid_halide は専用パス）───────────
    # Phase 128: イソシアン酸 H-N=C=O (isocyanic acid) の早期検出
    if _is_isocyanic_acid(graph):
        return "isocyanic acid"

    from .functional_group import detect_groups as _detect, principal_group as _pg
    _groups = _detect(graph)
    _pgrp = _pg(_groups)

    if _pgrp is not None and _pgrp.group_type == "diester":
        # Phase 127: 環状無水物 (oxa-dione) は環状パスへ (両 C が in_ring かつ ring O あり)
        _d_cs = [ai for ai in _pgrp.atom_indices
                 if get_atom(graph, ai).symbol == "C" and get_atom(graph, ai).in_ring]
        _d_ring_os = [ai for ai in _pgrp.atom_indices
                      if get_atom(graph, ai).symbol == "O" and get_atom(graph, ai).in_ring]
        if not (_d_cs and _d_ring_os):
            return _name_diester(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "ester":
        # Phase 30: ラクトン (cyclic ester) は環状パスへ
        _ester_c = _pgrp.atom_indices[0]
        _ring_o_in_ester = [ai for ai in _pgrp.atom_indices[1:]
                            if get_atom(graph, ai).in_ring]
        if not (get_atom(graph, _ester_c).in_ring and _ring_o_in_ester):
            return _name_ester(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "diacid_halide":
        return _name_diacid_halide(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "acid_halide":
        return _name_acid_halide(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "anhydride":
        return _name_anhydride(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "carbonate":
        return _name_carbonate(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "sulfonyl_chloride":
        return _name_sulfonyl_chloride(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "chloroformate":
        return _name_chloroformate(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "sulfonic_acid":
        return _name_sulfonic_acid(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "sulfinic_acid":
        return _name_sulfinic_acid(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "sulfonate_ester":
        return _name_sulfonate_sulfinate_ester(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "sulfinate_ester":
        return _name_sulfonate_sulfinate_ester(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type in ("sulfoxide", "sulfone"):
        return _name_sulfoxide_sulfone(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "sulfonamide":
        return _name_sulfonamide(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "sulfide":
        return _name_sulfide(graph, _pgrp, get_atom)

    # Phase 143/145: リン・ホウ素・ケイ素化合物
    if _pgrp is not None and _pgrp.group_type == "phosphate_ester":
        return _name_phosphate_ester(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "phosphonate_ester":
        return _name_phosphonate_ester(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "phosphinate_ester":
        return _name_phosphinate_ester(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "phosphonic_acid":
        return _name_phosphonic_acid(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "phosphinic_acid":
        return _name_phosphinic_acid(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "phosphane":
        return _name_phosphane(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "boronic_acid":
        return _name_boronic_acid(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "borinic_acid":
        return _name_borinic_acid(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "borane_org":
        return _name_organic_borane(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "silane_org":
        return _name_organic_silane(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "ammonium":
        return _name_ammonium(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "carbamate":
        # Phase 154: 環状カルバメート (oxazolidinone 等) は環状パスへ
        _carb_c = _pgrp.atom_indices[0]
        if not get_atom(graph, _carb_c).in_ring:
            return _name_carbamate(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "hydroperoxide":
        return _name_hydroperoxide(graph, _pgrp, get_atom)

    # Phase 87: チオ尿素 NC(=S)N (thioamide より先に確認)
    _thiourea = _name_thiourea_if_match(graph, get_atom)
    if _thiourea is not None:
        return _thiourea

    if _pgrp is not None and _pgrp.group_type == "thioamide":
        return _name_thioamide(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "peroxide":
        return _name_peroxide(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "thioester":
        return _name_thioester(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "disulfide":
        return _name_disulfide(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "nitroso":
        return _name_nitroso(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "azide":
        return _name_azide(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "isocyanate":
        return _name_isocyanate(graph, _pgrp, get_atom)

    # Phase 68–73: 新規官能基
    if _pgrp is not None and _pgrp.group_type == "isothiocyanate":
        return _name_isothiocyanate(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "cyanate":
        return _name_cyanate(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "thiocyanate":
        return _name_thiocyanate(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "carboxylate":
        return _name_carboxylate(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "dicarboxylate":
        return _name_dicarboxylate(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type in ("thioic_s_acid", "thioic_o_acid", "dithioic_acid"):
        return _name_thioic_acid(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "nitrate_ester":
        return _name_nitrate_ester(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "nitrite_ester":
        return _name_nitrite_ester(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "carbamic_acid":
        return _name_carbamic_acid(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "carbodiimide":
        return _name_carbodiimide(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "hydrazide":
        return _name_hydrazide(graph, _pgrp, get_atom)

    if _pgrp is not None and _pgrp.group_type == "peroxyacid":
        return _name_peroxyacid(graph, _pgrp, get_atom)

    # Phase 80: セミカルバゾン / チオセミカルバゾン
    if _pgrp is not None and _pgrp.group_type in (
        "semicarbazone", "aldsemicarbazone", "thiosemicarbazone", "aldthiosemicarbazone"
    ):
        return _name_semicarbazone(graph, _pgrp, get_atom)

    # Phase 79: N-置換ヒドラゾン
    if _pgrp is not None and _pgrp.group_type in ("aldhydrazone", "kethydrazone"):
        result = _name_substituted_hydrazone(graph, _pgrp, get_atom)
        if result is not None:
            return result

    # Phase 67: スルファミド / スルファミン酸 (炭素なし特別ケース)
    if _pgrp is not None and _pgrp.group_type == "sulfamide":
        return "sulfamide"
    if _pgrp is not None and _pgrp.group_type == "sulfamic_acid":
        return "sulfamic acid"

    # Phase 49: urea 特別ケース (NC(=O)N)
    if _is_urea(graph, get_atom):
        return "urea"

    # Phase 83: N-置換尿素 (CNC(=O)N → N-methylurea 等)
    _sub_urea = _name_substituted_urea_if_match(graph, get_atom)
    if _sub_urea is not None:
        return _sub_urea

    # Phase 125: アミン N-オキシド (R₃N→O)
    _n_oxide_name = _name_amine_n_oxide(graph, get_atom)
    if _n_oxide_name is not None:
        return _n_oxide_name

    # Phase 123: ジアゾ化合物 (C=N=N, diazo prefix)
    _diazo_name = _name_diazo_compound(graph, get_atom)
    if _diazo_name is not None:
        return _diazo_name

    # Phase 113: ヒドラジン (N-N 単結合、環外、カルボニル隣接なし)
    _hydrazine_name = _name_hydrazine_compound(graph, get_atom)
    if _hydrazine_name is not None:
        return _hydrazine_name

    # Phase 115: アゾ化合物 (N=N 二重結合)
    _azo_name = _name_azo_compound(graph, get_atom)
    if _azo_name is not None:
        return _azo_name

    # ─── 3. 環状 or 非環状の分岐 ─────────────────────────────────────
    if has_ring(graph):
        # Phase 28/32: アンカー C が環外にある場合 → acyclic path
        # Phase 28: anchor_c1=False (alcohol/amine/ketone 等) でアンカー C が環外
        # Phase 32: anchor_c1=True (carboxylic_acid 等) でアンカー C が環に隣接しない
        if _pgrp is not None and _pgrp.atom_indices:
            from .constants import FUNCTIONAL_GROUPS as _FGS
            _spec = _FGS.get(_pgrp.group_type)
            _anchor_idx = _pgrp.atom_indices[0]
            _anchor_in_ring = get_atom(graph, _anchor_idx).in_ring
            if _spec is not None and not _anchor_in_ring:
                _anchor_adj_ring = any(
                    get_atom(graph, nb).in_ring
                    for nb in graph.adjacency[_anchor_idx]
                )
                # Route to acyclic when:
                # - anchor_c1=False and anchor not in ring (Phase 28)
                # - anchor_c1=True and anchor not adjacent to ring (Phase 32)
                # Phase 40: alkene/alkyne で環の方が大きい場合は ring path を優先
                if not _spec.anchor_c1 or not _anchor_adj_ring:
                    if _pgrp.group_type in ("alkene", "alkyne"):
                        from .molecule_analyzer import non_ring_carbon_indices
                        _nr_c = len(non_ring_carbon_indices(graph))
                        _rings_tmp = find_rings(graph)
                        _max_ring = max((len(r) for r in _rings_tmp), default=0)
                        # ヘテロ芳香族環も考慮 (find_rings は純炭素環のみ)
                        if _max_ring == 0 and graph.ring_atom_sets:
                            _max_ring = max((len(r) for r in graph.ring_atom_sets), default=0)
                        if _max_ring > _nr_c:
                            pass  # ring path へ
                        else:
                            return _name_acyclic(graph, detect_groups, principal_group,
                                                 find_principal_chain, get_multiple_bond_locants,
                                                 collect_substituents, assign_stereochemistry,
                                                 assemble_name, get_atom)
                    else:
                        return _name_acyclic(graph, detect_groups, principal_group,
                                             find_principal_chain, get_multiple_bond_locants,
                                             collect_substituents, assign_stereochemistry,
                                             assemble_name, get_atom)

        # 3a. スピロ・架橋二環（Phase 15）
        poly_name = name_polycyclic(graph)
        if poly_name is not None:
            return poly_name

        # 3b. 縮合環保留名テーブル (Phase 133: 全環系対象、ヘテロ・炭素両方)
        fused_name = _try_fused_hetero_retained(graph)
        if fused_name is not None:
            return fused_name

        # 3c. ヘテロ環保留名チェック（Phase 12）
        if has_hetero_ring(graph):
            hetero_name = name_heterocycle(graph)
            if hetero_name is not None:
                return hetero_name

        return _name_cyclic(graph, find_rings, find_principal_ring,
                            collect_ring_substituents, assemble_ring_name,
                            detect_groups, principal_group, assign_stereochemistry,
                            get_atom)
    else:
        return _name_acyclic(graph, detect_groups, principal_group,
                             find_principal_chain, get_multiple_bond_locants,
                             collect_substituents, assign_stereochemistry,
                             assemble_name, get_atom)


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

    acid_carbons = _collect_acid_chain(graph, c1, ester_os, get_atom)
    n_acid = len(acid_carbons)

    stem = CHAIN_PREFIX.get(n_acid, f"C{n_acid}")
    _ene_de, _yne_de = _chain_multiple_bonds(graph, acid_carbons)
    if _ene_de or _yne_de:
        from .name_assembler import _format_multiple_bonds as _fmt_de
        acid_name = f"{stem}{_fmt_de(_ene_de, _yne_de)}edioate"
    else:
        acid_name = f"{stem}anedioate"

    # alkyl 部分の組み立て
    if alkyl1 == alkyl2:
        mult = MULTIPLIER.get(2, "di")
        return f"{mult}{alkyl1} {acid_name}"
    return f"{min(alkyl1, alkyl2)} {max(alkyl1, alkyl2)} {acid_name}"


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
    # c1 の O 隣接を除外して主鎖収集
    o_excl = {nb for nb in graph.adjacency[c1] if get_atom(graph, nb).symbol == "O"}
    chain = _collect_acid_chain(graph, c1, o_excl, get_atom)
    n = len(chain)
    # 保留ジアニオン名
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
    # O-側のインデックスを除外して主鎖収集
    o_idxs = {nb for nb in graph.adjacency[carbonyl_c] if get_atom(graph, nb).symbol == "O"}
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
        return f"{stem}{_fmt(ene, yne)}oate"
    return f"{stem}anoate"


def _name_thioic_acid(graph, pgrp, get_atom) -> str:
    """チオカルボン酸命名 (Phase 149)"""
    from .constants import CHAIN_PREFIX
    gtype = pgrp.group_type
    carbonyl_c = pgrp.atom_indices[0]
    chalcogen_idxs = {nb for nb in graph.adjacency[carbonyl_c]
                      if get_atom(graph, nb).symbol in ("O", "S")}
    acid_chain = _collect_acid_chain(graph, carbonyl_c, chalcogen_idxs, get_atom)
    n = len(acid_chain)
    stem = CHAIN_PREFIX.get(n, f"C{n}")
    ene, yne = _chain_multiple_bonds(graph, acid_chain)
    if ene or yne:
        from .name_assembler import _format_multiple_bonds as _fmt
        mb = _fmt(ene, yne)
    else:
        mb = ""
    if gtype == "thioic_s_acid":
        if not mb:
            return f"{stem}anethioic S-acid"
        return f"{stem}{mb}thioic S-acid"
    elif gtype == "thioic_o_acid":
        if not mb:
            return f"{stem}anethioic O-acid"
        return f"{stem}{mb}thioic O-acid"
    else:  # dithioic_acid
        if not mb:
            return f"{stem}anedithioic acid"
        return f"{stem}{mb}dithioic acid"


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

    return f"{alkyl_name} {acid_name}"


def _collect_acid_chain(graph, start_c, excluded_set, get_atom) -> list[int]:
    """カルボニル C から excluded_set を除いた炭素鎖を DFS で列挙。"""
    visited: set[int] = set(excluded_set)
    result: list[int] = []

    def dfs(idx: int) -> None:
        if idx in visited:
            return
        atom = get_atom(graph, idx)
        if atom.symbol != "C":
            return
        visited.add(idx)
        result.append(idx)
        for nb in graph.adjacency[idx]:
            dfs(nb)

    dfs(start_c)
    return result


def _chain_multiple_bonds(graph, chain: list[int]) -> tuple[list[int], list[int]]:
    """chain 内の連続原子ペアの多重結合ロカントを返す。(ene_locs, yne_locs)"""
    from .molecule_analyzer import get_bond_order
    ene: list[int] = []
    yne: list[int] = []
    for i in range(len(chain) - 1):
        bo = get_bond_order(graph, chain[i], chain[i + 1])
        if bo == 2.0:
            ene.append(i + 1)
        elif bo == 3.0:
            yne.append(i + 1)
    return ene, yne


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

    # 芳香環に直接結合したカルボニル → benzoyl 型 (Phase 107)
    for nb_idx in graph.adjacency[carbonyl_c]:
        if nb_idx in halogen_set:
            continue
        nb = get_atom(graph, nb_idx)
        if nb.symbol == "C" and nb.is_aromatic and nb.in_ring:
            # 純粋ベンゼン環かチェック
            ring_atoms = next(
                (set(rt) for rt in (graph.ring_atom_sets or []) if nb_idx in rt), set()
            )
            if (len(ring_atoms) == 6
                    and all(get_atom(graph, a).symbol == "C" for a in ring_atoms)):
                return f"benzoyl {halide_name}"
            # 他の縮合芳香環等は後続の通常パスへ
            break

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

    # 鎖上の置換基を収集 (Phase 72)
    locant_map = {c: i + 1 for i, c in enumerate(acid_carbons)}
    # 主官能基原子: carbonyl_c の =O, ハロゲン, carbonyl_c 自身
    principal_atoms = set(pgrp.atom_indices) | halogen_set
    for nb_idx in graph.adjacency[carbonyl_c]:
        if get_atom(graph, nb_idx).symbol == "O":
            principal_atoms.add(nb_idx)
    subs = collect_substituents(graph, acid_carbons, locant_map, list(principal_atoms))
    if not subs:
        # Phase 120: 置換基なし 1C/2C の保留名 (IUPAC 2013 PIN)
        if n_acid == 1 and halide_name == "chloride":
            return "formyl chloride"
        if n_acid == 2 and halide_name == "chloride":
            return "acetyl chloride"
        return base
    prefix = _build_prefix(subs)
    return f"{prefix}{base}"


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

    c1 = carbonyl_cs[0]
    halogen_nb_set: set[int] = set()
    for nb_idx in graph.adjacency[c1]:
        nb = get_atom(graph, nb_idx)
        if nb.symbol in halogen_syms:
            halogen_nb_set.add(nb_idx)

    acid_carbons = _collect_acid_chain(graph, c1, halogen_nb_set, get_atom)
    n = len(acid_carbons)
    stem = CHAIN_PREFIX.get(n, f"C{n}")

    halide_names_sorted = sorted(halide_names)
    if halide_names_sorted[0] == halide_names_sorted[1]:
        halide_str = f"di{halide_names_sorted[0]}"
    else:
        halide_str = f"{halide_names_sorted[0]} {halide_names_sorted[1]}"

    return f"{stem}anedioyl {halide_str}"


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
    if c_acid:
        aryl = _aryl_sulfonyl_prefix(graph, c_acid[0], s_idx, get_atom)
        if aryl is not None:
            acid_name = f"{aryl}{kind}"
        else:
            acid_chain = _collect_acid_chain(graph, c_acid[0], {s_idx}, get_atom)
            n = len(acid_chain)
            acid_stem = CHAIN_PREFIX.get(n, f"C{n}")
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

    return f"{ester_alkyl} {acid_name}"


def _name_sulfoxide_sulfone(graph, pgrp, get_atom) -> str:
    """
    スルホキシド・スルホン命名: {alkyl1} {alkyl2} sulfoxide/sulfone
    例: CS(=O)C → dimethyl sulfoxide
        CCS(=O)(=O)C → ethyl methyl sulfone
    """
    from .constants import CHAIN_PREFIX, MULTIPLIER
    from .substituent import _name_carbon_substituent

    s_idx = pgrp.atom_indices[0]
    kind = "sulfoxide" if pgrp.group_type == "sulfoxide" else "sulfone"

    c_neighbors = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if len(c_neighbors) < 2:
        return f"sulfoxide"

    names = sorted(_name_carbon_substituent(graph, c, {s_idx}) for c in c_neighbors)

    if names[0] == names[1]:
        mult = MULTIPLIER.get(2, "di")
        return f"{mult}{names[0]} {kind}"
    return f"{names[0]} {names[1]} {kind}"


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
    if aryl is not None:
        # N 上置換基は後続の N-subs ロジックで追記するため base だけ作る
        base = f"{aryl}sulfonamide"
    else:
        acid_chain = _collect_acid_chain(graph, c_on_s[0], {s_idx}, get_atom)
        stem = CHAIN_PREFIX.get(len(acid_chain), f"C{len(acid_chain)}")
        _ene_sn, _yne_sn = _chain_multiple_bonds(graph, acid_chain)
        if _ene_sn or _yne_sn:
            from .name_assembler import _format_multiple_bonds as _fmt_sn
            base = f"{stem}{_fmt_sn(_ene_sn, _yne_sn)}esulfonamide"
        else:
            base = f"{stem}anesulfonamide"

    # N 上の置換基
    n_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "N"]
    if not n_on_s:
        return base

    n_idx = n_on_s[0]
    c_on_n = [nb for nb in graph.adjacency[n_idx]
               if get_atom(graph, nb).symbol == "C"]
    if not c_on_n:
        return base

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in c_on_n]
    sub_counts = Counter(n_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        if cnt == 1:
            prefix_parts.append(f"N-{sub}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub}")

    prefix = "-".join(prefix_parts)
    return f"{prefix}{base}"


def _aryl_sulfonyl_prefix(graph, c_start, s_idx, get_atom) -> str | None:
    """
    c_start が 6員ベンゼン環の芳香族 C であれば環名 ("benzene", "4-methylbenzene" 等) を返す。
    それ以外は None。
    """
    atom = get_atom(graph, c_start)
    if not (atom.is_aromatic and atom.in_ring):
        return None
    ring_atoms = next(
        (rt for rt in (graph.ring_atom_sets or []) if c_start in rt), None
    )
    if ring_atoms is None or len(ring_atoms) != 6:
        return None
    if not all(get_atom(graph, a).symbol == "C" for a in ring_atoms):
        return None
    from .ring_handler import _assign_ring_locants, collect_ring_substituents, assemble_ring_name
    ring_chain = _assign_ring_locants(graph, list(ring_atoms), True, "alkane", [s_idx])
    substituents = collect_ring_substituents(graph, ring_chain, [s_idx])
    return assemble_ring_name(ring_chain, substituents, "alkane", None, [])


def _name_sulfonic_acid(graph, pgrp, get_atom) -> str:
    """
    スルホン酸命名: {stem}anesulfonic acid
    例: CS(=O)(=O)O → methanesulfonic acid
        c1ccc(S(=O)(=O)O)cc1 → benzenesulfonic acid
    """
    from .constants import CHAIN_PREFIX
    s_idx = pgrp.atom_indices[0]
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_s:
        return "sulfonic acid"
    aryl = _aryl_sulfonyl_prefix(graph, c_on_s[0], s_idx, get_atom)
    if aryl is not None:
        return f"{aryl}sulfonic acid"
    acid_chain = _collect_acid_chain(graph, c_on_s[0], {s_idx}, get_atom)
    stem = CHAIN_PREFIX.get(len(acid_chain), f"C{len(acid_chain)}")
    _ene_so, _yne_so = _chain_multiple_bonds(graph, acid_chain)
    if _ene_so or _yne_so:
        from .name_assembler import _format_multiple_bonds as _fmt_so
        return f"{stem}{_fmt_so(_ene_so, _yne_so)}esulfonic acid"
    return f"{stem}anesulfonic acid"


def _name_sulfinic_acid(graph, pgrp, get_atom) -> str:
    """
    スルフィン酸命名: {stem}anesulfinic acid (Phase 38)
    例: CS(=O)O → methanesulfinic acid
        c1ccc(S(=O)O)cc1 → benzenesulfinic acid
    """
    from .constants import CHAIN_PREFIX
    s_idx = pgrp.atom_indices[0]
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_s:
        return "sulfinic acid"
    aryl = _aryl_sulfonyl_prefix(graph, c_on_s[0], s_idx, get_atom)
    if aryl is not None:
        return f"{aryl}sulfinic acid"
    acid_chain = _collect_acid_chain(graph, c_on_s[0], {s_idx}, get_atom)
    stem = CHAIN_PREFIX.get(len(acid_chain), f"C{len(acid_chain)}")
    _ene_si, _yne_si = _chain_multiple_bonds(graph, acid_chain)
    if _ene_si or _yne_si:
        from .name_assembler import _format_multiple_bonds as _fmt_si
        return f"{stem}{_fmt_si(_ene_si, _yne_si)}esulfinic acid"
    return f"{stem}anesulfinic acid"


def _name_sulfonyl_chloride(graph, pgrp, get_atom) -> str:
    """スルホニルクロライド命名: {stem}anesulfonyl chloride (Phase 59/110)"""
    from .constants import CHAIN_PREFIX
    s_idx = pgrp.atom_indices[0]
    c_on_s = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_on_s:
        return "sulfonyl chloride"
    aryl = _aryl_sulfonyl_prefix(graph, c_on_s[0], s_idx, get_atom)
    if aryl is not None:
        return f"{aryl}sulfonyl chloride"
    acid_chain = _collect_acid_chain(graph, c_on_s[0], {s_idx}, get_atom)
    stem = CHAIN_PREFIX.get(len(acid_chain), f"C{len(acid_chain)}")
    _ene_sc, _yne_sc = _chain_multiple_bonds(graph, acid_chain)
    if _ene_sc or _yne_sc:
        from .name_assembler import _format_multiple_bonds as _fmt_sc
        return f"{stem}{_fmt_sc(_ene_sc, _yne_sc)}esulfonyl chloride"
    return f"{stem}anesulfonyl chloride"


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
    チオエーテル命名: {alkyl1} {alkyl2} sulfide
    例: CSC → dimethyl sulfide, CSCC → ethyl methyl sulfide
        CSc1ccccc1 → methyl phenyl sulfide
    """
    from .constants import CHAIN_PREFIX, MULTIPLIER
    from .substituent import _name_carbon_substituent

    def _alkyl_or_aryl(c_idx: int) -> str:
        atom = get_atom(graph, c_idx)
        if atom.is_aromatic and atom.in_ring:
            ring_neighbors = [nb for nb in graph.adjacency[c_idx]
                              if get_atom(graph, nb).in_ring]
            if len(ring_neighbors) >= 2:
                return "phenyl"
        return _name_carbon_substituent(graph, c_idx, {s_idx})

    s_idx = pgrp.atom_indices[0]
    c_neighbors = [nb for nb in graph.adjacency[s_idx] if get_atom(graph, nb).symbol == "C"]
    if len(c_neighbors) < 2:
        return "sulfide"
    names = sorted(_alkyl_or_aryl(c) for c in c_neighbors)
    if names[0] == names[1]:
        mult = MULTIPLIER.get(2, "di")
        return f"{mult}{names[0]} sulfide"
    return f"{names[0]} {names[1]} sulfide"


def _name_carbonate(graph, pgrp, get_atom) -> str:
    """
    炭酸エステル命名: {alkyl1} {alkyl2} carbonate (Phase 39)
    例: COC(=O)OC → dimethyl carbonate
        CCOC(=O)OCC → diethyl carbonate
        COC(=O)OCC → ethyl methyl carbonate (アルファベット順)
    """
    from .constants import MULTIPLIER
    from .substituent import _name_carbon_substituent
    from collections import Counter

    carbonyl_c = pgrp.atom_indices[0]
    single_o_list = [nb for nb in graph.adjacency[carbonyl_c]
                     if get_atom(graph, nb).symbol == "O"
                     and any(get_atom(graph, on).symbol == "C"
                             for on in graph.adjacency[nb] if on != carbonyl_c)]
    alkyl_names = []
    for o_idx in single_o_list:
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
    return " ".join(parts) + " carbonate"


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
    parts = []
    for sub in sorted(counts):
        n = counts[sub]
        if sub == "hydrogen":
            mult = MULTIPLIER.get(n, "") if n > 1 else ""
            parts.append(f"{mult}hydrogen")
        else:
            mult = MULTIPLIER.get(n, "") if n > 1 else ""
            parts.append(f"{mult}{sub}")
    return " ".join(parts) + " phosphate"


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


def _name_phosphonic_acid(graph, pgrp, get_atom) -> str:
    """ホスホン酸: {alkyl}phosphonic acid  (Phase 143)"""
    from .substituent import _name_carbon_substituent
    p_idx = pgrp.atom_indices[0]
    c_neighbors = [nb for nb in graph.adjacency[p_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_neighbors:
        return "phosphonic acid"
    sub = _name_carbon_substituent(graph, c_neighbors[0], {p_idx})
    return f"{sub}phosphonic acid"


def _name_phosphinic_acid(graph, pgrp, get_atom) -> str:
    """ホスフィン酸: di{alkyl}phosphinic acid  (Phase 143)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter
    p_idx = pgrp.atom_indices[0]
    c_neighbors = [nb for nb in graph.adjacency[p_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_neighbors:
        return "phosphinic acid"
    names = sorted(_name_carbon_substituent(graph, c, {p_idx}) for c in c_neighbors)
    counts = Counter(names)
    parts = []
    for sub in sorted(counts):
        n = counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    return "".join(parts) + "phosphinic acid"


def _name_phosphane(graph, pgrp, get_atom) -> str:
    """ホスファン: {alkyl(s)}phosphane  (Phase 143)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter
    p_idx = pgrp.atom_indices[0]
    c_neighbors = [nb for nb in graph.adjacency[p_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_neighbors:
        return "phosphane"
    names = sorted(_name_carbon_substituent(graph, c, {p_idx}) for c in c_neighbors)
    counts = Counter(names)
    parts = []
    for sub in sorted(counts):
        n = counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    return "".join(parts) + "phosphane"


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
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter
    b_idx = pgrp.atom_indices[0]
    c_neighbors = [nb for nb in graph.adjacency[b_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_neighbors:
        return "borinic acid"
    names = sorted(_name_carbon_substituent(graph, c, {b_idx}) for c in c_neighbors)
    counts = Counter(names)
    parts = []
    for sub in sorted(counts):
        n = counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    return "".join(parts) + "borinic acid"


def _name_organic_borane(graph, pgrp, get_atom) -> str:
    """有機ボラン: {alkyl(s)}borane  (Phase 143)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter
    b_idx = pgrp.atom_indices[0]
    c_neighbors = [nb for nb in graph.adjacency[b_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_neighbors:
        return "borane"
    names = sorted(_name_carbon_substituent(graph, c, {b_idx}) for c in c_neighbors)
    counts = Counter(names)
    parts = []
    for sub in sorted(counts):
        n = counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    return "".join(parts) + "borane"


def _name_ammonium(graph, pgrp, get_atom) -> str:
    """アンモニウムイオン: tetraRazanium  (Phase 146)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter
    n_idx = pgrp.atom_indices[0]
    c_neighbors = [nb for nb in graph.adjacency[n_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_neighbors:
        return "azanium"
    names = sorted(_name_carbon_substituent(graph, c, {n_idx}) for c in c_neighbors)
    counts = Counter(names)
    parts = []
    for sub in sorted(counts):
        n = counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    return "".join(parts) + "azanium"


def _name_organic_silane(graph, pgrp, get_atom) -> str:
    """有機シラン: {alkyl(s)}silane  (Phase 143)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    from collections import Counter
    si_idx = pgrp.atom_indices[0]
    c_neighbors = [nb for nb in graph.adjacency[si_idx] if get_atom(graph, nb).symbol == "C"]
    if not c_neighbors:
        return "silane"
    names = sorted(_name_carbon_substituent(graph, c, {si_idx}) for c in c_neighbors)
    counts = Counter(names)
    parts = []
    for sub in sorted(counts):
        n = counts[sub]
        mult = MULTIPLIER.get(n, "") if n > 1 else ""
        parts.append(f"{mult}{sub}")
    return "".join(parts) + "silane"


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
    n_c_nbrs = [nb for nb in graph.adjacency[n_idx]
                if nb != carbonyl_c and get_atom(graph, nb).symbol == "C"]
    if not n_c_nbrs:
        return base

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in n_c_nbrs]
    sub_counts = Counter(n_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        if cnt == 1:
            prefix_parts.append(f"N-{sub}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub}")
    n_prefix = "-".join(prefix_parts)
    return f"{alkyl_name} {n_prefix}carbamate"


def _name_peroxide(graph, pgrp, get_atom) -> str:
    """有機ペルオキシド命名: {dialkyl} peroxide (Phase 57)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    c1_idx = pgrp.atom_indices[0]
    o1_idx = pgrp.atom_indices[1]
    o2_idx = pgrp.atom_indices[2]
    c2_idx = pgrp.atom_indices[3]

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

    # 酸鎖 (カルボニル C から S を除外した方向)
    acid_chain = _collect_acid_chain(graph, carbonyl_c, {s_idx}, get_atom)
    n_acid = len(acid_chain)
    stem = CHAIN_PREFIX.get(n_acid, f"C{n_acid}")
    acid_name = f"{stem}anethioate"

    return f"S-{s_alkyl} {acid_name}"


def _name_disulfide(graph, pgrp, get_atom) -> str:
    """ジスルフィド命名: {dialkyl} disulfide (Phase 56)"""
    from .substituent import _name_carbon_substituent
    from .constants import MULTIPLIER
    s1_idx = pgrp.atom_indices[0]
    s2_idx = pgrp.atom_indices[1]
    c1_atoms = pgrp.atom_indices[2:]
    c2_atoms = []
    for nb in graph.adjacency[s2_idx]:
        if nb != s1_idx and get_atom(graph, nb).symbol == "C":
            c2_atoms.append(nb)
            break

    alkyl1 = _name_carbon_substituent(graph, c1_atoms[0], {s1_idx}) if c1_atoms else "methyl"
    alkyl2 = _name_carbon_substituent(graph, c2_atoms[0], {s2_idx}) if c2_atoms else "methyl"

    if alkyl1 == alkyl2:
        mult = MULTIPLIER.get(2, "di")
        return f"{mult}{alkyl1} disulfide"
    return f"{min(alkyl1, alkyl2)} {max(alkyl1, alkyl2)} disulfide"


def _name_nitroso(graph, pgrp, get_atom) -> str:
    """ニトロソ化合物: nitroso{alkane} (Phase 52)"""
    from .constants import CHAIN_PREFIX
    alkyl_c = pgrp.atom_indices[0]
    n_idx = pgrp.atom_indices[1]
    # 芳香環に直結している場合: nitroso + ring name
    c_atom = get_atom(graph, alkyl_c)
    if c_atom.is_aromatic and c_atom.in_ring:
        ring_atoms = set()
        for rt in (graph.ring_atom_sets or []):
            if alkyl_c in rt:
                ring_atoms.update(rt)
        if all(get_atom(graph, a).symbol == "C" for a in ring_atoms) and len(ring_atoms) == 6:
            return "nitrosobenzene"
        # ヘテロ環はフォールバック
        return "nitrosobenzene"
    # alkyl 側の炭素鎖を DFS 収集
    acid_chain = _collect_acid_chain(graph, alkyl_c, {n_idx}, get_atom)
    n_c = len(acid_chain)
    stem = CHAIN_PREFIX.get(n_c, f"C{n_c}")
    return f"nitroso{stem}ane"


def _name_azide(graph, pgrp, get_atom) -> str:
    """アジド: {alkyl} azide (Phase 53)"""
    from .substituent import _name_carbon_substituent
    alkyl_c = pgrp.atom_indices[0]
    n_idx = pgrp.atom_indices[1]
    alkyl_name = _name_carbon_substituent(graph, alkyl_c, {n_idx})
    return f"{alkyl_name} azide"


def _name_isocyanate(graph, pgrp, get_atom) -> str:
    """イソシアネート: {alkyl} isocyanate (Phase 54)"""
    from .substituent import _name_carbon_substituent
    alkyl_c = pgrp.atom_indices[0]
    n_idx = pgrp.atom_indices[1]
    alkyl_name = _name_carbon_substituent(graph, alkyl_c, {n_idx})
    return f"{alkyl_name} isocyanate"


def _name_isothiocyanate(graph, pgrp, get_atom) -> str:
    """イソチオシアネート: {alkyl} isothiocyanate (Phase 68)"""
    from .substituent import _name_carbon_substituent
    alkyl_c = pgrp.atom_indices[0]
    n_idx = pgrp.atom_indices[1]
    alkyl_name = _name_carbon_substituent(graph, alkyl_c, {n_idx})
    return f"{alkyl_name} isothiocyanate"


def _name_cyanate(graph, pgrp, get_atom) -> str:
    """シアン酸エステル: {alkyl} cyanate (Phase 69); HOCN → cyanic acid (Phase 128)"""
    from .substituent import _name_carbon_substituent
    cyano_c = pgrp.atom_indices[0]
    n_idx = pgrp.atom_indices[1] if len(pgrp.atom_indices) > 1 else None
    o_idx = pgrp.atom_indices[2] if len(pgrp.atom_indices) > 2 else None
    # O の隣の C が alkyl 部分
    if o_idx is not None:
        for nb in graph.adjacency[o_idx]:
            if nb != cyano_c and get_atom(graph, nb).symbol == "C":
                excluded = {o_idx} | ({n_idx} if n_idx else set())
                alkyl_name = _name_carbon_substituent(graph, nb, excluded)
                return f"{alkyl_name} cyanate"
        # O の隣に C がない (H のみ) → 遊離酸 HO-C≡N = cyanic acid
        return "cyanic acid"
    return "cyanate"


def _name_thiocyanate(graph, pgrp, get_atom) -> str:
    """チオシアン酸エステル: {alkyl} thiocyanate (Phase 69); HSC≡N → thiocyanic acid (Phase 128)"""
    from .substituent import _name_carbon_substituent
    cyano_c = pgrp.atom_indices[0]
    n_idx = pgrp.atom_indices[1] if len(pgrp.atom_indices) > 1 else None
    s_idx = pgrp.atom_indices[2] if len(pgrp.atom_indices) > 2 else None
    if s_idx is not None:
        for nb in graph.adjacency[s_idx]:
            if nb != cyano_c and get_atom(graph, nb).symbol == "C":
                excluded = {s_idx} | ({n_idx} if n_idx else set())
                alkyl_name = _name_carbon_substituent(graph, nb, excluded)
                return f"{alkyl_name} thiocyanate"
        # S の隣に C がない (H のみ) → 遊離酸 HS-C≡N = thiocyanic acid
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
        if cnt == 1:
            parts.append(f"N-{sub}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            parts.append(f"N,N-{mult}{sub}")
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

    if len(alkyl_names) == 2 and alkyl_names[0] == alkyl_names[1]:
        mult = MULTIPLIER.get(2, "di")
        return f"N,N'-{mult}{alkyl_names[0]}carbodiimide"

    names_sorted = sorted(alkyl_names)
    if len(names_sorted) == 2:
        return f"N-{names_sorted[0]}-N'-{names_sorted[1]}carbodiimide"
    return f"N-{names_sorted[0]}carbodiimide"


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
            if cnt == 1:
                parts.append(f"{tag}-{sub}")
            else:
                mult = MULTIPLIER.get(cnt, str(cnt))
                parts.append(f"{tag},{tag}-{mult}{sub}")
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
            if cnt == 1:
                parts.append(f"{tag}-{sub}")
            else:
                mult = MULTIPLIER.get(cnt, str(cnt))
                parts.append(f"{tag},{tag}-{mult}{sub}")
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

    if is_ald:
        return f"{stem}anal {suffix}"
    else:
        loc = chain.locant_map.get(hydrazone_c, 2)
        return f"{stem}an-{loc}-one {suffix}"


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

    if pgrp.group_type == "aldhydrazone":
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
        sub_parts.append(f"N-{mult}{sub_name}")
    n_prefix = ",".join(sub_parts) if sub_parts else ""

    if n_prefix:
        return f"{parent_name} {n_prefix}hydrazone"
    return f"{parent_name} hydrazone"


def _name_peroxyacid(graph, pgrp, get_atom) -> str:
    """
    ペルオキシ酸命名: {stem}aneperoxoic acid (Phase 77)
    例: CC(=O)OO → ethaneperoxoic acid
        CCC(=O)OO → propaneperoxoic acid
    """
    from .constants import CHAIN_PREFIX

    carbonyl_c = pgrp.atom_indices[0]
    # O1 (O-O-H の外側 O) と O2 を除外して炭素鎖を収集
    excluded: set[int] = {ai for ai in pgrp.atom_indices
                          if get_atom(graph, ai).symbol == "O"}

    acid_chain = _collect_acid_chain(graph, carbonyl_c, excluded, get_atom)
    stem = CHAIN_PREFIX.get(len(acid_chain), f"C{len(acid_chain)}")
    return f"{stem}aneperoxoic acid"


def _name_hydrazide(graph, pgrp, get_atom) -> str:
    """
    ヒドラジド命名: {stem}anehydrazide (Phase 75)
    例: CC(=O)NN → ethanehydrazide
        CCC(=O)NN → propanehydrazide
    """
    from .constants import CHAIN_PREFIX

    carbonyl_c = pgrp.atom_indices[0]
    n_idx = next((ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "N"), None)
    excluded: set[int] = set()
    if n_idx is not None:
        excluded.add(n_idx)
        # N-N の N2 も除外
        for nb in graph.adjacency[n_idx]:
            if nb != carbonyl_c and get_atom(graph, nb).symbol == "N":
                excluded.add(nb)

    acid_chain = _collect_acid_chain(graph, carbonyl_c, excluded, get_atom)
    stem = CHAIN_PREFIX.get(len(acid_chain), f"C{len(acid_chain)}")
    _ene_hy, _yne_hy = _chain_multiple_bonds(graph, acid_chain)
    if _ene_hy or _yne_hy:
        from .name_assembler import _format_multiple_bonds as _fmt_hy
        return f"{stem}{_fmt_hy(_ene_hy, _yne_hy)}ehydrazide"
    return f"{stem}anehydrazide"


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
                    if not has_cn_double:
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
        return f"1,2-{s1}-2-{s2}hydrazine"

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
        return f"azo{stem}ane"

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

    c_nbrs = [nb for nb in graph.adjacency[n_idx]
               if nb != carbonyl_c and get_atom(graph, nb).symbol == "C"]
    if not c_nbrs:
        return parent_name

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in c_nbrs]
    sub_counts = Counter(n_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        if cnt == 1:
            prefix_parts.append(f"N-{sub}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub}")
    prefix = "-".join(prefix_parts)
    return f"{prefix}{parent_name}"


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

    def _acid_stem_name(chain: list[int]) -> str:
        stem = CHAIN_PREFIX.get(len(chain), f"C{len(chain)}")
        _ene, _yne = _chain_multiple_bonds(graph, chain)
        if _ene or _yne:
            return f"{stem}{_fmt_anh(_ene, _yne)}oic"
        return f"{stem}anoic"

    acid1 = _acid_stem_name(chain1)
    acid2 = _acid_stem_name(chain2)

    acids = sorted([acid1, acid2])
    if acids[0] == acids[1]:
        return f"{acids[0]} anhydride"
    return f"{acids[0]} {acids[1]} anhydride"


def _name_cyclic(graph, find_rings, find_principal_ring,
                 collect_ring_substituents, assemble_ring_name,
                 detect_groups, principal_group, assign_stereochemistry,
                 get_atom) -> str:
    """環状化合物の命名フロー。"""
    # 官能基検出（環状分子用）
    groups = detect_groups(graph)
    pgrp = principal_group(groups)
    pgrp_type = pgrp.group_type if pgrp is not None else "alkane"

    # 環の検出
    rings = find_rings(graph)
    if not rings:
        # フォールバック: 純炭素環がない場合 (ヘテロ環のみ) は ring_atom_sets を使う
        if graph.ring_atom_sets:
            ring_list = sorted(graph.ring_atom_sets, key=len, reverse=True)
            rings = [list(ring_list[0])]
        else:
            raise ValueError("Ring detection failed despite in_ring=True atoms.")

    # Phase 40: alkene/alkyne で官能基アンカーが環外の場合 → alkane 扱いにして
    # 環外の C=C/C≡C を置換基として収集できるようにする
    if pgrp is not None and pgrp_type in ("alkene", "alkyne"):
        ring_atom_set = set()
        for r in rings:
            ring_atom_set.update(r)
        if not any(ai in ring_atom_set for ai in pgrp.atom_indices):
            pgrp = None
            pgrp_type = "alkane"

    # 主環の選択とロカント割り当て（主官能基アンカーを C1 にする）
    pgrp_atoms_list = pgrp.atom_indices if pgrp is not None else []
    ring_chain = find_principal_ring(graph, rings, pgrp_type, pgrp_atoms_list)

    # suffix ロカント (cyclohexanol, cyclohexanone, cyclohexanamine 等)
    from .constants import FUNCTIONAL_GROUPS as _FG
    suffix_locant: int | None = None
    if pgrp is not None and _FG.get(pgrp.group_type, _FG["alkane"]).needs_locant:
        for ai in pgrp.atom_indices:
            atom = get_atom(graph, ai)
            if atom.symbol == "C" and ai in ring_chain.locant_map:
                suffix_locant = ring_chain.locant_map[ai]
                break
    # 縮合多環 (naphthalene 等) の環外官能基ロカント (naphthalene-2-carboxylic acid 等)
    if (suffix_locant is None and pgrp is not None
            and ring_chain.base_name is not None):
        _spec_ex = _FG.get(pgrp.group_type)
        if _spec_ex is not None and _spec_ex.anchor_c1:
            _anchor_idx = pgrp.atom_indices[0] if pgrp.atom_indices else None
            if _anchor_idx is not None and _anchor_idx not in ring_chain.locant_map:
                for _nb in graph.adjacency[_anchor_idx]:
                    if _nb in ring_chain.locant_map:
                        suffix_locant = ring_chain.locant_map[_nb]
                        break

    # 置換基収集
    pgrp_atoms = pgrp.atom_indices if pgrp is not None else []
    substituents = collect_ring_substituents(graph, ring_chain, pgrp_atoms)

    # 立体化学
    from .chain_finder import PrincipalChain
    dummy_chain = PrincipalChain(
        atom_indices=ring_chain.ring_atoms,
        locant_map=ring_chain.locant_map,
    )
    stereo = assign_stereochemistry(graph, dummy_chain)

    # Phase 27/33: ポリオール / ジオン on ring (benzene 芳香環 or シクロアルカン)
    if pgrp is not None and pgrp_type in ("diol", "triol", "dione", "diamine", "triamine"):
        fg_c_atoms = {ai for ai in pgrp.atom_indices
                      if get_atom(graph, ai).symbol == "C"
                      and ai in ring_chain.locant_map}
        if fg_c_atoms:
            # 全回転 × 2方向で最小ロカントセットを探す
            ring_atoms = ring_chain.ring_atoms
            n_ring = len(ring_atoms)
            best_locs: tuple[int, ...] | None = None
            for start in range(n_ring):
                for direction in (1, -1):
                    if direction == 1:
                        order = [ring_atoms[(start + i) % n_ring] for i in range(n_ring)]
                    else:
                        order = [ring_atoms[(start - i) % n_ring] for i in range(n_ring)]
                    loc_map = {idx: i + 1 for i, idx in enumerate(order)}
                    locs = tuple(sorted(loc_map[a] for a in fg_c_atoms))
                    if best_locs is None or locs < best_locs:
                        best_locs = locs
            loc_str = ",".join(str(l) for l in (best_locs or ()))
            suffix_word = {
                "diol": "diol", "triol": "triol", "dione": "dione",
                "diamine": "diamine", "triamine": "triamine",
            }[pgrp_type]
            from .name_assembler import _build_prefix as _bp
            from .constants import CHAIN_PREFIX as _CP
            if ring_chain.is_aromatic and ring_chain.ring_size == 6:
                ring_base = "benzene"
            elif ring_chain.base_name is not None:
                ring_base = ring_chain.base_name
            else:
                _stem = _CP.get(ring_chain.ring_size, f"C{ring_chain.ring_size}")
                ring_base = f"cyclo{_stem}ane"
            prefix_part = _bp(substituents) if substituents else ""
            base_name = f"{ring_base}-{loc_str}-{suffix_word}"
            # 立体記述子 (Phase 66: ring diol/dione 立体化学)
            if stereo:
                combined = ",".join(d.strip("()") for d in stereo)
                stereo_str = f"({combined})-"
            else:
                stereo_str = ""
            full = f"{prefix_part}{base_name}" if prefix_part else base_name
            return f"{stereo_str}{full}"

    # N-置換アミド on 環 (benzamide → N-methylbenzamide など) の処理
    if pgrp is not None and pgrp_type == "amide":
        n_atom_idx = next(
            (ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "N"), None
        )
        if n_atom_idx is not None:
            h_count = sum(
                1 for nb in graph.adjacency[n_atom_idx]
                if get_atom(graph, nb).symbol == "H"
            )
            if h_count < 2:
                from .constants import MULTIPLIER
                from .substituent import _name_carbon_substituent
                from collections import Counter as _Counter
                carbonyl_c = next(
                    (ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "C"), None
                )
                base_amide_name = assemble_ring_name(
                    ring_chain=ring_chain,
                    substituents=substituents,
                    principal_grp_type=pgrp_type,
                    suffix_locant=suffix_locant,
                    stereo_descriptors=stereo,
                )
                excluded_n = {n_atom_idx}
                if carbonyl_c is not None:
                    excluded_n.add(carbonyl_c)
                c_nbrs_n = [
                    nb for nb in graph.adjacency[n_atom_idx]
                    if nb not in excluded_n and get_atom(graph, nb).symbol == "C"
                ]
                o_nbrs_n = [
                    nb for nb in graph.adjacency[n_atom_idx]
                    if get_atom(graph, nb).symbol == "O"
                ]
                n_hydroxy = [
                    "hydroxy" for o in o_nbrs_n
                    if any(get_atom(graph, nb2).symbol == "H"
                           for nb2 in graph.adjacency[o])
                ]
                if not c_nbrs_n and not n_hydroxy:
                    return base_amide_name
                from .name_assembler import _needs_bis_tris as _nbp2
                n_subs = [_name_carbon_substituent(graph, c, {n_atom_idx}) for c in c_nbrs_n]
                n_subs += n_hydroxy
                sub_counts = _Counter(n_subs)
                prefix_parts = []
                for sub in sorted(sub_counts):
                    cnt = sub_counts[sub]
                    sub_str = f"({sub})" if _nbp2(sub) else sub
                    if cnt == 1:
                        prefix_parts.append(f"N-{sub_str}")
                    else:
                        mult = MULTIPLIER.get(cnt, f"{cnt}")
                        prefix_parts.append(f"N,N-{mult}{sub_str}")
                prefix = "-".join(prefix_parts)
                return f"{prefix}{base_amide_name}"

    # 名前組み立て
    return assemble_ring_name(
        ring_chain=ring_chain,
        substituents=substituents,
        principal_grp_type=pgrp_type,
        suffix_locant=suffix_locant,
        stereo_descriptors=stereo,
    )


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
        if cnt == 1:
            prefix_parts.append(f"N-{sub}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub}")
    return "-".join(prefix_parts) + imine_base


def _name_secondary_tertiary_amide(graph, carbonyl_c: int, n_idx: int, get_atom) -> str:
    """
    二級・三級アミドの命名: N-alkyl-alkanamide 形式。
    例: CC(=O)NC → N-methylethanamide
        CC(=O)N(C)C → N,N-dimethylethanamide
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

    if not c_nbrs and not n_hydroxy_subs:
        return parent_name

    n_subs = [_name_carbon_substituent(graph, c, {n_idx}) for c in c_nbrs]
    n_subs += n_hydroxy_subs
    sub_counts = Counter(n_subs)
    prefix_parts = []
    from .name_assembler import _needs_bis_tris as _nbp
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        sub_str = f"({sub})" if _nbp(sub) else sub
        if cnt == 1:
            prefix_parts.append(f"N-{sub_str}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub_str}")

    prefix = "-".join(prefix_parts)
    return f"{prefix}{parent_name}"


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

        if not n_subs:
            return parent_name

        sub_counts = Counter(n_subs)
        prefix_parts = []
        for sub in sorted(sub_counts):
            cnt = sub_counts[sub]
            if cnt == 1:
                prefix_parts.append(f"N-{sub}")
            else:
                mult = MULTIPLIER.get(cnt, f"{cnt}")
                prefix_parts.append(f"N,N-{mult}{sub}")
        prefix = "-".join(prefix_parts)
        return f"{prefix}{parent_name}"

    # 非環アミン: 非環炭素のみでカウント（環境界で停止）
    def chain_c_atoms(start_c: int) -> list[int]:
        if get_atom(graph, start_c).in_ring:
            return []
        return _collect_acid_chain(graph, start_c, {n_idx}, get_atom)

    chains = [chain_c_atoms(c) for c in c_neighbors]
    lengths = [len(ch) for ch in chains]
    max_len = max(lengths) if lengths else 1

    # 最長鎖を親鎖とする（最初に見つかったものを選択）
    parent_pos = next((i for i, l in enumerate(lengths) if l == max_len), 0)
    parent_chain = chains[parent_pos]
    stem = CHAIN_PREFIX.get(max_len, f"C{max_len}")

    # 多重結合の検出: ene/yne がある場合は親鎖の locant を明示
    _ene_am, _yne_am = _chain_multiple_bonds(graph, parent_chain)
    if _ene_am or _yne_am:
        from .name_assembler import _format_multiple_bonds as _fmt_am
        parent_name = f"{stem}{_fmt_am(_ene_am, _yne_am)}-1-amine"
    elif max_len >= 3:
        parent_name = f"{stem}an-1-amine"
    else:
        parent_name = f"{stem}anamine"

    # N 上の置換基（親鎖以外）
    n_subs = [
        _name_carbon_substituent(graph, c, {n_idx})
        for i, c in enumerate(c_neighbors) if i != parent_pos
    ]

    if not n_subs:
        return parent_name

    sub_counts = Counter(n_subs)
    prefix_parts = []
    for sub in sorted(sub_counts):
        cnt = sub_counts[sub]
        if cnt == 1:
            prefix_parts.append(f"N-{sub}")
        else:
            mult = MULTIPLIER.get(cnt, f"{cnt}")
            prefix_parts.append(f"N,N-{mult}{sub}")

    prefix = "-".join(prefix_parts)
    return f"{prefix}{parent_name}"


def _name_acyclic(graph, detect_groups, principal_group,
                  find_principal_chain, get_multiple_bond_locants,
                  collect_substituents, assign_stereochemistry,
                  assemble_name, get_atom) -> str:
    """非環状化合物の命名フロー（Phase 1-2）。"""
    from .functional_group import FunctionalGroup

    # 官能基検出
    groups = detect_groups(graph)
    pgrp: FunctionalGroup | None = principal_group(groups)

    # 二級・三級アミドは専用パスで処理
    if pgrp is not None and pgrp.group_type == "amide":
        n_atom_idx = next((ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "N"), None)
        if n_atom_idx is not None:
            h_count = sum(1 for nb in graph.adjacency[n_atom_idx]
                          if get_atom(graph, nb).symbol == "H")
            if h_count < 2:
                carbonyl_c = pgrp.atom_indices[0]
                return _name_secondary_tertiary_amide(graph, carbonyl_c, n_atom_idx, get_atom)

    # 二級・三級アミンは専用パスで処理
    if pgrp is not None and pgrp.group_type == "amine":
        if get_atom(graph, pgrp.atom_indices[0]).symbol == "N":
            n_idx = pgrp.atom_indices[0]
            c_nbrs = pgrp.atom_indices[1:]
            return _name_secondary_tertiary_amine(graph, n_idx, c_nbrs, get_atom)

    # N-置換イミン (Phase 64)
    if pgrp is not None and pgrp.group_type == "imine":
        n_atom_idx = next((ai for ai in pgrp.atom_indices if get_atom(graph, ai).symbol == "N"), None)
        if n_atom_idx is not None:
            imine_c = pgrp.atom_indices[0]
            n_c_subs = [nb for nb in graph.adjacency[n_atom_idx]
                        if nb != imine_c and get_atom(graph, nb).symbol == "C"]
            if n_c_subs:
                return _name_n_substituted_imine(
                    graph, imine_c, n_atom_idx, n_c_subs,
                    find_principal_chain, get_multiple_bond_locants,
                    collect_substituents, assemble_name, get_atom,
                )

    # 主鎖探索・ロカント割り当て
    chain = find_principal_chain(graph, pgrp)

    # 多重結合ロカント
    mb = get_multiple_bond_locants(graph, chain)

    # suffix ロカント（複数官能基対応）
    from .constants import FUNCTIONAL_GROUPS as _FG
    suffix_locants_list: list[int] = []
    if pgrp is not None and _FG.get(pgrp.group_type, _FG["alkane"]).needs_locant:
        if pgrp.group_type in ("diol", "triol"):
            # ジェミナルジオール対応: 各 O から親 C のロカントを収集 (重複許可)
            for ai in pgrp.atom_indices:
                atom = get_atom(graph, ai)
                if atom.symbol == "O":
                    for nb in graph.adjacency[ai]:
                        if get_atom(graph, nb).symbol == "C" and nb in chain.locant_map:
                            suffix_locants_list.append(chain.locant_map[nb])
                            break
        else:
            for ai in pgrp.atom_indices:
                atom = get_atom(graph, ai)
                if atom.symbol == "C" and ai in chain.locant_map:
                    loc = chain.locant_map[ai]
                    if loc not in suffix_locants_list:
                        suffix_locants_list.append(loc)
        suffix_locants_list.sort()

    suffix_locant = suffix_locants_list[0] if suffix_locants_list else None
    suffix_locants = suffix_locants_list if len(suffix_locants_list) > 1 else None

    # 置換基収集
    pgrp_atom_indices = pgrp.atom_indices if pgrp is not None else []
    substituents = collect_substituents(
        graph, chain.atom_indices, chain.locant_map, pgrp_atom_indices,
    )

    # 立体化学
    stereo = assign_stereochemistry(graph, chain)

    # 名前組み立て
    pgrp_type = pgrp.group_type if pgrp is not None else "alkane"
    return assemble_name(
        chain_length=chain.length,
        principal_group_type=pgrp_type,
        multiple_bonds=mb,
        substituents=substituents,
        stereo_descriptors=stereo,
        suffix_locant=suffix_locant,
        suffix_locants=suffix_locants,
    )
