"""
IUPAC 名の組み立て。

主鎖の語幹・suffix・ロカント・置換基接頭辞を組み合わせて
最終的な IUPAC 名文字列を生成する。

IUPAC 2013 Blue Book P-31 / P-44 / P-65 の規則に従う。

ロカント省略ルール:
  - ethene (C=C, 2C): ロカント省略 → "ethene" (not "eth-1-ene")
  - ethyne (C#C, 2C): ロカント省略 → "ethyne" (not "eth-1-yne")
  - ethanal (CC=O): C1固定、ロカント省略 → "ethanal"
  - ethanoic acid: C1固定、ロカント省略 → "ethanoic acid"
  - ethan-1-ol: IUPAC 2013 では常にロカント表記 → "ethan-1-ol"
  - propan-2-one: ロカント表記 → "propan-2-one"
"""

from __future__ import annotations

import re as _re

from .constants import CHAIN_PREFIX, SUFFIX_MAP, MULTIPLIER

# 複合置換基（ロカント・コンマ・既存倍数詞を含む）用の倍数詞
_BIS_MULTIPLIER: dict[int, str] = {
    2: "bis",
    3: "tris",
    4: "tetrakis",
    5: "pentakis",
    6: "hexakis",
    7: "heptakis",
    8: "octakis",
}


def _needs_bis_tris(name: str) -> bool:
    """
    置換基名が bis/tris を必要とするか判定する。

    以下のいずれかに該当すると True:
      - 名前に数字（ロカント）を含む: "1-methylethyl" → bis(1-methylethyl)
      - 名前がすでに倍数詞で始まる: "dimethyl..." → bis(dimethyl...)
      - アリール基を含む複合置換基: "phenylmethyl" → (phenylmethyl)
    """
    # Already parenthesized (e.g., '(carboxyoxy)') — no additional wrapping needed
    if name.startswith("(") and name.endswith(")"):
        return False
    if _re.search(r"[0-9]", name):
        return True
    for mult in ("di", "tri", "tetra", "penta", "hexa", "hepta", "octa",
                 "nona", "deca", "bis", "tris"):
        if name.startswith(mult) and len(name) > len(mult):
            return True
    # アリール含有複合置換基 (phenylmethyl 等) は括弧が必要
    for aryl in ("phenyl", "naphthyl", "furyl", "thienyl", "pyridyl", "indolyl",
                 "quinolyl", "benzofuryl"):
        if aryl in name and name != aryl:
            return True
    # carboxy/hydroxy/oxo 含有複合置換基 (carboxymethyl 等) は括弧が必要
    for pfx in ("carboxy", "hydroxy", "oxo", "amino", "imino", "nitroso"):
        if pfx in name and name != pfx:
            return True
    # alkylsulfinyl/alkylsulfonyl/alkylsulfanyl 等 — 常に括弧が必要
    for sfx in ("sulfinyl", "sulfonyl", "sulfanyl"):
        if name.endswith(sfx) and name != sfx:
            return True
    return False


def assemble_name(
    chain_length: int,
    principal_group_type: str,
    multiple_bonds: dict[str, list[int]],
    substituents: list[tuple[int, str]],
    stereo_descriptors: list[str],
    suffix_locant: int | None = None,
    suffix_locants: list[int] | None = None,
) -> str:
    """
    IUPAC 名を組み立てる。

    Args:
        chain_length: 主鎖炭素数
        principal_group_type: 'alkane', 'alcohol', 'ketone', ...
        multiple_bonds: {'ene': [ロカントリスト], 'yne': [...]}
        substituents: [(ロカント, 置換基名)] のリスト
        stereo_descriptors: 立体記述子のリスト
        suffix_locant: アルコール (-ol) やケトン (-one) の官能基ロカント

    Returns:
        IUPAC 名文字列
    """
    stem = CHAIN_PREFIX.get(chain_length)
    if stem is None:
        raise ValueError(f"Unsupported chain length: {chain_length}")

    suffix = SUFFIX_MAP[principal_group_type]

    # ─── Phase 120/128: 保留名の早期返却 ───────────────────────────
    # acetophenone: 2炭素ケトン + C1 に phenyl 置換基 1個 (IUPAC 2013 PIN)
    if (principal_group_type == "ketone"
            and chain_length == 2
            and (suffix_locant is None or suffix_locant == 1)
            and len(substituents) == 1
            and substituents[0][1] == "phenyl"):
        return "acetophenone"

    # guanidine: アミジン 1C + 'amino' 置換基 (IUPAC 2013 PIN P-66.4.1.3)
    if (principal_group_type == "amidine"
            and chain_length == 1
            and len(substituents) == 1
            and substituents[0][1] == "amino"):
        return "guanidine"

    # ─── 1. 接頭辞部分（substituents を整理）───────────────────────
    # 1炭素鎖ではロカントは常に 1 (自明) → 省略 ("phenyl" not "1-phenyl")
    effective_subs = substituents
    if chain_length == 1:
        effective_subs = [(None, name) for _, name in substituents]
    # Phase 35: エテン (2炭素鎖 alkene) の単一置換基はロカント省略
    # 例: ClCH=CH₂ → "chloroethene" (not "1-chloroethene")
    # 条件: chain_length=2, ene suffix, 置換基が1種1個のみ
    elif (chain_length == 2
          and principal_group_type == "alkene"
          and len(substituents) == 1):
        effective_subs = [(None, name) for _, name in substituents]
    # Phase 114: エタン (2炭素鎖 alkane) の単一置換基はロカント省略
    # エタンは対称なので C1=C2 → "chloroethane" (not "1-chloroethane")
    elif (chain_length == 2
          and principal_group_type == "alkane"
          and len(substituents) == 1):
        effective_subs = [(None, name) for _, name in substituents]
    # Phase 213: 2炭素鎖 alkane で全6位置が同一置換基 → ロカント省略 (hexafluoroethane 等)
    elif (chain_length == 2
          and principal_group_type == "alkane"
          and len(substituents) == 6
          and len({n for _, n in substituents}) == 1):
        effective_subs = [(None, name) for _, name in substituents]
    prefix_part = _build_prefix(effective_subs)

    # ─── 2. 語幹 + suffix の組み立て ───────────────────────────────
    name_body = _build_name_body(stem, suffix, multiple_bonds, suffix_locant, chain_length,
                                suffix_locants=suffix_locants,
                                has_substituents=bool(effective_subs))

    # ─── 3. 立体記述子 ─────────────────────────────────────────────
    # 複数不斉中心: (1R)-(2S) ではなく (1R,2S) にまとめる
    if stereo_descriptors:
        combined = ",".join(d.strip("()") for d in stereo_descriptors)
        stereo_part = f"({combined})"
    else:
        stereo_part = ""

    # ─── 4. 結合 ───────────────────────────────────────────────────
    # IUPAC: 接頭辞と親鎖名は直接結合（ハイフンなし）
    # 例: "2-methyl" + "propane" → "2-methylpropane"
    #     "4-chloro-2-methyl" + "hexane" → "4-chloro-2-methylhexane"
    if prefix_part:
        result = f"{prefix_part}{name_body}"
    else:
        result = name_body

    if stereo_part:
        result = f"{stereo_part}-{result}"

    return result


# ─── 接頭辞構築 ─────────────────────────────────────────────────────

def _build_prefix(substituents: list[tuple[int, str]]) -> str:
    """
    置換基リストを IUPAC 形式の接頭辞文字列に変換する。

    - アルファベット順にソート（倍数詞は無視）
    - 同じ置換基は倍数詞でまとめる
    - 例: [(2, 'methyl'), (4, 'chloro')] → '4-chloro-2-methyl'
    """
    if not substituents:
        return ""

    # 置換基名でグループ化
    by_name: dict[str, list[int]] = {}
    for loc, name in substituents:
        by_name.setdefault(name, []).append(loc)

    def alpha_key(name: str) -> str:
        """倍数詞を除いたベース名でアルファベットソート。"""
        for mult in ["di", "tri", "tetra", "penta", "hexa", "hepta", "octa"]:
            if name.startswith(mult):
                return name[len(mult):]
        return name

    sorted_names = sorted(by_name.keys(), key=alpha_key)

    # ロカントが全て None か (chain_length=1 のメタノン等)
    all_loc_none = all(loc is None for locs in by_name.values() for loc in locs)

    parts: list[str] = []
    prev_parens = False  # 直前の部品が括弧付きだったか (no-locant モード用)
    for name in sorted_names:
        locs = sorted(l for l in by_name[name] if l is not None)
        n_total = len(by_name[name])
        if not locs:
            # ロカントなし (1炭素鎖等)
            n = n_total
            if n > 1 and _needs_bis_tris(name):
                mult = _BIS_MULTIPLIER.get(n, f"({n}×)")
                parts.append(f"{mult}({name})")
                prev_parens = True
            elif n > 1:
                mult = MULTIPLIER.get(n, "")
                if prev_parens:
                    parts.append(f"{mult}({name})")
                    prev_parens = True
                else:
                    parts.append(f"{mult}{name}")
                    prev_parens = False
            elif _needs_bis_tris(name) or prev_parens:
                # 複合名 or 直前が括弧付き → 括弧で囲む (IUPAC P-14.5.7)
                parts.append(f"({name})")
                prev_parens = True
            else:
                parts.append(name)
                prev_parens = False
            continue
        prev_parens = False
        n = len(locs)
        loc_str = ",".join(str(l) for l in locs)
        if n > 1 and _needs_bis_tris(name):
            # 複合置換基が複数: bis/tris + 括弧
            mult = _BIS_MULTIPLIER.get(n, f"({n}×)")
            parts.append(f"{loc_str}-{mult}({name})")
        elif _needs_bis_tris(name):
            # 複合置換基が 1 つ: 括弧 (IUPAC P-14.5.7)
            mult = MULTIPLIER.get(n, "")
            parts.append(f"{loc_str}-{mult}({name})")
        else:
            mult = MULTIPLIER.get(n, "")
            parts.append(f"{loc_str}-{mult}{name}")

    # no-locant モード: 括弧がセパレータ役を担うのでハイフン不要
    sep = "" if all_loc_none else "-"
    return sep.join(parts)


# ─── 語幹 + suffix 構築 ─────────────────────────────────────────────

def _build_name_body(
    stem: str,
    suffix: str,
    multiple_bonds: dict[str, list[int]],
    suffix_locant: int | None,
    chain_length: int,
    suffix_locants: list[int] | None = None,
    has_substituents: bool = False,
) -> str:
    """
    語幹と suffix を組み合わせて名前の本体部分を作る。

    命名例:
      hexane         → 'hexane'
      hex-1-ene      → 'hex-1-ene'
      hexan-1-ol     → 'hexan-1-ol'
      propan-2-one   → 'propan-2-one'
      ethanal        → 'ethanal'
      ethanoic acid  → 'ethanoic acid'
    """
    ene_locs = multiple_bonds.get("ene", [])
    yne_locs = multiple_bonds.get("yne", [])

    has_ene = bool(ene_locs)
    has_yne = bool(yne_locs)
    has_multiple_bond = has_ene or has_yne

    if suffix == "ane":
        # アルカン: ロカントなし、多重結合なし
        return f"{stem}ane"

    if suffix == "ene":
        # アルケン (enyne も含む)
        if yne_locs:
            # en-yne 化合物: _format_multiple_bonds で en+yn を組み合わせ、末尾 'e' を付加
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}e"
        if chain_length == 2:
            # ethene: ロカント省略
            return f"{stem}ene"
        n_ene = len(ene_locs)
        loc_str = _locant_list_str(ene_locs)
        if n_ene == 1:
            return f"{stem}-{loc_str}-ene"
        # 2以上の C=C → diene/triene (allene, conjugated diene 等)
        # IUPAC: prop/but/pent... + "a" + "-locs-" + di/tri + "ene"
        _IENE_MULT = {2: "di", 3: "tri", 4: "tetra", 5: "penta"}
        mult = _IENE_MULT.get(n_ene, f"{n_ene}")
        return f"{stem}a-{loc_str}-{mult}ene"

    if suffix == "yne":
        # アルキン (diyne/triyne 対応)
        if chain_length == 2:
            return f"{stem}yne"
        n_yne = len(yne_locs)
        loc_str = _locant_list_str(yne_locs)
        if n_yne == 1:
            return f"{stem}-{loc_str}-yne"
        mult = _BOND_MULT.get(n_yne, f"{n_yne}")
        return f"{stem}a-{loc_str}-{mult}yne"

    if suffix == "al":
        # アルデヒド: C1固定、ロカント省略
        # ethanal, propanal, butanal ...
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}al"
        return f"{stem}anal"

    if suffix == "dial":
        # ジアルデヒド: C1,Cn 固定、ロカント省略 (propanedial 等)
        # "dial" は子音始まり → 'e' 保持
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}edial"
        return f"{stem}anedial"

    if suffix == "oic acid":
        # カルボン酸: C1固定、ロカント省略
        # Phase 120: 1炭素鎖は "formic acid", 2炭素鎖は "acetic acid" (IUPAC 2013 PIN)
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}oic acid"
        if chain_length == 1:
            return "formic acid"
        if chain_length == 2:
            return "acetic acid"
        return f"{stem}anoic acid"

    if suffix == "ol":
        # アルコール: 1炭素鎖はロカント省略 (methanol 等)
        loc = suffix_locant if suffix_locant is not None else 1
        if has_multiple_bond:
            # Phase 34: 2炭素鎖 (ethenol, ethynol) はロカント省略
            if chain_length == 2:
                if has_ene and not has_yne:
                    return f"{stem}enol"
                if has_yne and not has_ene:
                    return f"{stem}ynol"
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}-{loc}-ol"
        if chain_length <= 2 and loc == 1:
            # Phase 117: 2炭素鎖loc=1は "ethanol" (IUPAC 2013 保留名 PIN)
            return f"{stem}anol"
        return f"{stem}an-{loc}-ol"

    if suffix == "one":
        # ケトン: 1炭素鎖はロカント省略 (diphenylmethanone)
        if chain_length == 1:
            return f"{stem}anone"
        loc = suffix_locant if suffix_locant is not None else 2
        if has_multiple_bond:
            # Phase 129: 2炭素鎖 ene-one → "ethen-{loc}-one" (ene ロカント省略)
            if chain_length == 2 and ene_locs == [1] and not yne_locs:
                if loc == 1:
                    return f"{stem}enone"   # ethenone (locant 1 unambiguous)
                return f"{stem}en-{loc}-one"
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}-{loc}-one"
        return f"{stem}an-{loc}-one"

    if suffix == "amine":
        # アミン: ロカント表記 (propan-2-amine 等); 1-2炭素鎖 C1 はロカント省略
        loc = suffix_locant if suffix_locant is not None else 1
        if has_multiple_bond:
            # 2炭素 ene-amine: ロカント省略 → "ethenamine"
            if chain_length == 2 and ene_locs == [1] and not yne_locs and loc == 1:
                return f"{stem}enamine"
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}-{loc}-amine"
        if chain_length <= 2 and loc == 1:
            return f"{stem}anamine"
        return f"{stem}an-{loc}-amine"

    if suffix == "imine":
        loc = suffix_locant if suffix_locant is not None else 1
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}-{loc}-imine"
        if chain_length <= 2 and loc == 1:
            return f"{stem}animine"
        return f"{stem}an-{loc}-imine"

    if suffix == "one oxime":
        # ケトオキシム: ロカント表記 (propan-2-one oxime 等)
        loc = suffix_locant if suffix_locant is not None else 2
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}-{loc}-one oxime"
        return f"{stem}an-{loc}-one oxime"

    if suffix == "al oxime":
        # アルドキシム: C1 固定、ロカント省略 (ethanal oxime 等)
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}al oxime"
        return f"{stem}anal oxime"

    if suffix == "nitrile":
        # ニトリル: ニトリル C は常に C1 → ロカント省略
        # "nitrile" は子音始まりのため "en/yn" の後に 'e' を保持する
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}enitrile"
        return f"{stem}anenitrile"

    if suffix == "dinitrile":
        # ジニトリル: C1 と Cn に固定、ロカント省略 (butanedinitrile)
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}edinitrile"
        return f"{stem}anedinitrile"

    if suffix == "thiol":
        # チオール: -thiol は子音始まりのため e を保持; 1-2炭素鎖 C1 はロカント省略
        loc = suffix_locant if suffix_locant is not None else 1
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}e-{loc}-thiol"
        if chain_length <= 2 and loc == 1:
            return f"{stem}anethiol"
        return f"{stem}ane-{loc}-thiol"

    if suffix == "selenol":
        loc = suffix_locant if suffix_locant is not None else 1
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}e-{loc}-selenol"
        if chain_length <= 2 and loc == 1:
            return f"{stem}aneselenol"
        return f"{stem}ane-{loc}-selenol"

    if suffix == "tellurol":
        loc = suffix_locant if suffix_locant is not None else 1
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}e-{loc}-tellurol"
        if chain_length <= 2 and loc == 1:
            return f"{stem}anetellurol"
        return f"{stem}ane-{loc}-tellurol"

    if suffix == "dithiol":
        # ジチオール: -dithiol は子音始まりのため e を保持 (ethane-1,2-dithiol)
        if chain_length == 1:
            return f"{stem}anedithiol"
        locs = sorted(suffix_locants) if suffix_locants else [1, 2]
        loc_str = _locant_list_str(locs)
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}e-{loc_str}-dithiol"
        return f"{stem}ane-{loc_str}-dithiol"

    if suffix == "thial":
        # チオアルデヒド: C1固定、ロカント省略 (ethanethial, propanethial)
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}ethial"
        return f"{stem}anethial"

    if suffix == "thione":
        # チオケトン: IUPAC 2013 — ロカント前の e を省略 (propan-2-thione)
        loc = suffix_locant if suffix_locant is not None else 2
        if has_multiple_bond:
            if chain_length == 2 and ene_locs == [1] and not yne_locs:
                return f"{stem}ene-{loc}-thione"
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}e-{loc}-thione"
        return f"{stem}an-{loc}-thione"

    # ─── 複数官能基 suffix ────────────────────────────────────────────

    if suffix == "dioic acid":
        # "dioic acid" は子音始まり → "en/yn" の後に 'e' を保持する
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}edioic acid"
        return f"{stem}anedioic acid"

    if suffix == "diol":
        if chain_length == 1:
            return f"{stem}anediol"
        locs = sorted(suffix_locants) if suffix_locants else [1, 2]
        loc_str = _locant_list_str(locs)
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}e-{loc_str}-diol"
        return f"{stem}ane-{loc_str}-diol"

    if suffix == "triol":
        locs = sorted(suffix_locants) if suffix_locants else [1, 2, 3]
        loc_str = _locant_list_str(locs)
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}e-{loc_str}-triol"
        return f"{stem}ane-{loc_str}-triol"

    if suffix in ("dione", "trione"):
        # -dione/-trione は子音始まりのため e を保持 (pentane-2,4-dione 等)
        default_locs = [2, 4] if suffix == "dione" else [2, 3, 4]
        locs = sorted(suffix_locants) if suffix_locants else default_locs
        loc_str = _locant_list_str(locs)
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}e-{loc_str}-{suffix}"
        return f"{stem}ane-{loc_str}-{suffix}"

    if suffix in ("diamine", "triamine"):
        # ジアミン / トリアミン: 各 N のロカントを付加
        if chain_length == 1:
            return f"{stem}ane{suffix}"
        locs = sorted(suffix_locants) if suffix_locants else [1, 2]
        loc_str = _locant_list_str(locs)
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}-{loc_str}-{suffix}"
        return f"{stem}ane-{loc_str}-{suffix}"

    if suffix == "amide":
        # C1 固定、ロカント省略; Phase 120: 1C→formamide, 2C→acetamide (IUPAC 2013 PIN)
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}amide"
        if chain_length == 1:
            return "formamide"
        if chain_length == 2:
            return "acetamide"
        return f"{stem}anamide"

    if suffix == "diamide":
        # ジアミド: C1, Cn 固定、ロカント省略 (ethanediamide = oxamide)
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}ediamide"
        return f"{stem}anediamide"

    if suffix == "imidamide":
        # Phase 70: アミジン C1 固定 (ethanimidamide)
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}imidamide"
        return f"{stem}animidamide"

    if suffix == "thioamide":
        # Phase 41: チオアミド C1 固定 (ethanethioamide)
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}thioamide"
        return f"{stem}anethioamide"

    if suffix == "one hydrazone":
        # Phase 43: ケトヒドラゾン (propan-2-one hydrazone 等)
        loc = suffix_locant if suffix_locant is not None else 2
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}-{loc}-one hydrazone"
        return f"{stem}an-{loc}-one hydrazone"

    if suffix == "al hydrazone":
        # Phase 43: アルドヒドラゾン (ethanal hydrazone 等)
        if has_multiple_bond:
            mb = _format_multiple_bonds(ene_locs, yne_locs)
            return f"{stem}{mb}al hydrazone"
        return f"{stem}anal hydrazone"

    # フォールバック
    return f"{stem}{suffix}"


_BOND_MULT = {1: "", 2: "di", 3: "tri", 4: "tetra", 5: "penta"}


def _format_multiple_bonds(ene_locs: list[int], yne_locs: list[int]) -> str:
    """
    多重結合の中間部分を構築する。suffix の直前に挿入する部分。
    例: ene=[1]      → '-1-en'
        ene=[3,5]    → 'a-3,5-dien'   (diene: leading '-' → 'a')
        ene=[1], yne=[3] → '-1-en-3-yn'
        ene=[1,3], yne=[5] → 'a-1,3-dien-5-yn'
    """
    needs_a = len(ene_locs) > 1 or len(yne_locs) > 1
    parts = []
    if ene_locs:
        loc_str = _locant_list_str(ene_locs)
        mult = _BOND_MULT.get(len(ene_locs), f"{len(ene_locs)}")
        parts.append(f"-{loc_str}-{mult}en")
    if yne_locs:
        loc_str = _locant_list_str(yne_locs)
        mult = _BOND_MULT.get(len(yne_locs), f"{len(yne_locs)}")
        parts.append(f"-{loc_str}-{mult}yn")
    result = "".join(parts)
    if needs_a and result.startswith("-"):
        result = "a" + result  # "a" + "-2,4-dien" → "a-2,4-dien"; stem+"a-2,4-dien" → "penta-2,4-dien"
    return result


def _locant_list_str(locants: list[int]) -> str:
    """ロカントリストを '1' or '1,3' などの文字列に変換。"""
    return ",".join(str(l) for l in sorted(locants))


def fix_enclosing_marks(name: str) -> str:
    """Apply IUPAC 2013 P-16.3.4 alternating enclosing marks.

    () groups whose content contains () are upgraded to [].
    [] groups whose content contains [] are upgraded to {}.
    Simple groups like '(E)-', '(1R,2S)' are left unchanged.
    """

    def _upgrade(s: str, open_in: str, close_in: str,
                 open_out: str, close_out: str) -> str:
        result: list[str] = []
        i = 0
        while i < len(s):
            if s[i] != open_in:
                result.append(s[i])
                i += 1
                continue
            depth = 1
            j = i + 1
            while j < len(s) and depth > 0:
                if s[j] == open_in:
                    depth += 1
                elif s[j] == close_in:
                    depth -= 1
                j += 1
            content = s[i + 1: j - 1]
            if open_in in content:
                result.append(open_out + content + close_out)
            else:
                result.append(open_in + content + close_in)
            i = j
        return "".join(result)

    name = _upgrade(name, "(", ")", "[", "]")
    name = _upgrade(name, "[", "]", "{", "}")
    return name
