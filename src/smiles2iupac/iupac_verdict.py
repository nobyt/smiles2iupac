"""IUPAC 2013 規則に基づく命名差異の分類ロジック。

本モジュールは PubChem IUPAC名と本ライブラリの出力名を比較し、
どちらが IUPAC 2013 Blue Book に適合しているかを判定する。

Verdict 一覧:
  pubchem_wrong    — PubChem が IUPAC 2013 に従っていない (我々の名前が PIN)
  our_retained     — 我々の名前が保留名 (保留); PubChem が PIN
  pubchem_retained — PubChem が保留名を使用; 我々の名前が系統的 PIN
  both_pin         — 両方が有効な PIN (IUPAC 2013 が複数を許容)
  tautomer         — 互変異性体の命名差 (どちらも正当)
  needs_review     — 自動分類不能 (手動確認が必要)
"""
from __future__ import annotations

import re
from typing import Optional

# ---------------------------------------------------------------------------
# 保留名データ (IUPAC 2013 P-31.1.3 / P-65.1.2 / P-102.4 等)
# key: 保留名 (小文字), value: None (PIN は PubChem 名を使用)
# ---------------------------------------------------------------------------

# 単独の化合物名が保留名
_RETAINED_EXACT: frozenset[str] = frozenset({
    # ケトン類 (P-65.1.2.2)
    "acetophenone",      # PIN: 1-phenylethanone
    "benzophenone",      # PIN: diphenylmethanone
    "acetone",           # PIN: propan-2-one
    # アゾ・アミン類 (P-62.4)
    "azobenzene",        # PIN: diphenyldiazene
    "diphenylamine",     # PIN: N-phenylaniline
    "triphenylamine",    # PIN: N,N-diphenylaniline
    # 縮合環 (P-31.1.3.4 / P-31.1.3.2)
    "biphenyl",          # PIN: 1,1'-biphenyl
    "indoline",          # PIN: 2,3-dihydro-1H-indole
    "phthalimide",       # PIN: isoindole-1,3-dione
    "fluorene",          # PIN: 9H-fluorene
    "chromane",          # PIN: 3,4-dihydro-2H-chromene
    "coumarin",          # PIN: 2H-chromen-2-one
    "xanthene",          # PIN: 9H-xanthene
    "thioxanthene",      # PIN: 9H-thioxanthene
    "carbazole",         # PIN: 9H-carbazole
    "phenoxazine",       # PIN: 10H-phenoxazine
    "phenothiazine",     # PIN: 10H-phenothiazine
    # プリン塩基 (P-31.1.3.4)
    "xanthine",          # PIN: 3,7-dihydropurine-2,6-dione
    "adenine",
    "guanine",
    "caffeine",
    "theophylline",
    # その他の慣用名
    "acetaldehyde",      # PIN: ethanal
    "acrolein",          # PIN: prop-2-enal
    "ethanohydrazide",   # 慣用名 (acetohydrazide)
})

# 化合物名中に現れる保留親骨格 (部分一致)
_RETAINED_PARENTS: tuple[str, ...] = (
    "biphenyl",
    "indoline",
    "chromane",
    "coumarin",
    "xanthene",
    "thioxanthene",
    "fluorene",
    "phthalimide",
    "phenoxazine",
    "phenothiazine",
    "carbazole",
)

# 窒素複素環の エノール(-ol) vs ケト(-one) 互変異性体パターン
# エノール形の親骨格名 (小文字, 部分一致用)
_HETERO_OL_PARENTS = re.compile(
    r'\b(?:pyridin|quinolin|isoquinolin|phthalazin|pyrimidin|quinazolin'
    r'|quinoxalin|pyridazin|naphthyridin|acridin|phenanthridin'
    r'|cinnolin|perimidin|pteridin)\w*-\d+(?:\(\d+H\))?-ol\b',
    re.IGNORECASE,
)

# PubChem が保留縮合環成分を使う場合のパターン
# 例: benzo[f]benzimidazole (benzimidazole は保留), naphtho[2,3-d]imidazole は系統的
_PUBCHEM_RETAINED_FUSED = re.compile(
    r'benzo\[[^\]]+\]benz(?:imidazole|triazole|oxazole|thiazole)',
    re.IGNORECASE,
)
_OUR_SYSTEMATIC_FUSED = re.compile(
    r'naphtho\[|[a-z]+o\[\d,\d-[a-z]\](?:imidazole|triazole|oxazole|thiazole)',
    re.IGNORECASE,
)


def _norm(s: str) -> str:
    return " ".join(s.strip().lower().split())


def _verdict(
    verdict: str,
    note: Optional[str],
    iupac_source: Optional[str],
    iupac_pin: Optional[str] = None,
) -> dict:
    return {
        "iupac_verdict": verdict,
        "iupac_note": note,
        "iupac_source": iupac_source,
        "iupac_pin": iupac_pin,
    }


def classify_mismatch(our_name: str, pubchem_name: str) -> dict:
    """IUPAC 2013 規則に基づき、命名差異を分類する。

    Args:
        our_name: 本ライブラリが出力した名前 (テストの expected 値)。
        pubchem_name: PubChem から取得した IUPACName。

    Returns:
        dict with keys: iupac_verdict, iupac_note, iupac_source, iupac_pin
    """
    our = our_name.strip()
    pub = pubchem_name.strip()

    # ------------------------------------------------------------------ #
    # pubchem_wrong: PubChem が IUPAC 2013 に従っていない                  #
    # ------------------------------------------------------------------ #

    # 1. 非標準ハイフン: PubChem が括弧閉じと次の置換基の間にハイフンを挿入
    #    例: (4-chlorophenyl)-phenylmethanone → (4-chlorophenyl)(phenyl)methanone
    #    パターン A: ハイフン除去で一致
    if re.sub(r'\)-', ')', pub) == our:
        return _verdict(
            "pubchem_wrong",
            "PubChem が置換基グループ間に非標準のハイフンを挿入している",
            "IUPAC 2013 P-16.3.3",
        )
    # パターン B: PubChem がハイフンを挿入 + 2番目の括弧なし置換基を括弧付きにすると一致
    #    例: (4-chlorophenyl)-phenylmethanone → (4-chlorophenyl)(phenyl)methanone
    pub_no_hyph = re.sub(r'\)-', ')', pub)
    # )-X(... を )(X)(... に変換: カッコなし単語に括弧追加
    pub_add_paren = re.sub(r'\)([a-z][a-z0-9]*)', r')(\1)', pub_no_hyph)
    if _norm(pub_add_paren) == _norm(our):
        return _verdict(
            "pubchem_wrong",
            "PubChem が非標準ハイフンを挿入し、置換基括弧を省略している (P-16.3.3, P-16.3.4.1)",
            "IUPAC 2013 P-16.3.3, P-16.3.4.1",
        )

    # 2. 括弧の省略: 局所番号または置換を含む置換基前置語に括弧が必要
    #    IUPAC 2013 P-16.3.4.1: 局所番号を含む / 置換されている前置語は括弧必須
    if '(' in our and '(' not in pub:
        contents = re.findall(r'\(([^)]+)\)', our)
        for content in contents:
            has_locant = bool(re.search(r'\d', content))
            # 置換されているグループ: ハロ、スルホニル等で始まる
            is_substituted = bool(re.search(
                r'^(?:trifluoro|difluoro|chloro|bromo|fluoro|iodo'
                r'|methanesulfonyl|methylsulfonyl|ethane|methane)',
                content, re.IGNORECASE,
            ))
            if has_locant or is_substituted:
                # PubChem 名が括弧なしの同内容を含む場合
                if _norm(content) in _norm(pub):
                    return _verdict(
                        "pubchem_wrong",
                        f"PubChem が '{content}' の括弧を省略。"
                        "局所番号または置換を含む前置語には括弧が必須 (IUPAC 2013 P-16.3.4.1)",
                        "IUPAC 2013 P-16.3.4.1",
                    )

    # 3. 括弧除去 + α の変換で PubChem 名と一致するか確認 (組み合わせ違反)
    if '(' in our:
        our_no_paren = re.sub(r'\(([^)]+)\)', r'\1', our)
        # 3a. 括弧除去のみ
        if _norm(our_no_paren) == _norm(pub):
            return _verdict(
                "pubchem_wrong",
                "PubChem が IUPAC 2013 P-16.3.4.1 で必須の括弧を省略している",
                "IUPAC 2013 P-16.3.4.1",
            )
        # 3b. 括弧除去 + ケトン局所番号 (-1-one → -one) 除去
        our_np_nk = re.sub(r'an-1-one\b', 'anone', our_no_paren)
        if _norm(our_np_nk) == _norm(pub):
            return _verdict(
                "pubchem_wrong",
                "PubChem が括弧 (P-16.3.4.1) とケトンの局所番号 (P-65.1.2.3) を省略している",
                "IUPAC 2013 P-16.3.4.1, P-65.1.2.3",
            )
        # 3c. 括弧除去 + N-yl 局所番号除去: N-(pyridin-4-yl) → N-pyridin-4-ylの後の -yl 統合
        our_np_na = re.sub(r'an-1-one\b', 'anone',
                           re.sub(r'\)([a-z])', r'\1', our_no_paren))
        if _norm(our_np_na) == _norm(pub):
            return _verdict(
                "pubchem_wrong",
                "PubChem が括弧 (P-16.3.4.1) を省略している",
                "IUPAC 2013 P-16.3.4.1",
            )

    # 1c. ハイフン後の非括弧置換基に括弧追加: (X)-phenylmethanone → (X)(phenyl)methanone
    #     PubChem が括弧外置換基の前にハイフンを挿入し、かつ括弧を省略するパターン
    if ')-' in pub and '(' in pub:
        # )-word(methanone|etc) を )(word)(methanone|etc) に変換して比較
        pub_fixed = re.sub(
            r'\)-([a-z]+)(methanone|amine|amide|diazene|sulfide|disulfide|oxide)',
            r')(\1)\2',
            pub,
        )
        if _norm(pub_fixed) == _norm(our):
            return _verdict(
                "pubchem_wrong",
                "PubChem が非標準ハイフンを挿入し置換基括弧を省略 (P-16.3.3, P-16.3.4.1)",
                "IUPAC 2013 P-16.3.3, P-16.3.4.1",
            )

    # 4. アルケニル末端局所番号の省略: but-2-en-1-yl → but-2-enyl
    if re.search(r'\w+en-\d+-yl\b', our) and not re.search(r'\w+en-\d+-yl\b', pub):
        our_norm = re.sub(r'(en)-\d+-(yl)', r'\1\2', our)
        if _norm(our_norm) == _norm(pub):
            return _verdict(
                "pubchem_wrong",
                "PubChem がアルケニル基の末端局所番号を省略。IUPAC 2013 では必須",
                "IUPAC 2013 P-31.1.2.2",
            )

    # 5. E/Z 立体記述子の局所番号省略: (2E) → (E) — 単独または複合
    if re.search(r'\(\d+[EZ]\)', our) and re.search(r'\([EZ]\)', pub):
        our_no_ez_loc = re.sub(r'\((\d+)([EZ])\)', r'(\2)', our)
        if _norm(our_no_ez_loc) == _norm(pub):
            return _verdict(
                "pubchem_wrong",
                "PubChem が E/Z 立体記述子の局所番号を省略。IUPAC 2013 P-93.4 では必須",
                "IUPAC 2013 P-93.4",
            )
        # 複合: E/Z 局所番号 + アルケニル末端局所番号の両方を省略
        our_no_ez_no_en = re.sub(r'(en)-\d+-(yl)', r'\1\2', our_no_ez_loc)
        if _norm(our_no_ez_no_en) == _norm(pub):
            return _verdict(
                "pubchem_wrong",
                "PubChem が E/Z 局所番号 (P-93.4) とアルケニル末端局所番号 (P-31.1.2.2) を省略",
                "IUPAC 2013 P-93.4, P-31.1.2.2",
            )

    # 6. ケトン名の局所番号省略: ethan-1-one → ethanone
    if re.search(r'an-1-one\b', our) and not re.search(r'an-1-one\b', pub):
        our_norm = re.sub(r'an-1-one\b', 'anone', our)
        if _norm(our_norm) == _norm(pub):
            return _verdict(
                "pubchem_wrong",
                "PubChem がケトン名の局所番号-1を省略。IUPAC 2013 P-65.1.2.3 では必須",
                "IUPAC 2013 P-65.1.2.3",
            )

    # ------------------------------------------------------------------ #
    # our_retained: 我々の名前が IUPAC 2013 保留名                         #
    # ------------------------------------------------------------------ #

    # 7. 窒素複素芳香環のエノール(-ol) vs ケト(-one) 互変異性体
    #    IUPAC 2013 P-31.1.3: 優勢互変異性体 (通常ケト形) が PIN
    if _HETERO_OL_PARENTS.search(our) and re.search(
        r'\b(?:1H|2H|3H)-.*-one\b|(?:1H|2H|3H)-\w+-\d+-one\b', pub, re.IGNORECASE
    ):
        return _verdict(
            "our_retained",
            f"'{our}' はエノール互変異性体形の名前 (保留形)。"
            f"IUPAC 2013 PIN はより安定なケト形 '{pub}'",
            "IUPAC 2013 P-31.1.3",
            iupac_pin=pub,
        )

    # 7b. より広いエノール/ケト: -ol と -one の対応
    if re.search(r'\w+-\d+(?:\(\d+H\))?-ol\b', our) and re.search(r'-one\b', pub):
        # ヘテロ原子を含む化合物のみ (ベンゼン環の -ol は phenol として別)
        if re.search(r'[nso]', our, re.IGNORECASE):
            return _verdict(
                "our_retained",
                f"'{our}' はエノール形の名前 (保留形)。"
                f"IUPAC 2013 PIN はケト形 '{pub}'",
                "IUPAC 2013 P-31.1.3",
                iupac_pin=pub,
            )

    # 8. チオール/チオン互変異性体 (硫黄類似体)
    if (re.search(r'\w+-\d+(?:\(\d+H\))?-thiol\b', our)
            and re.search(r'-thione\b', pub)):
        return _verdict(
            "our_retained",
            f"'{our}' はチオール互変異性体形 (保留形)。"
            f"IUPAC 2013 PIN はチオン形 '{pub}'",
            "IUPAC 2013 P-31.1.3",
            iupac_pin=pub,
        )

    # 9. 既知の保留慣用名 (完全一致)
    our_lower = our.lower()
    if our_lower in _RETAINED_EXACT:
        return _verdict(
            "our_retained",
            f"'{our}' は IUPAC 2013 保留名。PIN は '{pub}'",
            "IUPAC 2013 P-31.1.3",
            iupac_pin=pub,
        )

    # 10. 保留親骨格を含む化合物 (語末/語境界一致)
    #     4-methylbiphenyl のように語中に埋め込まれる場合も考慮 (word boundary 不使用)
    for parent in _RETAINED_PARENTS:
        # 末尾一致 OR 非アルファ文字で区切られた出現
        if our_lower.endswith(parent) or re.search(
            r'(?<![a-z])' + re.escape(parent) + r'(?![a-z])', our_lower
        ):
            return _verdict(
                "our_retained",
                f"我々の名前が保留親骨格 '{parent}' を使用。"
                f"PubChem は系統的 PIN '{pub}' を使用",
                "IUPAC 2013 P-31.1.3",
                iupac_pin=pub,
            )

    # 11. アミノ酸 L-/D- 表記 (P-102.4)
    if re.match(r'^[LD]-', our):
        return _verdict(
            "our_retained",
            "L-/D- アミノ酸名は IUPAC 2013 保留表記。"
            f"PubChem は (R/S) を用いた系統的 PIN '{pub}' を使用",
            "IUPAC 2013 P-102.4",
            iupac_pin=pub,
        )

    # ------------------------------------------------------------------ #
    # pubchem_retained: PubChem が保留縮合環成分を使用                     #
    # ------------------------------------------------------------------ #

    # 12. PubChem が benzo[f]benzimidazole 等の保留縮合環を使用
    if _PUBCHEM_RETAINED_FUSED.search(pub) and _OUR_SYSTEMATIC_FUSED.search(our):
        return _verdict(
            "pubchem_retained",
            f"PubChem が保留縮合環成分を含む '{pub}' を使用。"
            f"我々の系統的融合命名 '{our}' がより IUPAC 2013 に適合",
            "IUPAC 2013 P-31.1.3",
            iupac_pin=our,
        )

    # ------------------------------------------------------------------ #
    # tautomer: 指示水素位置の差 (互変異性体)                               #
    # ------------------------------------------------------------------ #

    # 13. 指示水素 (1H-, 2H-, 3H-, (3H)- 等) を除くと一致
    def _strip_ih(s: str) -> str:
        s = re.sub(r'\b\d+H-', '', s)
        s = re.sub(r'\(\d+H\)-?', '', s)
        return _norm(s)

    if _strip_ih(our) == _strip_ih(pub):
        return _verdict(
            "tautomer",
            "名前の差異は指示水素の位置のみ — 互変異性体の命名",
            "IUPAC 2013 P-31.1.3",
        )

    # 14. 指示水素記法の差: 2(3H)-one vs 3H-...-2-one
    our_ih = re.sub(r'-(\d+)\((\d+)H\)-(\w+)', r'-\3(\2H)-\1', our)
    if _norm(our_ih) == _norm(pub) or _norm(our) == _norm(
        re.sub(r'-(\d+)\((\d+)H\)-(\w+)', r'-\3(\2H)-\1', pub)
    ):
        return _verdict(
            "tautomer",
            "指示水素の記法が異なるが同一互変異性体を表す",
            "IUPAC 2013 P-31.1.3",
        )

    # ------------------------------------------------------------------ #
    # needs_review: 自動分類不能                                            #
    # ------------------------------------------------------------------ #

    return _verdict("needs_review", None, None)
