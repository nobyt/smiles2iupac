# smiles2iupac リファクタリング計画

Claude Sonnet が 1 ステップずつ実行できるように書かれた改善計画。
各ステップは独立してコミット可能で、上から順に実施することを想定している
(Step 1 のテスト高速化が以降の全ステップの回帰確認を支えるため、順序を守ること)。

## 現状の定量データ(2026-07 時点)

| 項目 | 値 | 問題 |
|---|---|---|
| `heterocycle_handler.py` | 5,628 行 | 最大のホットスポット |
| `_apply_hetero_suffixes()` | **1,900 行**(1関数) | ring 名の文字列比較 ~896 箇所、ハードコード名 112 個 |
| `functional_group.py detect_groups()` | **1,883 行**(1関数) | 全官能基検出が 1 関数に直列 |
| `group_namers.py` | 7,264 行 / 159 関数 | 分割は済んでいるがモジュールが巨大 |
| `__init__.py _smiles_to_iupac_raw()` | 235 行、handler 直列呼び出し 28 個 | **呼び出し順序がバグ源**(Phase 853 で実証) |
| `_FUSED_LOCANT_MAP` | 547 キー、**None ロカント 3,342 個** | None は置換基を**無言で捨てる**(Phase 849/851 のバグ原因) |
| Phase 番号サフィックス付き識別子 | 43 個(`_oxo_base_849` 等) | 意味が読み取れない命名 |
| テストファイル | `test_phase*.py` が 842 個 | フルスイート **35 分**(12,531 件) |

## 安全ルール(全ステップ共通)

1. **テストが唯一の仕様**。`tests/` の期待値は絶対に変更しない。
   期待値を変えたくなったらそのステップを中止してユーザーに報告する。
2. 1 ステップ = 1 コミット。コミット前に必ず回帰テストを実行する。
3. テスト実行は `uv run pytest`(システム Python は 3.8 で動かない)。
4. 動作を変えないリファクタリングでは、変更前後で対象 SMILES 群の出力が
   一致することをスクリプトで確認する(各ステップに検証コマンドを記載)。
5. `N-methylpyrrolidine` 等の "N" ロカントは IUPAC 2013 P-31.1.3.4 準拠で
   **正しい**。`1-methyl...` に「修正」しないこと(過去に誤判断しかけた実績あり)。

---

## Step 1: テスト並列化(pytest-xdist 導入)

**目的**: フルスイート 35 分 → 数分に短縮し、以降のリファクタリングの
回帰確認コストを下げる。

**手順**:
1. `uv add --group dev pytest-xdist`
2. `pyproject.toml` の `[tool.pytest.ini_options]`(なければ追加)に
   `addopts = "-n auto"` を設定する。
   個別デバッグ時は `uv run pytest tests/test_phaseNNN.py -n 0` で無効化できる。
3. `uv run pytest tests/ -q` を実行し、**12,531 passed / 2 skipped** が
   維持されることと所要時間を確認する。

**検証**: パス数が変わらないこと。並列化で落ちるテストが出た場合は
テスト間の状態共有(キャッシュファイル等)が原因なので、`-n auto` を外して
そのテストを特定し、ユーザーに報告する。

**コミット例**: `Refactor step 1: parallelize test suite with pytest-xdist`

---

## Step 2: `_FUSED_LOCANT_MAP` の None ロカント監査

**目的**: None ロカントは「その原子の置換基を無言で名前から落とす」ことを
意味し、Phase 849(アクリジン N-10)・Phase 851(ベンゾチアゾール N-3)の
バグ原因だった。3,342 個の None のうち「ヘテロ原子(N/O/S)なのに None」の
ものを洗い出し、既知の危険箇所を可視化する。

**手順**:
1. `tools/audit_locant_map.py` を新規作成する。内容:
   - `_FUSED_LOCANT_MAP` の各キー(canonical SMILES)を RDKit でパースし、
     None が割り当てられた原子のうち元素が N/O/S のものを列挙する。
   - 出力形式: `SMILES  atom_idx  element  → None`
2. 実行して結果を `docs/locant_map_audit.md` に保存する。
3. **このステップではマップを修正しない**(修正は個別 Phase として、
   テストとセットで行う)。監査ツールとレポートのみコミットする。

**検証**: `uv run python tools/audit_locant_map.py` がエラーなく完走すること。

**コミット例**: `Refactor step 2: add locant-map audit tool (flags silent substituent drops)`

---

## Step 3: Phase 849–853 の oxo/thione ブロックをデータテーブル化

**目的**: `heterocycle_handler.py` の `_apply_hetero_suffixes` 末尾にある
Phase 849–851 の oxo ブロック(~L2632 以降)は
`if full_base == "quinoline": ... elif full_base == "quinazoline": ...` の
長い分岐で、リング追加のたびに elif が増える。宣言的テーブルに置き換える。

**手順**:
1. モジュールレベルに変換テーブルを定義する:
   ```python
   # (ring name, oxo locant) → N locant that suppresses indicated-H
   _KETO_SUFFIX_TABLE: dict[tuple[str, int], int] = {
       ("quinoline", 2): 1,
       ("quinoline", 4): 1,
       ("isoquinoline", 1): 2,
       ("isoquinoline", 3): 2,
       ("acridine", 9): 10,
       ("phenanthridine", 6): 5,
       ("phthalazine", 1): 2,
       ("quinazoline", 4): 3,
       ("quinazoline", 2): 1,
       ("cinnoline", 3): 2,
       ("cinnoline", 4): 1,
       ("quinoxaline", 2): 1,
       ("1H-benzimidazole", 2): 3,
       ("1,3-benzoxazole", 2): 3,
       ("1,3-benzothiazole", 2): 3,
   }
   ```
2. Phase 849 ブロックの if/elif チェーンを、このテーブル参照 1 本にする:
   ```python
   _n_loc = _KETO_SUFFIX_TABLE.get((full_base, _oxo_loc))
   if _n_loc is not None:
       if _n_sub(_n_loc):
           _oxo_base = f"{_stem(full_base)}-{_oxo_loc}-one"
       else:
           _oxo_base = f"{_stem(full_base)}-{_oxo_loc}({_n_loc}H)-one"
   ```
   `_stem()` は語尾の e を落とすヘルパ(`quinoline` → `quinolin`、
   `1,3-benzoxazole` → `1,3-benzoxazol`)。既存の出力文字列
   (例: `quinazolin-4(3H)-one`)と**完全一致**することを必ず確認する。
3. 検証スクリプト: Phase 849/850/851/853 の全テスト SMILES について
   変更前後の出力を比較する。
   `uv run pytest tests/test_phase849.py tests/test_phase850.py tests/test_phase851.py tests/test_phase853.py -n 0 -q`
4. フルスイート実行 → コミット。

**コミット例**: `Refactor step 3: table-drive keto suffix rules (phases 849-851)`

---

## Step 4: 汎用 `(nH)` 除去を正規表現化

**目的**: Phase 852 の indicated-H 除去は `range(20)` ループで
`f"({i}H)"` を試している(`heterocycle_handler.py` ~L2620)。
正規表現 1 回に置き換え、`(1H,3H)` 型を誤って触らないことをテストで固定する。

**手順**:
1. 置き換え:
   ```python
   _IH_RE = re.compile(r"\((\d+)H\)")
   m = _IH_RE.search(base_with_suffix)
   if m and any(l == int(m.group(1)) for l, _ in other):
       base_with_suffix = base_with_suffix.replace(m.group(0), "")
   ```
   ※ `(1H,3H)` は `\(\d+H\)` にマッチしないので現行と同じく除去されない。
2. `tests/test_phase852.py` を `-n 0` で実行 → フルスイート → コミット。

**コミット例**: `Refactor step 4: regex-based indicated-H stripping`

---

## Step 5: `_apply_hetero_suffixes`(1,900 行)の分割

**目的**: 1 関数 1,900 行・ring 名比較 896 箇所は保守限界。
外部から見た入出力(`full_base`, `substituents` → 名前文字列)を変えずに、
内部を 3 つのサブ関数に分割する。

**手順**(サブステップごとにテスト・コミット):
1. **5a**: hydroxy ブロック(~L838–1600)を
   `_hydroxy_to_one_suffix(full_base, substituents) -> str | None` に抽出。
   元の関数は戻り値が None でなければそれを返すだけにする。
2. **5b**: sulfanyl ブロック(~L1737–2630)を
   `_sulfanyl_to_thione_suffix(...)` に同様に抽出。
3. **5c**: Step 3 で作った oxo テーブル処理を `_keto_to_one_suffix(...)` に抽出。
4. 各サブステップで: 該当 Phase テスト(ファイル名は docstring の
   Phase 番号で grep できる)→ フルスイート → コミット。
5. **注意**: 変数のスコープ移動のみ行い、ロジック・文字列は 1 文字も
   変えない。diff は「移動 + def 行 + return 行」だけになるのが理想。

**コミット例**: `Refactor step 5a: extract _hydroxy_to_one_suffix from _apply_hetero_suffixes`

---

## Step 6: `__init__.py` ディスパッチチェーンのレジストリ化

**目的**: `_smiles_to_iupac_raw` は 28 個のハンドラを直列に呼ぶ構造で、
**呼び出し順序自体が仕様**になっている(Phase 853 は「PGRP_DISPATCH が
fused-hetero チェックより先に走る」順序バグだった)。順序を明示した
レジストリに変え、順序変更を diff で追えるようにする。

**手順**:
1. 各ハンドラを `(graph, get_atom) -> str | None` のシグネチャに揃える
   薄いラッパを作り、モジュールレベルのリストに順序どおり登録する:
   ```python
   _EARLY_HANDLERS: list[Callable[..., str | None]] = [
       _name_thiourea_if_match,
       _carbonohydrazide_handler,
       _fused_hetero_before_pgrp,   # Phase 853
       _pgrp_dispatch_handler,
       ...
   ]
   ```
2. `_smiles_to_iupac_raw` 本体は
   `for h in _EARLY_HANDLERS: r = h(graph, get_atom); if r: return r`
   のループにする。
3. **重要**: `_pgrp` を必要とするハンドラがあるため、`_pgrp` は
   ループ前に 1 回計算してハンドラに渡す(クロージャまたは引数追加)。
   遅延インポートしている関数はラッパ内で import する。
4. 1 ハンドラずつ移す必要はない。ただし移行後にフルスイートを必ず実行。

**コミット例**: `Refactor step 6: ordered handler registry in _smiles_to_iupac_raw`

---

## Step 7: `detect_groups`(1,883 行)の分割

**目的**: 全官能基検出が 1 関数。Step 5 と同じ「抽出のみ」方式で、
検出単位ごとの関数(`_detect_esters`, `_detect_thio_acids`, ...)に分割する。

**手順**:
1. 関数内のコメント区切り(`# C(=O)-SH: thioic S-acid` 等)を境界として
   10〜20 個のブロックに分ける。
2. 各ブロックを `_detect_<name>(graph, get_atom, groups: list) -> None`
   (groups に append する形)として抽出し、`detect_groups` は
   それらを順に呼ぶだけにする。**検出順序を維持すること**。
3. 3〜5 ブロック抽出するごとにフルスイート → コミット。

**コミット例**: `Refactor step 7a: extract acid/ester detectors from detect_groups`

---

## Step 8: Phase サフィックス変数のリネーム

**目的**: `_oxo_base_849`、`_ih_key_852` など 43 個の Phase 番号付き識別子を
意味のある名前(`_keto_suffix_base`、`_indicated_h_key`)に改名する。
Phase 由来の説明はコメントに残す。

**手順**:
1. `grep -rnE '_[a-z_]+_[0-9]{3}\b' src/` で一覧を出す。
2. ファイルごとに機械的にリネーム(スコープはすべて関数ローカルなので
   衝突リスクは低いが、同一関数内での重複名に注意)。
3. フルスイート → コミット。Step 3/5 実施後なら対象は大きく減っているはず。

**コミット例**: `Refactor step 8: rename phase-suffixed identifiers to semantic names`

---

## Step 9: テストの整理(任意・低優先)

**目的**: `test_phase*.py` が 842 ファイル。実行には支障ないが、
「この化合物クラスのテストはどこか」が探しにくい。

**方針**(実施する場合):
- 新規テストは今後 `tests/domains/test_<domain>.py`
  (例: `test_fused_thiones.py`)に追加し、phase ファイルは凍結する。
- 既存ファイルの統合・移動は**行わない**(git 履歴とテスト名の
  トレーサビリティが失われるため)。README にこの方針を明記する。

**コミット例**: `Refactor step 9: document test organization policy`

---

## 実施しないこと(明示)

- テスト期待値の変更(仕様変更にあたる)
- `_FUSED_LOCANT_MAP` の None の一括修正(1 件ずつ Phase として、
  PubChem 照合とテストを添えて行う)
- `group_namers.py` のファイル分割(159 関数あるが相互依存が薄く、
  現状の実害が小さい。Step 5–7 完了後に再評価)

## 進め方の目安

| ステップ | 期待効果 | リスク | 優先度 |
|---|---|---|---|
| 1. pytest-xdist | 開発速度が数倍 | 低 | ★★★ |
| 2. locant 監査 | 将来バグの可視化 | なし(読み取りのみ) | ★★★ |
| 3. keto テーブル化 | Phase 追加が 1 行に | 低 | ★★★ |
| 4. regex 化 | 微小 | 低 | ★★ |
| 5. suffix 関数分割 | 最大ファイルの保守性 | 中(移動ミス) | ★★ |
| 6. handler レジストリ | 順序バグの再発防止 | 中 | ★★ |
| 7. detect_groups 分割 | 同上 | 中 | ★ |
| 8. リネーム | 可読性 | 低 | ★ |
| 9. テスト方針 | ドキュメントのみ | なし | ★ |
