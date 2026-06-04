# smiles2iupac

SMILES 文字列を **IUPAC 2013 推奨名** に変換するライブラリ＆CLI ツールです。

[English README](README.md)

SMILES のパース・分子グラフ構築に [RDKit](https://www.rdkit.org/) を使用し、
命名ロジックはスクラッチ実装（3,890 超のテストで検証済み）。

---

## 対応範囲

| カテゴリ | 例 |
|---|---|
| アルカン・アルケン・アルキン | butane, but-2-ene, but-2-yne |
| ハロゲン化物・アルコール | 2-chloropropane, propan-1-ol |
| カルボン酸・エステル・アミド | acetic acid, ethyl acetate, acetamide |
| 環状化合物 | cyclohexane, benzene, naphthalene |
| 立体化学 | (R)-alanine, (E)-but-2-ene |
| 複素環 | pyridine, furan, imidazole |
| 含窒素官能基 | amine, nitrile, amidine, hydrazide |
| 含硫黄官能基 | thiol, sulfonic acid, sulfoxide |
| 含リン・ケイ素化合物 | phosphate ester, trimethylsilanol |
| その他 | エーテル、過酸化物、イソシアナートなど |

---

## インストール

### 必要環境

- Python 3.11 以上
- [uv](https://docs.astral.sh/uv/) または pip
- RDKit 2023.9 以上

### uv（推奨）

```bash
# リポジトリをクローン
git clone https://github.com/yourname/smiles2iupac.git
cd smiles2iupac

# 仮想環境の作成と依存関係のインストール
uv sync

# CLI をパスに追加して使用する場合（オプション）
uv tool install .
```

`uv sync` 後は `.venv/bin/smiles2iupac` で CLI が使えます。

### pip

```bash
git clone https://github.com/yourname/smiles2iupac.git
cd smiles2iupac

pip install .
```

> **RDKit に関する注意**: pip でインストールされる `rdkit` パッケージはほとんどの環境で動作しますが、
> 問題が発生する場合は conda 経由のインストールが安定しています。
>
> ```bash
> conda install -c conda-forge rdkit
> ```

---

## CLI の使い方

### 基本

```bash
smiles2iupac "CC(=O)O"
# acetic acid

smiles2iupac "c1ccccc1"
# benzene

smiles2iupac "C[C@@H](N)C(=O)O"
# D-alanine
```

### 複数入力（stdin）

1 行に 1 つの SMILES を渡すと、対応する IUPAC 名を 1 行ずつ出力します。

```bash
printf 'CC\nCCC\nCCCC\n' | smiles2iupac -
# ethane
# propane
# butane

smiles2iupac - < smiles.txt
```

### ヘルプ

```bash
smiles2iupac --help
smiles2iupac --version
```

### 使用例

```bash
smiles2iupac "CC(C)C"
# 2-methylpropane

smiles2iupac "c1ccc(O)cc1"
# phenol

smiles2iupac "CC(=O)OCC"
# ethyl acetate

smiles2iupac "c1ccc(cc1)C(=O)O"
# benzoic acid

smiles2iupac "ClC(Cl)Cl"
# trichloromethane
```

---

## ライブラリとしての使い方

```python
from smiles2iupac import smiles_to_iupac

print(smiles_to_iupac("CC"))           # ethane
print(smiles_to_iupac("c1ccccc1"))     # benzene
print(smiles_to_iupac("CC(=O)O"))      # acetic acid
print(smiles_to_iupac("C[C@@H](N)C(=O)O"))  # D-alanine
```

---

## 開発

```bash
# テストの実行
uv run pytest

# 特定フェーズのテストのみ実行
uv run pytest tests/test_phase1.py -v
```

---

## ライセンス

MIT
