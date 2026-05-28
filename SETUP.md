# smiles2iupac セットアップガイド

## 必要環境
- Python 3.11+
- uv
- RDKit (conda 経由が最も安定)

## インストール手順

### 1. リポジトリのセットアップ

```bash
cd smiles2iupac

# uv で仮想環境を作成 + 依存関係インストール
uv sync
```

### 2. RDKit のインストール (conda 推奨)

RDKit は conda-forge から取得するのが最も安定しています。

```bash
# conda 環境を使う場合
conda create -n smiles2iupac python=3.11
conda activate smiles2iupac
conda install -c conda-forge rdkit

# pip のみの場合 (Python 3.11 以下推奨)
pip install rdkit
```

### 3. 開発用インストール

```bash
# プロジェクトルートから
uv pip install -e ".[dev]"
# または
pip install -e ".[dev]"
```

## テストの実行

### RDKit 不要のテスト (name_assembler のユニットテスト)

```bash
cd smiles2iupac
PYTHONPATH=src python -m pytest tests/test_name_assembler.py -v
```

### 全テスト (RDKit 必要)

```bash
uv run pytest tests/ -v
# または
python -m pytest tests/ -v
```

## 使用例

```python
from smiles2iupac import smiles_to_iupac

# アルカン
print(smiles_to_iupac("CCCCCC"))        # hexane
print(smiles_to_iupac("CC(C)C"))        # 2-methylpropane
print(smiles_to_iupac("CCC(C)CC"))      # 3-methylpentane

# アルコール
print(smiles_to_iupac("CCO"))           # ethan-1-ol
print(smiles_to_iupac("CC(O)CCC"))      # pentan-2-ol

# ハロゲン化物
print(smiles_to_iupac("CCCCl"))         # 1-chloropropane
print(smiles_to_iupac("CC(Cl)C"))       # 2-chloropropane

# カルボン酸
print(smiles_to_iupac("CC(=O)O"))       # ethanoic acid
print(smiles_to_iupac("CCC(=O)O"))      # propanoic acid

# アルデヒド
print(smiles_to_iupac("CC=O"))          # ethanal

# ケトン
print(smiles_to_iupac("CC(=O)C"))       # propan-2-one

# アルケン
print(smiles_to_iupac("CC=CC"))         # but-2-ene

# アルキン
print(smiles_to_iupac("CC#CC"))         # but-2-yne

# 三環縮合多環系
print(smiles_to_iupac("c1ccc2cc3ccccc3cc2c1"))    # anthracene
print(smiles_to_iupac("c1ccc2ccc3ccccc3c2c1"))    # phenanthrene
print(smiles_to_iupac("c1ccc2c(C)c3ccccc3cc2c1")) # 9-methylanthracene
```

## 対応状況

| Phase | 対象 | 状態 |
|-------|------|------|
| Phase 1 | 直鎖・分岐アルカン | ✅ |
| Phase 1 | アルコール | ✅ |
| Phase 1 | ハロゲン化物 | ✅ |
| Phase 1 | アルケン・アルキン | ✅ |
| Phase 2 | カルボン酸 | ✅ |
| Phase 2 | アルデヒド | ✅ |
| Phase 2 | ケトン | ✅ |
| Phase 2 | 立体化学 (R/S, E/Z) | ✅ |
| Phase 3 | シクロアルカン | ✅ |
| Phase 3 | 芳香族 (ベンゼン誘導体) | ✅ |
| Phase 4 | ベンゼン官能基誘導体 (phenol, benzoic acid 等) | ✅ |
| Phase 5 | ナフタレン / 置換ナフタレン | ✅ |
| Phase 6 | 三環縮合多環系 (anthracene, phenanthrene) | ✅ |
| Phase 7 | 含窒素官能基 (amine, nitrile, nitro prefix) | ✅ |
| 改善 | 環外アルデヒド → carbaldehyde 形式 | ✅ |
| 改善 | 複合置換基の bis/tris 接頭辞 | ✅ |

## 注意事項

- 複合置換基 (isopropyl 等) は `1-methylethyl` 形式で命名されます
- 立体化学の R/S・E/Z は RDKit の CIPCode を利用して取得しています
