"""molrs InChI: コーパス全 SMILES の RDKit InChI / InChIKey / AuxInfo をダンプする。

molrs 側の inchi モジュールの差分検証 (式層・正準番号 /N:・フル文字列・キー) に
使う。AuxInfo の /N: は公式の正準番号付け (最難関 I3 の単独ゲート)。

出力: MOLRS_ROOT/corpus/inchi_dump.jsonl.gz (1 行 = 1 分子)
  {"s": smiles,
   "inchi": "InChI=1S/...",       # 標準 InChI (取得失敗時は "")
   "key": "XXXX-...",             # InChIKey
   "formula": "C2H4O2",           # InChI 式層 (先頭セグメント)
   "n": [[comp0 の元原子 idx (1 始まり canonical 順)], ...]}  # AuxInfo /N:
"""
from __future__ import annotations

import gzip
import json
import os
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import inchi
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

REPO_ROOT = Path(__file__).resolve().parent.parent
MOLRS_ROOT = Path(os.environ.get("MOLRS_ROOT", Path.home() / "ghq/github.com/nobyt/molrs"))
CORPUS = MOLRS_ROOT / "corpus" / "corpus.jsonl"
OUT = MOLRS_ROOT / "corpus" / "inchi_dump.jsonl.gz"


def parse_aux_numbering(aux: str) -> list[list[int]]:
    """AuxInfo の /N:comp1;comp2;... を [[canonical順の元 idx(1始まり)], ...] に。"""
    if not aux:
        return []
    for seg in aux.split("/"):
        if seg.startswith("N:"):
            body = seg[2:]
            comps = []
            for comp in body.split(";"):
                if comp:
                    comps.append([int(x) for x in comp.split(",")])
            return comps
    return []


def dump_one(smiles: str) -> dict:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return {"s": smiles, "inchi": "", "key": "", "formula": "", "n": []}
    try:
        res = inchi.MolToInchiAndAuxInfo(m)
        ic = res[0] or ""
        aux = res[-1] or ""
    except Exception:  # noqa: BLE001
        ic, aux = "", ""
    key = inchi.InchiToInchiKey(ic) if ic else ""
    formula = ""
    if ic.startswith("InChI=1S/"):
        formula = ic[len("InChI=1S/"):].split("/", 1)[0]
    return {
        "s": smiles,
        "inchi": ic,
        "key": key or "",
        "formula": formula,
        "n": parse_aux_numbering(aux),
    }


def main() -> int:
    rows = [json.loads(l) for l in CORPUS.open()]
    n_ok = 0
    with gzip.open(OUT, "wt", encoding="utf-8") as f:
        for row in rows:
            rec = dump_one(row["smiles"])
            if rec["inchi"]:
                n_ok += 1
            f.write(json.dumps(rec, separators=(",", ":")) + "\n")
    print(f"dumped {len(rows)} molecules ({n_ok} with InChI) -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
