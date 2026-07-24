"""Dump Python smiles2iupac names for the corpus (IUPAC naming gate for molrs)."""
import gzip, json, os, sys
from pathlib import Path
REPO=Path(__file__).resolve().parent.parent
MOLRS=Path(os.environ.get("MOLRS_ROOT", Path.home()/"ghq/github.com/nobyt/molrs"))
sys.path.insert(0, str(REPO/"src"))
from smiles2iupac import smiles_to_iupac
rows=[json.loads(l) for l in (MOLRS/"corpus/corpus.jsonl").open()]
out=MOLRS/"corpus/names.jsonl.gz"
n_ok=0
with gzip.open(out,"wt") as f:
    for r in rows:
        smi=r["smiles"]
        try:
            name=smiles_to_iupac(smi)
            n_ok+=1
        except Exception as e:
            name=""
        f.write(json.dumps({"s":smi,"name":name},separators=(",",":"))+"\n")
print(f"named {n_ok}/{len(rows)} -> {out}")
