# SMILES 목록을 Python 타이퍼로 라벨링해 JSON 으로 출력 (JS 대조용)
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from rdkit import Chem  # noqa: E402
from dreiding_typer import DreidingFF, TypingError, assign_labels  # noqa: E402

ff = DreidingFF()
out = {}
for smi in open(sys.argv[1]).read().split():
    m = Chem.MolFromSmiles(smi)
    if m is None:
        out[smi] = {"parse_error": True}
        continue
    m = Chem.AddHs(m)
    try:
        labels, warnings = assign_labels(m, ff)
        out[smi] = {"labels": labels, "warnings": warnings}
    except TypingError as e:
        out[smi] = {"error": str(e)}
json.dump(out, sys.stdout, ensure_ascii=False, indent=0)
