// SMILES 목록을 JS 타이퍼로 라벨링해 JSON 으로 출력 (Python 대조용)
const fs = require("fs");
const path = require("path");
const T = require("../src/typer.js");
const initRDKit = require("@rdkit/rdkit");

(async () => {
  const R = await initRDKit();
  const ff = new T.DreidingFF(JSON.parse(fs.readFileSync(path.join(__dirname, "../build/tables.json"))));
  const out = {};
  for (const smi of fs.readFileSync(process.argv[2], "utf8").split(/\s+/).filter(Boolean)) {
    const m = R.get_mol(smi);
    if (!m || !m.is_valid()) { out[smi] = { parse_error: true }; continue; }
    m.add_hs_in_place();
    const g = T.graphFromRDKit(m);
    m.delete();
    try {
      const r = T.assignLabels(g, ff);
      out[smi] = { labels: r.labels, warnings: r.warnings };
    } catch (e) {
      if (!(e instanceof T.TypingError)) throw e;
      out[smi] = { error: e.message };
    }
  }
  process.stdout.write(JSON.stringify(out, null, 0));
})();
