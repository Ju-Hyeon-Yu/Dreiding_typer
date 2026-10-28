// 화면 동작 — RDKit.js 준비, 입력 처리, 결과(경고·타입 카드·원자 표·뷰어·.lt) 표시
(function () {
  "use strict";
  const $ = id => document.getElementById(id);
  const T = window.DreidingTyper, P = window.DreidingPipeline;
  const tables = JSON.parse($("dreiding-tables").textContent);
  const ff = new T.DreidingFF(tables);
  const viewer = new window.MolViewer($("view"), $("tip"));
  let R = null, file = null, result = null;

  $("srcInfo").textContent = `moltemplate ${tables.moltemplate || ""}, ${Object.keys(tables.labels).length}개 라벨`.replace("  ", " ");

  // ------------------------------------------------------------ RDKit.js 준비 (압축 wasm 풀기)
  async function loadRDKit() {
    const b64 = $("rdkit-wasm").textContent.trim();
    const bin = Uint8Array.from(atob(b64), c => c.charCodeAt(0));
    if (typeof DecompressionStream === "undefined") throw new Error("이 브라우저는 압축 해제를 지원하지 않습니다. 최신 Edge 나 Chrome 을 쓰십시오");
    const stream = new Blob([bin]).stream().pipeThrough(new DecompressionStream("gzip"));
    const wasm = new Uint8Array(await new Response(stream).arrayBuffer());
    // 이 RDKit.js 빌드는 wasmBinary 옵션을 읽지 않으므로 instantiateWasm 으로 직접 넘긴다 (file:// 에서 fetch 불가)
    return window.initRDKitModule({
      instantiateWasm(imports, receive) {
        WebAssembly.instantiate(wasm, imports).then(r => receive(r.instance, r.module));
        return {};
      },
    });
  }
  loadRDKit().then(mod => {
    R = mod;
    $("status").textContent = `준비되었습니다 (RDKit ${R.version()}).`;
    $("status").classList.add("ready");
    $("run").disabled = false;
  }).catch(e => {
    $("status").textContent = "RDKit 엔진을 불러오지 못했습니다.";
    showMsgs([], `RDKit 엔진을 불러오지 못했습니다: ${e.message}`);
  });

  // ------------------------------------------------------------ 입력
  function setFile(f) {
    if (!f) { file = null; $("fileChip").hidden = true; $("drop").hidden = false; $("file").value = ""; return; }
    f.text().then(text => {
      file = { name: f.name, text };
      $("fileName").textContent = `${f.name} (${(f.size / 1024).toFixed(1)} KB)`;
      $("fileChip").hidden = false;
      $("drop").hidden = true;
    });
  }
  $("drop").addEventListener("click", () => $("file").click());
  $("drop").addEventListener("keydown", e => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); $("file").click(); } });
  $("file").addEventListener("change", e => setFile(e.target.files[0]));
  $("clearFile").addEventListener("click", () => setFile(null));
  ["dragenter", "dragover"].forEach(ev => document.addEventListener(ev, e => { e.preventDefault(); $("drop").classList.add("over"); }));
  ["dragleave", "drop"].forEach(ev => document.addEventListener(ev, e => { e.preventDefault(); if (ev === "drop" || e.target === document.documentElement) $("drop").classList.remove("over"); }));
  document.addEventListener("drop", e => { const f = e.dataTransfer && e.dataTransfer.files[0]; if (f) setFile(f); });
  $("smiles").addEventListener("keydown", e => { if (e.key === "Enter" && !$("run").disabled) run(); });
  document.querySelectorAll("[data-example]").forEach(b => b.addEventListener("click", () => {
    $("smiles").value = b.dataset.example;
    setFile(null);
    if (R) run();
  }));
  $("molName").addEventListener("input", () => {
    const v = $("molName").value.trim();
    if (/^[A-Za-z_][A-Za-z0-9_]*$/.test(v) && $("outName").dataset.touched !== "1") $("outName").value = v === "MOL" ? "single_molecule.lt" : `${v}.lt`;
  });
  $("outName").addEventListener("input", () => { $("outName").dataset.touched = "1"; });

  // ------------------------------------------------------------ 실행
  function run() {
    if (!R) return;
    try {
      result = P.run(R, ff, { smiles: $("smiles").value, name: $("molName").value, file });
    } catch (e) {
      if (!(e instanceof T.TypingError) && !(e instanceof window.DreidingIO.InputError)) console.error(e);
      result = null;
      showMsgs([], e.message);
      return;
    }
    showMsgs(result.warnings, null);
    render(result);
  }
  $("run").addEventListener("click", run);

  function showMsgs(warnings, error) {
    const box = $("msgs");
    box.textContent = "";
    if (error) {
      const d = document.createElement("div");
      d.className = "msg err";
      d.textContent = error;
      box.append(d);
    }
    for (const w of warnings) {
      const d = document.createElement("div");
      d.className = "msg warn";
      d.textContent = w;
      box.append(d);
    }
  }

  function ranges(nums) {
    const out = [];
    for (let k = 0; k < nums.length; k++) {
      let j = k;
      while (j + 1 < nums.length && nums[j + 1] === nums[j] + 1) j++;
      out.push(j - k >= 2 ? `${nums[k]}–${nums[j]}` : j > k ? `${nums[k]}, ${nums[j]}` : `${nums[k]}`);
      k = j;
    }
    return out.join(", ");
  }

  function render(r) {
    // 타입 순서: 원소 기호 순, 같은 원소 안에서는 처음 나온 순
    const firstSeen = new Map();
    r.rows.forEach((row, i) => { if (!firstSeen.has(row.label)) firstSeen.set(row.label, i); });
    const types = [...firstSeen.keys()].sort((a, b) => {
      const ea = r.rows[firstSeen.get(a)].sym, eb = r.rows[firstSeen.get(b)].sym;
      const rank = e => (e === "C" ? 0 : e === "H" ? 1 : 2);
      return rank(ea) - rank(eb) || ea.localeCompare(eb) || firstSeen.get(a) - firstSeen.get(b);
    });
    const flagged = new Set(r.flagged);
    viewer.setMolecule({ atoms: r.rows, bonds: r.bonds, flagged }, types);
    $("empty").hidden = true;
    $("toolbar").hidden = false;
    $("legend").hidden = flagged.size === 0;

    const grid = $("typeGrid");
    grid.textContent = "";
    for (const t of types) {
      const idx = r.rows.map((row, i) => (row.label === t ? i + 1 : 0)).filter(Boolean);
      const b = document.createElement("button");
      b.type = "button";
      b.className = "tcard";
      b.setAttribute("aria-pressed", "false");
      b.dataset.type = t;
      const sw = document.createElement("span");
      sw.className = "sw";
      sw.style.background = viewer.typeColor.get(t);
      const nm = document.createElement("span");
      nm.className = "name";
      nm.textContent = t;
      const ct = document.createElement("span");
      ct.className = "count";
      ct.textContent = `${idx.length}개`;
      const at = document.createElement("span");
      at.className = "atoms";
      at.textContent = `원자 ${ranges(idx)}`;
      b.append(sw, nm, ct, at);
      b.title = `질량 ${tables.labels[t]} — 눌러서 이 타입만 강조`;
      b.addEventListener("click", () => toggleType(t, b));
      grid.append(b);
    }
    $("typeSummary").textContent = `원자 ${r.rows.length}개, 타입 ${types.length}종, 조각 ${r.nFrags}개입니다. 타입을 누르면 해당 원자만 강조합니다.`;
    $("clearSel").hidden = true;

    const tb = $("atomRows");
    tb.textContent = "";
    r.rows.forEach((row, i) => {
      const tr = document.createElement("tr");
      if (flagged.has(i)) tr.className = "flag";
      const cells = [[i + 1, "num"], [row.name, ""], [row.label, "lab"], [row.q.toFixed(6), "num"], [row.frag, "num"],
        [row.xyz[0].toFixed(3), "num"], [row.xyz[1].toFixed(3), "num"], [row.xyz[2].toFixed(3), "num"]];
      for (const [v, cls] of cells) { const td = document.createElement("td"); td.textContent = v; if (cls) td.className = cls; tr.append(td); }
      tb.append(tr);
    });
    $("atomSummary").textContent = `원자별 목록 (${r.rows.length}개)`;
    $("atomDetails").hidden = false;

    $("lt").value = r.lt;
    $("ltTitle").textContent = `${outName()} (${r.name})`;
    $("copy").disabled = false;
    $("download").disabled = false;
    $("copied").textContent = "";
  }

  function toggleType(t, btn) {
    if (viewer.selected.has(t)) viewer.selected.delete(t); else viewer.selected.add(t);
    btn.setAttribute("aria-pressed", String(viewer.selected.has(t)));
    $("clearSel").hidden = viewer.selected.size === 0;
    viewer.draw();
  }
  $("clearSel").addEventListener("click", () => {
    viewer.selected.clear();
    document.querySelectorAll(".tcard").forEach(b => b.setAttribute("aria-pressed", "false"));
    $("clearSel").hidden = true;
    viewer.draw();
  });

  // ------------------------------------------------------------ 뷰어 도구
  document.querySelectorAll("[data-label]").forEach(b => b.addEventListener("click", () => {
    viewer.labelMode = b.dataset.label;
    document.querySelectorAll("[data-label]").forEach(x => x.setAttribute("aria-pressed", String(x === b)));
    viewer.draw();
  }));
  document.querySelectorAll("[data-color]").forEach(b => b.addEventListener("click", () => {
    viewer.colorMode = b.dataset.color;
    document.querySelectorAll("[data-color]").forEach(x => x.setAttribute("aria-pressed", String(x === b)));
    viewer.draw();
  }));
  $("resetView").addEventListener("click", () => viewer.resetView());
  $("png").addEventListener("click", () => {
    const src = $("view"), c = document.createElement("canvas");
    c.width = src.width; c.height = src.height;
    const ctx = c.getContext("2d");
    ctx.fillStyle = getComputedStyle(document.querySelector(".stage")).backgroundColor;
    ctx.fillRect(0, 0, c.width, c.height);
    ctx.drawImage(src, 0, 0);
    c.toBlob(b => save(b, outName().replace(/\.lt$/, "") + ".png"));
  });

  // ------------------------------------------------------------ 출력
  function outName() {
    let n = $("outName").value.trim() || "single_molecule.lt";
    if (!/\.lt$/i.test(n)) n += ".lt";
    return n.replace(/[\\/:*?"<>|]/g, "_");
  }
  function save(blob, name) {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = name;
    document.body.append(a);
    a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
  }
  $("download").addEventListener("click", () => {
    if (result) save(new Blob([result.lt], { type: "text/plain;charset=utf-8" }), outName());
  });
  $("copy").addEventListener("click", async () => {
    if (!result) return;
    let ok = false;
    try { await navigator.clipboard.writeText(result.lt); ok = true; } catch (e) {
      const ta = $("lt");
      ta.focus(); ta.select();
      try { ok = document.execCommand("copy"); } catch (e2) { ok = false; }
    }
    $("copied").textContent = ok ? "복사했습니다." : "복사하지 못했습니다. 아래 칸을 직접 선택해 복사하십시오.";
  });
})();
