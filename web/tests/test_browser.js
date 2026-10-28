// 빌드된 단일 HTML 을 file:// 로 열어 (네트워크 차단) 입력별 결과·복사·내려받기를 확인하는 Playwright 테스트
// 사용: node web/tests/test_browser.js [dreiding_typer.html]   (playwright 필요, 없으면 건너뜀)
const fs = require("fs");
const os = require("os");
const path = require("path");
const assert = require("assert");
const { execFileSync } = require("child_process");

let chromium;
try { ({ chromium } = require("playwright")); } catch (e) { console.log("playwright 가 없어 브라우저 테스트를 건너뜁니다"); process.exit(0); }

const HTML = path.resolve(process.argv[2] || path.join(__dirname, "../../dreiding_typer.html"));
const body = s => s.split("\n").slice(1).join("\n").replace(/^# 전하: .*$/m, "");

(async () => {
  const tmp = fs.mkdtempSync(path.join(os.tmpdir(), "dreiding-browser-"));
  execFileSync("python3", [path.join(__dirname, "make_cases.py"), tmp]);
  const cases = JSON.parse(fs.readFileSync(path.join(tmp, "cases.json")));
  const browser = await chromium.launch();
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 1000 }, acceptDownloads: true });
  await ctx.grantPermissions(["clipboard-read", "clipboard-write"]);
  const blocked = [];
  await ctx.route("**/*", r => (r.request().url().startsWith("file:") ? r.continue() : (blocked.push(r.request().url()), r.abort())));
  const page = await ctx.newPage();
  const errors = [];
  page.on("pageerror", e => errors.push(e.message));
  page.on("console", m => { if (m.type() === "error") errors.push(m.text()); });
  const load = async f => {
    await page.setInputFiles("#file", f);
    await page.waitForFunction(() => !document.getElementById("fileChip").hidden);
  };
  let failed = 0;
  const test = async (name, fn) => {
    try { await fn(); console.log("ok   ", name); } catch (e) { failed++; console.log("FAIL ", name, "\n     ", e.message.split("\n")[0]); }
  };

  await test("오프라인(file://)으로 RDKit 엔진 로딩", async () => {
    await page.goto("file://" + HTML);
    await page.waitForFunction(() => document.getElementById("status").classList.contains("ready"), null, { timeout: 60000 });
  });

  await test("SMILES 만으로 .lt 생성과 타입 카드", async () => {
    await page.click("text=아세트아닐라이드");
    const lt = await page.inputValue("#lt");
    assert.ok(lt.includes('import "dreiding.lt"') && lt.includes("MOL inherits DREIDING {"));
    assert.ok(lt.includes("@atom:N_R_b1_d2_hd") && lt.includes("@atom:C_2_b1"));
    assert.ok(await page.locator(".tcard").count() >= 6);
  });

  await test(`mol2 입력 .lt 가 Python CLI 와 같음 (${cases.length}개 계)`, async () => {
    for (const c of cases) {
      await page.fill("#smiles", c.smiles);
      await load(c.mol2);
      await page.click("#run");
      assert.strictEqual(body(await page.inputValue("#lt")), body(fs.readFileSync(c.lt, "utf8")), c.smiles);
      await page.click("#clearFile");
    }
  });

  await test("pdb·xyz 입력과 전하 경고", async () => {
    for (const key of ["pdb", "pdb_noconect", "xyz"]) {
      await page.fill("#smiles", cases[0].smiles);
      await load(cases[0][key]);
      await page.click("#run");
      assert.ok((await page.locator(".msg.warn").allTextContents()).some(t => t.includes("전하 정보가 없어")), key);
      await page.click("#clearFile");
    }
  });

  await test("구조가 다르면 오류 표시", async () => {
    await page.fill("#smiles", "CCO");
    await load(cases[0].mol2);
    await page.click("#run");
    assert.ok((await page.locator(".msg.err").textContent()).includes("원자/결합 수가 다릅니다"));
    await page.click("#clearFile");
  });

  await test("타입 강조·마우스 정보·복사·내려받기", async () => {
    await page.fill("#smiles", cases[0].smiles);
    await load(cases[0].mol2);
    await page.click("#run");
    await page.locator(".tcard", { hasText: "N_R_d2" }).click();
    assert.strictEqual(await page.locator(".tcard[aria-pressed=true]").count(), 1);
    await page.click("#clearSel");
    const box = await page.locator("#view").boundingBox();
    let shown = false;
    for (let gx = 0.1; gx <= 0.9 && !shown; gx += 0.02) {
      for (let gy = 0.2; gy <= 0.8 && !shown; gy += 0.04) {
        await page.mouse.move(box.x + box.width * gx, box.y + box.height * gy);
        shown = await page.locator("#tip").isVisible();
      }
    }
    assert.ok(shown && /원자 \d+/.test(await page.textContent("#tip")), "마우스 정보");
    await page.click("#copy");
    assert.strictEqual(await page.evaluate(() => navigator.clipboard.readText()), await page.inputValue("#lt"));
    const [dl] = await Promise.all([page.waitForEvent("download"), page.click("#download")]);
    assert.strictEqual(dl.suggestedFilename(), "single_molecule.lt");
    assert.strictEqual(fs.readFileSync(await dl.path(), "utf8"), await page.inputValue("#lt"));
  });

  await test("모바일 폭에서 가로 스크롤 없음", async () => {
    await page.setViewportSize({ width: 390, height: 900 });
    assert.ok(await page.evaluate(() => document.documentElement.scrollWidth) <= 390);
  });

  await test("외부 네트워크 요청 없음, 페이지 오류 없음", async () => {
    assert.deepStrictEqual(blocked, []);
    assert.deepStrictEqual(errors, []);
  });

  await browser.close();
  fs.rmSync(tmp, { recursive: true, force: true });
  console.log(failed ? `\n${failed} failed` : "\nall passed");
  process.exit(failed ? 1 : 0);
})();
