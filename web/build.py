# 규칙표·JS 모듈·RDKit.js(압축 wasm)를 템플릿에 끼워 넣어 오프라인 단일 HTML 을 만드는 빌드 스크립트
"""사용 (저장소 최상위에서):
    cd web && npm install && cd ..
    python3 web/build.py                 # → dreiding_typer.html

dreiding.lt 는 설치된 moltemplate 의 것을 읽는다(--dreiding-lt 로 바꿀 수 있음).
"""
import argparse
import base64
import gzip
import os
import subprocess
import sys

WEB = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(WEB)


def main():
    ap = argparse.ArgumentParser(description="오프라인 단일 HTML 빌드")
    ap.add_argument("-o", "--out", default=os.path.join(ROOT, "dreiding_typer.html"))
    ap.add_argument("--dreiding-lt", default=None, help="규칙을 읽을 dreiding.lt (기본: 설치된 moltemplate)")
    ap.add_argument("--rdkit-dist", default=os.path.join(WEB, "node_modules", "@rdkit", "rdkit", "dist"))
    args = ap.parse_args()

    js_path = os.path.join(args.rdkit_dist, "RDKit_minimal.js")
    wasm_path = os.path.join(args.rdkit_dist, "RDKit_minimal.wasm")
    if not os.path.exists(wasm_path):
        sys.exit(f"{wasm_path} 가 없습니다. web 폴더에서 npm install 을 먼저 하십시오")

    cmd = [sys.executable, os.path.join(WEB, "build_tables.py")] + ([args.dreiding_lt] if args.dreiding_lt else [])
    tables = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    with open(wasm_path, "rb") as f:
        wasm_b64 = base64.b64encode(gzip.compress(f.read(), 9, mtime=0)).decode("ascii")

    def src(name, check=True):
        with open(os.path.join(WEB, "src", name), encoding="utf-8") as f:
            text = f.read()
        if check and "</script" in text.lower():
            sys.exit(f"{name} 에 </script 문자열이 있어 인라인할 수 없습니다")
        return text

    with open(js_path, encoding="utf-8") as f:
        rdkit_js = f.read()
    parts = {
        "/*TABLES*/": tables,
        "/*RDKIT_WASM_GZ_B64*/": wasm_b64,
        "/*RDKIT_JS*/": rdkit_js,
        "/*TYPER_JS*/": src("typer.js"),
        "/*IO_JS*/": src("io.js"),
        "/*GEOM_JS*/": src("geom.js"),
        "/*PIPELINE_JS*/": src("pipeline.js"),
        "/*VIEWER_JS*/": src("viewer.js"),
        "/*APP_JS*/": src("app.js"),
    }
    html = src("index.html", check=False)
    for key, val in parts.items():
        if html.count(key) != 1:
            sys.exit(f"템플릿에 {key} 자리가 정확히 하나 있어야 합니다")
        html = html.replace(key, val)
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"{args.out}: {os.path.getsize(args.out) / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
