#!/usr/bin/env python3
"""Reproduce the deterministic outputs of the certified PATH25 estimators.

Usage (from any directory):
    python -I -B reproduce_h6.py --out <new directory outside this archive>

The script first checks MANIFEST.sha256.txt. It then runs, one at a time
and without modifying them, the stored programs below, writing every new
file into --out, and compares each deterministic output byte by byte with
the stored one:

  1. verificar_h6_u2.py             -> RESULTADO_H6_U2.json
  2. estimador_h6_u3.py  (k <= 8)   -> ESTIMACION_H6_U3.json
  3. contrastar_h6_u3.py            -> RESULTADO_H6_U3.json
  4. estimador_h6_u4.py  (k <= 5)   -> ESTIMACION_H6_U4.json
  5. contrastar_h6_u4.py            -> RESULTADO_H6_U4.json
  6. control_v14_fenomenologico.py  -> CONTROL_V14_FENOMENOLOGICO.json

Steps 3 and 5 run on the stored, sealed estimate, as in the original
order (estimate first, reference afterwards); their result records the
seal time, so they are compared with the stored results. A fresh estimate
from steps 2 and 4 is byte-identical to the stored one when reproduction
succeeds. Seal files, resource records and times are not deterministic
and are not compared. Python standard library only.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
COMP = os.path.join(ROOT, "qec_revision_matematica", "05_COMPROBACIONES")
D2 = os.path.join(COMP, "H6_U2_CIRCUITO_PATH25_MULTICOMPONENTE_20261112")
D3 = os.path.join(COMP, "H6_U3_ESTIMADOR_CERTIFICADO_CIRCUITO_PATH25_20261112")
D4 = os.path.join(COMP, "H6_U4_ESTIMADOR_CERTIFICADO_FENOMENOLOGICO_PATH25_20261112")
D5 = os.path.join(COMP, "H6_U5_RECURSOS_G4_Y_MATRIZ_FINAL_20261112")


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def check_manifest():
    bad = []
    with open(os.path.join(ROOT, "MANIFEST.sha256.txt"), encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            digest, rel = line.split("  ", 1)
            p = os.path.join(ROOT, *rel.split("/"))
            if not os.path.isfile(p) or sha(p) != digest:
                bad.append(rel)
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="new directory outside the extracted archive")
    args = ap.parse_args()
    out = os.path.abspath(args.out)
    if os.path.exists(out):
        raise SystemExit("--out already exists; give a new directory.")
    try:
        inside = os.path.commonpath([out, ROOT]) == ROOT
    except ValueError:  # different drives on Windows
        inside = False
    if inside:
        raise SystemExit("--out must be outside the extracted archive.")
    bad = check_manifest()
    if bad:
        raise SystemExit("Manifest mismatch: " + ", ".join(bad[:10]))
    os.makedirs(out)
    o = lambda name: os.path.join(out, name)
    steps = [
        ("verificar_h6_u2", os.path.join(D2, "verificar_h6_u2.py"),
         ["--output", o("RESULTADO_H6_U2.json"), "--tiempos", o("TIEMPOS_U2.txt")],
         [(o("RESULTADO_H6_U2.json"), os.path.join(D2, "RESULTADO_H6_U2.json"))]),
        ("estimador_h6_u3", os.path.join(D3, "estimador_h6_u3.py"),
         ["--output", o("ESTIMACION_H6_U3.json"), "--sello", o("SELLO_U3.json"),
          "--recursos", o("RECURSOS_ESTIMADOR_U3.json")],
         [(o("ESTIMACION_H6_U3.json"), os.path.join(D3, "ESTIMACION_H6_U3.json"))]),
        ("contrastar_h6_u3", os.path.join(D3, "contrastar_h6_u3.py"),
         ["--estimacion", os.path.join(D3, "ESTIMACION_H6_U3.json"),
          "--sello", os.path.join(D3, "SELLO_ESTIMACION.json"),
          "--recursos-estimador", os.path.join(D3, "RECURSOS_ESTIMADOR.json"),
          "--output", o("RESULTADO_H6_U3.json"), "--recursos", o("RECURSOS_CONTRASTE_U3.json")],
         [(o("RESULTADO_H6_U3.json"), os.path.join(D3, "RESULTADO_H6_U3.json"))]),
        ("estimador_h6_u4", os.path.join(D4, "estimador_h6_u4.py"),
         ["--output", o("ESTIMACION_H6_U4.json"), "--sello", o("SELLO_U4.json"),
          "--recursos", o("RECURSOS_ESTIMADOR_U4.json")],
         [(o("ESTIMACION_H6_U4.json"), os.path.join(D4, "ESTIMACION_H6_U4.json"))]),
        ("contrastar_h6_u4", os.path.join(D4, "contrastar_h6_u4.py"),
         ["--estimacion", os.path.join(D4, "ESTIMACION_H6_U4.json"),
          "--sello", os.path.join(D4, "SELLO_ESTIMACION.json"),
          "--recursos-estimador", os.path.join(D4, "RECURSOS_ESTIMADOR.json"),
          "--output", o("RESULTADO_H6_U4.json"), "--recursos", o("RECURSOS_CONTRASTE_U4.json")],
         [(o("RESULTADO_H6_U4.json"), os.path.join(D4, "RESULTADO_H6_U4.json"))]),
        ("control_v14_fenomenologico", os.path.join(D5, "control_v14_fenomenologico.py"),
         ["--output", o("CONTROL_V14_FENOMENOLOGICO.json")],
         [(o("CONTROL_V14_FENOMENOLOGICO.json"), os.path.join(D5, "CONTROL_V14_FENOMENOLOGICO.json"))]),
    ]
    summary = {"python": sys.version.split()[0], "steps": [], "status": "PASS"}
    for label, script, argv, compare in steps:
        t = time.perf_counter()
        pr = subprocess.run([sys.executable, "-I", "-B", script] + argv, capture_output=True, text=True)
        row = {"step": label, "returncode": pr.returncode, "seconds": round(time.perf_counter() - t, 2),
               "stdout_tail": pr.stdout.strip()[-200:], "identical": []}
        if pr.returncode != 0:
            row["stderr_tail"] = pr.stderr.strip()[-600:]
            summary["status"] = "FAIL"
        for new, stored in compare:
            same = os.path.isfile(new) and sha(new) == sha(stored)
            row["identical"].append({"file": os.path.basename(stored), "byte_identical": same})
            if not same:
                summary["status"] = "FAIL"
        summary["steps"].append(row)
        print("%-28s rc=%d  %6.2f s  %s" % (label, pr.returncode, row["seconds"],
                                         " ".join("%s=%s" % (c["file"], "identical" if c["byte_identical"]
                                                             else "DIFFERENT") for c in row["identical"])))
        if summary["status"] != "PASS":
            break
    with open(o("REPRODUCTION_SUMMARY.json"), "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(summary, indent=1, sort_keys=True) + "\n")
    print(summary["status"])
    sys.exit(0 if summary["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
