#!/usr/bin/env python3
"""H6-U5, parte A. Conductor de mediciones: vuelve a ejecutar, de uno en uno
y SIN modificarlos, los programas sellados de las unidades 2, 3 y 4, con
salidas nuevas en la carpeta `mediciones`. Comprueba la reproducibilidad
byte a byte de las salidas deterministas y resume tiempos y memoria.

Uso: python -I -B medir_h6_u5.py --output <nuevo.json>
Los tiempos y la memoria no son deterministas: este JSON no se sella.
"""
import argparse
import hashlib
import json
import os
import statistics
import subprocess
import sys
import time

T0 = time.perf_counter()
HERE = os.path.dirname(os.path.abspath(__file__))
COMP = os.path.dirname(HERE)
D2 = os.path.join(COMP, "H6_U2_CIRCUITO_PATH25_MULTICOMPONENTE_20261112")
D3 = os.path.join(COMP, "H6_U3_ESTIMADOR_CERTIFICADO_CIRCUITO_PATH25_20261112")
D4 = os.path.join(COMP, "H6_U4_ESTIMADOR_CERTIFICADO_FENOMENOLOGICO_PATH25_20261112")
MED = os.path.join(HERE, "mediciones")
PROGRAMS = {
    os.path.join(D2, "verificar_h6_u2.py"): "66ab5ab9e275c2a3daca1246b128a531ed9d89ee7d40fdeb274210d085da28d0",
    os.path.join(D3, "estimador_h6_u3.py"): "6b6812fb1a315ae40f0a75d3b8f12d94fcd3bfdc417bc72936e0cabf79682ac2",
    os.path.join(D3, "contrastar_h6_u3.py"): "e3205070b8952d3129c879fe64aefc16ad9bb987e9cd14041d5337fd02f1c7ea",
    os.path.join(D4, "estimador_h6_u4.py"): "620f8bea31171c1f2f7b7bd67f2415a2350018fb79cf1ab32f342bba893546b2",
    os.path.join(D4, "contrastar_h6_u4.py"): "ddd1c23069772afb5acc3a6498b094ea5e06bd46d0dbee78f54d8481131fa549",
}
SEALED = {
    "RESULTADO_H6_U2": "e23905fb361a5a5dea8cc5414461182b1faf909b33d3b0979215bf9ff1347250",
    "ESTIMACION_H6_U3": "349658833be721cbe617ac7f0948d5b25184a46b2fc2b6d62a23790eafad8500",
    "RESULTADO_H6_U3": "fe714cddb5dc5b9c71b27ccf0e53e51749940e1b6bf7bbcccb523c4efaf6a3ff",
    "ESTIMACION_H6_U4": "654367777c23c83cdc6cccbbdb756544f85cca57ad0e86a680d6deec473304e9",
    "RESULTADO_H6_U4": "027973eb8c4a6cdacb2f1478a505ed509dfcef2cd951f8fe6bb7ea353902cb6c",
}
REPS = 5
MAX_RUNS = 30
WALL_CAP = 580.0
CHECKS = []
RUNS = []


class Stop(Exception):
    pass


def check(name, cond, detail=None):
    row = {"name": name, "pass": bool(cond)}
    if detail is not None:
        row["detail"] = detail
    CHECKS.append(row)
    if not cond:
        raise Stop(name)


def sha_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def launch(label, script, argv):
    if len(RUNS) >= MAX_RUNS or time.perf_counter() - T0 > WALL_CAP:
        raise Stop("techo de ejecuciones o de tiempo")
    t = time.perf_counter()
    pr = subprocess.run([sys.executable, "-I", "-B", script] + argv, capture_output=True, text=True)
    wall = time.perf_counter() - t
    RUNS.append({"etiqueta": label, "programa": os.path.basename(script), "rc": pr.returncode,
                 "segundos_proceso": round(wall, 4), "salida": pr.stdout.strip()[-200:]})
    if pr.returncode != 0:
        CHECKS.append({"name": "ejecucion " + label, "pass": False, "detail": (pr.stdout + pr.stderr)[-600:]})
        raise Stop("ejecucion fallida: " + label)
    return wall


def stats(vals):
    return {"mediana": round(statistics.median(vals), 4), "min": round(min(vals), 4), "max": round(max(vals), 4),
            "n": len(vals)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    if os.path.exists(args.output) or os.path.exists(MED):
        raise SystemExit("La salida o la carpeta de mediciones ya existen.")
    os.makedirs(MED)
    out = {"unit": "H6-U5", "programa": "mediciones",
           "nota": "tiempos de reloj y memoria no deterministas; los contadores salen de las salidas selladas"}
    status = "PASS"
    try:
        run(out)
    except Stop as ex:
        status = "FAIL"
        out["stopped_at"] = str(ex)
    out["status"] = status
    out["checks"] = CHECKS
    out["ejecuciones"] = RUNS
    out["ejecuciones_totales"] = len(RUNS)
    out["segundos_totales"] = round(time.perf_counter() - T0, 2)
    with open(args.output, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(out, ensure_ascii=False, indent=1, sort_keys=True) + "\n")
    print("%s mediciones ejecuciones=%d controles=%d fallidos=%d t=%.1fs"
          % (status, len(RUNS), len(CHECKS), sum(1 for c in CHECKS if not c["pass"]), out["segundos_totales"]))
    if status != "PASS":
        print("DETENIDO EN:", out.get("stopped_at"))
        sys.exit(1)


def run(out):
    check("P0 los cinco programas sellados no han cambiado", all(sha_file(p) == h for p, h in PROGRAMS.items()))
    j = lambda *a: os.path.join(*a)
    sealed3 = json.loads(open(j(D3, "ESTIMACION_H6_U3.json"), encoding="utf-8").read())
    sealed4 = json.loads(open(j(D4, "ESTIMACION_H6_U4.json"), encoding="utf-8").read())
    summary = {}

    def estimator_group(tag, script, kmax, sealed):
        walls, inner, mem, cat, gas = [], [], [], [], []
        same = True
        for i in range(REPS):
            base = j(MED, "%s_k%d_%d" % (tag, kmax, i))
            walls.append(launch("%s k=%d #%d" % (tag, kmax, i), script,
                                ["--output", base + ".json", "--sello", base + "_sello.json",
                                 "--recursos", base + "_rec.json", "--kmax", str(kmax)]))
            rec = json.loads(open(base + "_rec.json", encoding="utf-8").read())
            tm = rec["tiempos"]
            inner.append(tm["total"])
            mem.append(rec["pico_memoria_bytes"])
            cat.append(sum(v for k, v in tm.items() if k.startswith("catalogo_")))
            gas.append(sum(v for k, v in tm.items() if k.startswith("gas_")))
            got = json.loads(open(base + ".json", encoding="utf-8").read())
            same = same and got["status"] == "PASS" and got["tabla"] == sealed["tabla"][:kmax]
        check("M2 %s con corte k=%d: las filas coinciden con la tabla sellada en las %d repeticiones"
              % (tag, kmax, REPS), same)
        summary["%s_k%d" % (tag, kmax)] = {
            "segundos_proceso": stats(walls), "segundos_internos": stats(inner),
            "segundos_catalogo": stats(cat), "segundos_gas": stats(gas),
            "pico_memoria_MiB": stats([m / 1048576 for m in mem])}

    estimator_group("u3_estimador", j(D3, "estimador_h6_u3.py"), 2, sealed3)
    estimator_group("u4_estimador", j(D4, "estimador_h6_u4.py"), 3, sealed4)
    estimator_group("u4_estimador", j(D4, "estimador_h6_u4.py"), 4, sealed4)

    def contrast_group(tag, ddir, script, sealed_name, phases):
        walls, ref, mem = [], [], []
        same = True
        for i in range(REPS):
            base = j(MED, "%s_%d" % (tag, i))
            walls.append(launch("%s #%d" % (tag, i), j(ddir, script),
                                ["--estimacion", j(ddir, sealed_name[0]), "--sello", j(ddir, "SELLO_ESTIMACION.json"),
                                 "--recursos-estimador", j(ddir, "RECURSOS_ESTIMADOR.json"),
                                 "--output", base + ".json", "--recursos", base + "_rec.json"]))
            rec = json.loads(open(base + "_rec.json", encoding="utf-8").read())
            ref.append(sum(rec["tiempos"][p] for p in phases))
            mem.append(rec["pico_memoria_bytes"])
            same = same and sha_file(base + ".json") == SEALED[sealed_name[1]]
        check("M1 %s: el resultado del contraste es identico byte a byte al sellado en las %d repeticiones"
              % (tag, REPS), same)
        summary[tag] = {"segundos_proceso": stats(walls), "segundos_referencia_exacta": stats(ref),
                        "pico_memoria_MiB": stats([m / 1048576 for m in mem])}

    contrast_group("u3_contraste", D3, "contrastar_h6_u3.py", ("ESTIMACION_H6_U3.json", "RESULTADO_H6_U3"),
                   ("referencia_preparacion", "referencia_enumeracion"))
    contrast_group("u4_contraste", D4, "contrastar_h6_u4.py", ("ESTIMACION_H6_U4.json", "RESULTADO_H6_U4"),
                   ("referencia_tabla_decoder", "referencia_enumeradores"))

    # reejecuciones completas
    w2 = launch("u2 verificador completo", j(D2, "verificar_h6_u2.py"),
                ["--output", j(MED, "u2_completo.json"), "--tiempos", j(MED, "u2_completo_tiempos.txt")])
    check("M1 unidad 2: RESULTADO reproducido byte a byte", sha_file(j(MED, "u2_completo.json")) == SEALED["RESULTADO_H6_U2"])
    b3 = j(MED, "u3_estimador_k8_completo")
    w3 = launch("u3 estimador completo", j(D3, "estimador_h6_u3.py"),
                ["--output", b3 + ".json", "--sello", b3 + "_sello.json", "--recursos", b3 + "_rec.json"])
    check("M1 unidad 3: ESTIMACION reproducida byte a byte", sha_file(b3 + ".json") == SEALED["ESTIMACION_H6_U3"])
    b4 = j(MED, "u4_estimador_k5_completo")
    w4 = launch("u4 estimador completo", j(D4, "estimador_h6_u4.py"),
                ["--output", b4 + ".json", "--sello", b4 + "_sello.json", "--recursos", b4 + "_rec.json"])
    check("M1 unidad 4: ESTIMACION reproducida byte a byte", sha_file(b4 + ".json") == SEALED["ESTIMACION_H6_U4"])
    r3 = json.loads(open(b3 + "_rec.json", encoding="utf-8").read())
    r4 = json.loads(open(b4 + "_rec.json", encoding="utf-8").read())
    summary["reejecuciones_completas"] = {
        "u2_verificador_segundos_proceso": round(w2, 3),
        "u3_estimador_k8": {"segundos_proceso": round(w3, 3), "segundos_internos": r3["tiempos"]["total"],
                            "pico_memoria_MiB": round(r3["pico_memoria_bytes"] / 1048576, 1)},
        "u4_estimador_k5": {"segundos_proceso": round(w4, 3), "segundos_internos": r4["tiempos"]["total"],
                            "pico_memoria_MiB": round(r4["pico_memoria_bytes"] / 1048576, 1)},
    }
    out["resumen"] = summary


if __name__ == "__main__":
    main()
