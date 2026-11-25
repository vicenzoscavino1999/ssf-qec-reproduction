#!/usr/bin/env python3
"""H6-U4. Contraste de la estimacion sellada con una referencia exacta global.

Solo arranca si el sello existe y coincide. Referencia: para cada sindrome
observado s, H_X R = 0 equivale a D = sindrome terminal de s, y entonces E
recorre C(s) + ker H_X. Se acumulan enumeradores enteros por (peso de D,
peso de E) de las clases C(s)+S_X y C(s)+L+S_X: valen para todo (p,q).
Anclajes: fracciones exactas publicadas en el paquete de observables y
espectros q=0 del piloto PATH25. SOLO DESPUES se abre la estimacion.

Uso: python -I -B contrastar_h6_u4.py --estimacion <json> --sello <json>
            --recursos-estimador <json> --output <nuevo.json> --recursos <nuevo.json>
"""
import argparse
import csv
import hashlib
import importlib.util
import io
import json
import os
import sys
import time
import zipfile
from decimal import Decimal, localcontext
from fractions import Fraction as F

T0 = time.perf_counter()
HERE = os.path.dirname(os.path.abspath(__file__))
COMP = os.path.dirname(HERE)
QEC = os.path.join(os.path.dirname(os.path.dirname(COMP)), "qec")
U2 = os.path.join(COMP, "H6_U2_CIRCUITO_PATH25_MULTICOMPONENTE_20261112", "verificar_h6_u2.py")
PILOT_REF = os.path.join(COMP, "PILOTO_CONECTADO_PATH25_20261028", "attempt_01", "REFERENCE.json")
QEM_ZIP = os.path.join(QEC, "SSF_QEM_extension_v2.zip")
QEM_MEMBER = "qem_extension_v2/results/path25/path25_metrics.csv"
PIN = {
    U2: "66ab5ab9e275c2a3daca1246b128a531ed9d89ee7d40fdeb274210d085da28d0",
    PILOT_REF: "397a90d292e6cea5ad70f6d182b59d884e74584a2bf96d203e68b38f2b311dc7",
    QEM_ZIP: "21e4327ae7818dcc8228d93a3861b43b1164bfbeb62e07e4390b96917b0ceb5f",
}
PAIRS = ((F(1, 50), F(1, 100)), (F(1, 100), F(1, 200)), (F(1, 100), F(0)))
N, M = 25, 12
LX = 15
CHECKS = []
TIMES = {}


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


def key(pair):
    return "%s|%s" % (pair[0], pair[1])


def fr(d):
    a, b = d["exacto"].split("/")
    return F(int(a), int(b))


def show(x):
    with localcontext() as ctx:
        ctx.prec = 30
        return {"exacto": "%d/%d" % (x.numerator, x.denominator),
                "decimal": format(Decimal(x.numerator) / Decimal(x.denominator), ".12E")}


def peak_memory_bytes():
    try:
        import ctypes
        from ctypes import wintypes

        class PMC(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                        ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]
        k32 = ctypes.windll.kernel32
        k32.GetCurrentProcess.restype = wintypes.HANDLE
        k32.K32GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(PMC), wintypes.DWORD]
        pmc = PMC()
        pmc.cb = ctypes.sizeof(PMC)
        if k32.K32GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb):
            return int(pmc.PeakWorkingSetSize)
    except Exception:
        pass
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--estimacion", required=True)
    ap.add_argument("--sello", required=True)
    ap.add_argument("--recursos-estimador", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--recursos", required=True)
    args = ap.parse_args()
    for p in (args.output, args.recursos):
        if os.path.exists(p):
            raise SystemExit("La salida ya existe; use un nombre nuevo: " + p)
    out = {"unit": "H6-U4", "programa": "contraste"}
    status = "PASS"
    try:
        run(args, out)
    except Stop as ex:
        status = "FAIL"
        out["stopped_at"] = str(ex)
    out["status"] = status
    out["checks_total"] = len(CHECKS)
    out["checks_failed"] = sum(1 for c in CHECKS if not c["pass"])
    out["checks"] = CHECKS
    with open(args.output, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(out, ensure_ascii=False, indent=1, sort_keys=True) + "\n")
    TIMES["total"] = round(time.perf_counter() - T0, 4)
    with open(args.recursos, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps({"nota": "tiempos de reloj en segundos y memoria; no deterministas; no son resultado",
                            "tiempos": TIMES, "pico_memoria_bytes": peak_memory_bytes()},
                           ensure_ascii=False, indent=1) + "\n")
    print("%s contraste controles=%d fallidos=%d t=%.2fs" % (status, out["checks_total"], out["checks_failed"],
                                                           TIMES["total"]))
    if status != "PASS":
        print("DETENIDO EN:", out.get("stopped_at"))
        sys.exit(1)


def run(args, out):
    # ------------------------------------------------ 0. sello (sin abrir el contenido)
    seal = json.loads(open(args.sello, encoding="utf-8").read())
    est_sha = sha_file(args.estimacion)
    check("S1 la estimacion existe, esta sellada, su hash coincide con el sello y el estimador no ha cambiado",
          seal["sha256"] == est_sha and seal["status"] == "PASS"
          and seal["estimador_sha256"] == sha_file(os.path.join(HERE, "estimador_h6_u4.py")),
          {"sellado": seal["sellado"], "sha256": est_sha, "k_completados": seal.get("k_completados")})
    check("S2 archivos de anclaje sin cambios (verificador de la unidad 2, espectros del piloto, paquete de "
          "observables)", all(sha_file(p) == h for p, h in PIN.items()))
    out["sello"] = seal

    # ------------------------------------------------ 1. referencia exacta global
    t = time.perf_counter()
    spec = importlib.util.spec_from_file_location("u2_sellado", U2)
    u2 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(u2)
    hx, hz, cols, actions, foot = u2.build_masks()
    T1 = u2.transitions_bits(foot)
    COR, END, TR, TOUCH = u2.resolve(T1, actions)
    bz = u2.rank_basis(hz)
    check("R0 H_Z tiene rango 12 y L=15 es un logico: sindrome nulo y fuera de S_X",
          len(bz) == 12 and u2.lin(cols, LX) == 0 and not u2.in_span(LX, bz))
    SX = [0]
    for g in hz:
        SX += [x ^ g for x in SX]
    check("R0b S_X tiene 4096 elementos distintos", len(SX) == 4096 and len(set(SX)) == 4096)
    TIMES["referencia_tabla_decoder"] = round(time.perf_counter() - t, 4)

    t = time.perf_counter()
    AS = [[0] * (N + 1) for _ in range(M + 1)]     # exito:   peso de D x peso de E
    AL = [[0] * (N + 1) for _ in range(M + 1)]     # fallo silencioso
    memo = {}
    evals = 0
    for s in range(1 << M):
        c = COR[s]
        tw = bin(END[s]).count("1")
        hs = memo.get(c)
        if hs is None:
            h1, h2 = [0] * (N + 1), [0] * (N + 1)
            cl = c ^ LX
            for g in SX:
                h1[(c ^ g).bit_count()] += 1
                h2[(cl ^ g).bit_count()] += 1
            evals += 2 * len(SX)
            hs = (h1, h2)
            memo[c] = hs
        rs, rl = AS[tw], AL[tw]
        for w in range(N + 1):
            rs[w] += hs[0][w]
            rl[w] += hs[1][w]
    TIMES["referencia_enumeradores"] = round(time.perf_counter() - t, 4)
    ref_cost = {"sindromes_en_tabla": 1 << M, "scores_tabla_completa": (1 << M) * len(actions),
                "evaluaciones_de_peso": evals, "correcciones_distintas": len(memo)}

    def evaluate(p, q):
        succ = sil = F(0)
        for tw in range(M + 1):
            qf = q ** tw * (1 - q) ** (M - tw)
            if qf == 0:
                continue
            for w in range(N + 1):
                pf = p ** w * (1 - p) ** (N - w)
                succ += AS[tw][w] * pf * qf
                sil += AL[tw][w] * pf * qf
        return {"P_0": succ + sil, "P_sil": sil, "Q_sil": sil / (succ + sil)}

    exact = {key(pr): evaluate(*pr) for pr in PAIRS}
    out["referencia_exacta"] = {k: {n: show(v) for n, v in d.items()} for k, d in exact.items()}

    # anclaje 1: fracciones exactas publicadas (paquete de observables)
    with zipfile.ZipFile(QEM_ZIP) as zf:
        rows = list(csv.DictReader(io.StringIO(zf.read(QEM_MEMBER).decode("utf-8"))))
    pub = {}
    for r in rows:
        pub.setdefault("%s|%s" % (r["p"], r["q"]), r)

    def frs(sv):
        a, _, b = sv.partition("/")
        return F(int(a), int(b or 1))

    bad = []
    for k_, d in exact.items():
        r = pub.get(k_)
        if r is None or frs(r["P0"]) != d["P_0"] or frs(r["Psil"]) != d["P_sil"] or frs(r["Qsil"]) != d["Q_sil"]:
            bad.append(k_)
    check("R1 la referencia coincide como fraccion exacta con los valores publicados (P_0, P_sil, Q_sil) en los "
          "tres pares de tasas", not bad, bad)
    # anclaje 2: espectros q=0 del piloto
    sp = json.loads(open(PILOT_REF, encoding="utf-8").read())["spectra"]
    check("R2 con D=0 los enumeradores reproducen los espectros por peso del piloto (sindrome nulo y fallo "
          "silencioso)",
          [AS[0][w] + AL[0][w] for w in range(N + 1)] == sp["zero_syndrome"] and AL[0] == sp["silent_failure"])

    # ------------------------------------------------ 2. AHORA se abre la estimacion
    est = json.loads(open(args.estimacion, encoding="utf-8").read())
    check("S3 la estimacion declara que no leyo referencia ni resultados de otras unidades",
          est["lee_referencia_o_resultados_de_otras_unidades"] is False and est["status"] == "PASS")
    table = est["tabla"]
    out["k_completados"] = est["k_completados"]
    out["estimador_truncado_por_techo_de_tiempo"] = est["truncado_por_techo_de_tiempo"]
    rows_out = []
    bad_contain, bad_err = [], []
    for row in table:
        k = row["k"]
        for pr in PAIRS:
            cell = row["tasas"][key(pr)]
            for nm in ("P_0", "P_sil", "Q_sil"):
                lo, hi = fr(cell[nm]["L"]), fr(cell[nm]["U"])
                ex = exact[key(pr)][nm]
                ok = lo <= ex <= hi
                if not ok:
                    bad_contain.append([k, key(pr), nm])
                item = {"k": k, "par": key(pr), "cantidad": nm, "L": cell[nm]["L"]["decimal"],
                        "U": cell[nm]["U"]["decimal"], "exacto": show(ex)["decimal"], "contenido": ok}
                if cell[nm]["valor_informado"] is not None:
                    val = fr(cell[nm]["valor_informado"])
                    gar = fr(cell[nm]["error_relativo_garantizado"])
                    real = abs(val - ex) / ex
                    item["error_garantizado"] = show(gar)["decimal"]
                    item["error_real"] = show(real)["decimal"]
                    if real > gar:
                        bad_err.append([k, key(pr), nm])
                rows_out.append(item)
    check("E1 los valores exactos de P_0, P_sil y Q_sil estan dentro de sus intervalos para todo k completado y "
          "los tres pares", not bad_contain, bad_contain)
    check("E3 el error real del valor informado no supera el error garantizado (todo k con L>0)",
          not bad_err, bad_err)
    mono = True
    for pr in PAIRS:
        Ls = [fr(r["tasas"][key(pr)]["P_sil"]["L"]) for r in table]
        ep = [fr(r["tasas"][key(pr)]["eps_k"]) for r in table]
        mono = mono and all(Ls[i] <= Ls[i + 1] for i in range(len(Ls) - 1)) \
            and all(ep[i] >= ep[i + 1] for i in range(len(ep) - 1))
    check("E4 L_k no decrece y eps_k no crece con k", mono)
    ks = est["k_estrella"]
    inf = est["estimacion_informada"]
    main_key = key(PAIRS[0])
    e5 = ks[main_key] is not None and \
        fr(inf[main_key]["certificado"]["P_sil"]["error_relativo_garantizado"]) <= F(1, 10)
    CHECKS.append({"name": "E5 en el par principal existe k* y cumple el objetivo: error garantizado de P_sil "
                           "<= 10 % (no alcanzarlo se informa; no detiene)", "pass": bool(e5), "detail": ks})
    out["objetivo_alcanzado_en_el_par_principal"] = bool(e5)

    # ------------------------------------------------ 3. tablas de salida
    out["contraste_por_k"] = rows_out
    informed = {}
    for pr in PAIRS:
        kk = key(pr)
        if ks[kk] is None:
            informed[kk] = None
            continue
        cert = inf[kk]["certificado"]
        informed[kk] = {"k": ks[kk], "polimeros": inf[kk]["polimeros"], "familias": inf[kk]["familias"],
                        "estados_del_gas": inf[kk]["estados_del_gas"], "terminos_del_gas": inf[kk]["terminos_del_gas"],
                        "recursos_acumulados": inf[kk]["recursos_acumulados"], "eps_k": cert["eps_k"]["decimal"]}
        for nm in ("P_0", "P_sil", "Q_sil"):
            val, gar = fr(cert[nm]["valor_informado"]), fr(cert[nm]["error_relativo_garantizado"])
            ex = exact[kk][nm]
            informed[kk][nm] = {"L": cert[nm]["L"]["decimal"], "U": cert[nm]["U"]["decimal"],
                                "valor_informado": show(val)["decimal"], "exacto": show(ex)["decimal"],
                                "error_garantizado": show(gar)["decimal"],
                                "error_real": show(abs(val - ex) / ex)["decimal"],
                                "anchura_relativa": cert[nm]["anchura_relativa"]["decimal"]}
    out["estimacion_informada_frente_a_referencia"] = informed
    rec_est = json.loads(open(args.recursos_estimador, encoding="utf-8").read())
    tm = rec_est["tiempos"]
    per_k = []
    cum = tm.get("preparacion", 0.0)
    for row in table:
        k = row["k"]
        cum += tm.get("catalogo_tamano_%d" % k, 0.0)
        per_k.append({"k": k, "historias_locales": row["acumulado"]["ejecuciones_locales"],
                      "sindromes_distintos": row["acumulado"]["sindromes_distintos"],
                      "scores": row["acumulado"]["scores"], "pruebas_conexidad": row["acumulado"]["pruebas_conexidad"],
                      "decoraciones": row["acumulado"]["decoraciones"], "polimeros": row["polimeros"],
                      "familias": row["familias"], "historias_retenidas": row["historias_retenidas"],
                      "estados_del_gas": row["estados_del_gas"], "terminos_del_gas": row["terminos_del_gas"],
                      "consultas_por_subconjunto": row["consultas_por_subconjunto"],
                      "bits_max": max(c["bits_max"] for c in row["tasas"].values()),
                      "segundos_no_deterministas": round(cum + tm.get("gas_k_%d" % k, 0.0), 4)})
    out["recursos_para_G4"] = {
        "estimador_por_k": per_k, "estimador_preparacion": est["preparacion"],
        "estimador_pico_memoria_bytes_no_determinista": rec_est.get("pico_memoria_bytes"),
        "referencia_global": ref_cost, "historias_totales": 1 << (N + M),
        "nota": "tiempos y memoria de la referencia en el archivo de recursos del contraste",
    }


if __name__ == "__main__":
    main()
