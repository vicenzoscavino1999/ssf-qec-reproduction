#!/usr/bin/env python3
"""H6-U3. Contraste de la estimacion sellada con la referencia exacta.

Solo arranca si SELLO_ESTIMACION.json existe y coincide con la estimacion.
Calcula la referencia por la via global de la unidad 2 (propagacion Pauli
literal de las 524288 asignaciones, decoder por tabla, clase por rango de
H_Z), usando las funciones del verificador sellado de la unidad 2, la
compara con los recuentos guardados en RESULTADO_H6_U2.json y SOLO DESPUES
abre la estimacion.

Uso: python -I -B contrastar_h6_u3.py --estimacion <json> --sello <json>
            --recursos-estimador <json> --output <nuevo.json> --recursos <nuevo.json>
"""
import argparse
import hashlib
import importlib.util
import json
import os
import sys
import time
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import product

T0 = time.perf_counter()
HERE = os.path.dirname(os.path.abspath(__file__))
U2 = os.path.join(os.path.dirname(HERE), "H6_U2_CIRCUITO_PATH25_MULTICOMPONENTE_20261112")
PIN = {
    "verificar_h6_u2.py": "66ab5ab9e275c2a3daca1246b128a531ed9d89ee7d40fdeb274210d085da28d0",
    "RESULTADO_H6_U2.json": "e23905fb361a5a5dea8cc5414461182b1faf909b33d3b0979215bf9ff1347250",
}
RATES = (F(1, 100), F(1, 1000))
N, M = 25, 12
LOCS = (
    ("dwait", 1, 0), ("dwait", 2, 0), ("dwait", 3, 0), ("dwait", 16, 0),
    ("dwait", 9, 0), ("dwait", 10, 0), ("cnot", 10, 8), ("read", 1),
)
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
    out = {"unit": "H6-U3", "programa": "contraste"}
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
          seal["sha256"] == est_sha and seal["status"] == "PASS" and seal["kmax"] == 8
          and seal["estimador_sha256"] == sha_file(os.path.join(HERE, "estimador_h6_u3.py")),
          {"sellado": seal["sellado"], "sha256": est_sha})
    check("S2 verificador y resultado de la unidad 2 sin cambios (hashes fijados)",
          all(sha_file(os.path.join(U2, k)) == v for k, v in PIN.items()))
    out["sello"] = seal

    # ------------------------------------------------ 1. referencia exacta global
    t = time.perf_counter()
    spec = importlib.util.spec_from_file_location("u2_sellado", os.path.join(U2, "verificar_h6_u2.py"))
    u2 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(u2)
    hx, hz, cols, actions, foot = u2.build_masks()
    T1 = u2.transitions_bits(foot)
    COR1, END1, TR1, TOUCH1 = u2.resolve(T1, actions)
    bz = u2.rank_basis(hz)
    layer = u2.greedy_schedule(hx)
    L = max(bin(x).count("1") for x in hx) + max(bin(x).count("1") for x in cols) - 1
    TIMES["referencia_preparacion"] = round(time.perf_counter() - t, 4)

    t = time.perf_counter()
    NL = len(LOCS)
    alph = [4 if l[0] == "dwait" else 16 if l[0] == "cnot" else 8 for l in LOCS]
    prog = []
    for k, l in enumerate(LOCS):
        if l[0] == "dwait" and l[2] == 0:
            prog.append((1, k, 1 << l[1]))
    for tt in range(1, L + 1):
        prog.append((0, tuple((1 << v, 1 << c) for (c, v), lt in sorted(layer.items()) if lt == tt)))
        for k, l in enumerate(LOCS):
            if l[0] == "cnot" and layer[(l[2], l[1])] == tt:
                prog.append((2, k, 1 << l[1], 1 << l[2]))
    reads = [(k, 1 << l[1]) for k, l in enumerate(LOCS) if l[0] == "read"]
    n_wait = sum(1 for l in LOCS if l[0] == "dwait")
    sig_of = []
    for l in LOCS:
        if l[0] == "dwait":
            sig_of.append([0, 1, 1, 1])
        elif l[0] == "cnot":
            sig_of.append([0] + [n_wait + 1] * 15)
        else:
            sig_of.append([((n_wait + 1) * 2) * ((1 if r & 3 else 0) + 2 * (r >> 2)) for r in range(8)])
    NSIG = (n_wait + 1) * 8
    c_tot, c_zero, c_sil = [0] * NSIG, [0] * NSIG, [0] * NSIG
    memo = {}
    done = 0
    rng = range(NL)
    for cats in product(*[range(x) for x in alph]):
        xd = zd = xa = za = 0
        for st in prog:
            op = st[0]
            if op == 0:
                for vb, cb in st[1]:
                    if xd & vb:
                        xa ^= cb
                    if za & cb:
                        zd ^= vb
            elif op == 1:
                p = cats[st[1]]
                if p:
                    if p != 3:
                        xd ^= st[2]
                    if p != 1:
                        zd ^= st[2]
            else:
                q = cats[st[1]]
                if q:
                    pc, pt = q >> 2, q & 3
                    if pc == 1 or pc == 2:
                        xd ^= st[2]
                    if pc == 2 or pc == 3:
                        zd ^= st[2]
                    if pt == 1 or pt == 2:
                        xa ^= st[3]
                    if pt == 2 or pt == 3:
                        za ^= st[3]
        s = xa
        for k, cb in reads:
            r = cats[k]
            if r:
                p = r & 3
                if p == 1 or p == 2:
                    s ^= cb
                if r >> 2:
                    s ^= cb
        R = xd ^ COR1[s]
        fl = memo.get(R)
        if fl is None:
            z0 = u2.lin(cols, R) == 0
            fl = (z0, z0 and not u2.in_span(R, bz))
            memo[R] = fl
        sg = 0
        for k in rng:
            sg += sig_of[k][cats[k]]
        c_tot[sg] += 1
        if fl[0]:
            c_zero[sg] += 1
            if fl[1]:
                c_sil[sg] += 1
        done += 1
    TIMES["referencia_enumeracion"] = round(time.perf_counter() - t, 4)
    u2res = json.loads(open(os.path.join(U2, "RESULTADO_H6_U2.json"), encoding="utf-8").read())
    rc = u2res["recuentos_por_tipo"]
    check("R1 la referencia recorre las 524288 asignaciones y sus recuentos por tipo son los guardados en la unidad 2",
          done == 524288 and c_tot == rc["total"] and c_zero == rc["filtro_fisico_nulo"]
          and c_sil == rc["fallo_silencioso"])

    def sig_prob(sig, e):
        nw = sig % (n_wait + 1)
        cn = (sig // (n_wait + 1)) % 2
        rp = (sig // ((n_wait + 1) * 2)) % 2
        rb = (sig // ((n_wait + 1) * 4)) % 2
        return ((e / 3) ** nw * (1 - e) ** (n_wait - nw) * (e / 15 if cn else 1 - e)
                * (e / 3 if rp else 1 - e) * (e if rb else 1 - e))

    exact = {}
    for e in RATES:
        sp = [sig_prob(g, e) for g in range(NSIG)]
        tot = sum(c_tot[g] * sp[g] for g in range(NSIG))
        P0 = sum(c_zero[g] * sp[g] for g in range(NSIG))
        Ps = sum(c_sil[g] * sp[g] for g in range(NSIG))
        check("R2 e=%s: las probabilidades de la referencia suman 1" % e, tot == 1)
        exact[str(e)] = {"P_0": P0, "P_sil": Ps, "Q_sil": Ps / P0}
    out["referencia_exacta"] = {e: {k: show(v) for k, v in d.items()} for e, d in exact.items()}
    ref_cost = {"asignaciones_propagadas": done, "puertas_CNOT_aplicadas": done * 42,
                "scores_tabla_completa": (1 << M) * len(actions), "sindromes_en_tabla": 1 << M,
                "residuos_distintos": len(memo)}

    # ------------------------------------------------ 2. AHORA se abre la estimacion
    est = json.loads(open(args.estimacion, encoding="utf-8").read())
    check("S3 la estimacion declara que no leyo referencia ni resultados de la unidad 2",
          est["lee_referencia_o_resultados_de_la_unidad_2"] is False and est["status"] == "PASS")
    table = est["tabla"]
    rows = []
    bad_contain, bad_err = [], []
    for row in table:
        k = row["k"]
        for e in RATES:
            cell = row["tasas"][str(e)]
            for nm in ("P_0", "P_sil", "Q_sil"):
                lo, hi = fr(cell[nm]["L"]), fr(cell[nm]["U"])
                ex = exact[str(e)][nm]
                ok = lo <= ex <= hi
                if not ok:
                    bad_contain.append([k, str(e), nm])
                item = {"k": k, "e": str(e), "cantidad": nm, "L": cell[nm]["L"]["decimal"],
                        "U": cell[nm]["U"]["decimal"], "exacto": show(ex)["decimal"], "contenido": ok}
                if cell[nm]["valor_informado"] is not None:
                    val = fr(cell[nm]["valor_informado"])
                    gar = fr(cell[nm]["error_relativo_garantizado"])
                    real = abs(val - ex) / ex
                    item["error_garantizado"] = show(gar)["decimal"]
                    item["error_real"] = show(real)["decimal"]
                    if real > gar:
                        bad_err.append([k, str(e), nm])
                rows.append(item)
    check("E1 los valores exactos de P_0, P_sil y Q_sil estan dentro de sus intervalos para todo k=1..8 y las dos "
          "tasas", not bad_contain, bad_contain)
    last = table[-1]
    check("E2 en k=8, eps_8=0 y los tres valores coinciden exactamente con la referencia",
          all(fr(last["tasas"][str(e)]["eps_k"]) == 0
              and all(fr(last["tasas"][str(e)][nm]["L"]) == exact[str(e)][nm] == fr(last["tasas"][str(e)][nm]["U"])
                      for nm in ("P_0", "P_sil", "Q_sil")) for e in RATES))
    check("E3 el error real del valor informado no supera el error garantizado (todo k con L>0, las tres cantidades)",
          not bad_err, bad_err)
    mono = True
    for e in RATES:
        Ls = [fr(r["tasas"][str(e)]["P_sil"]["L"]) for r in table]
        ep = [fr(r["tasas"][str(e)]["eps_k"]) for r in table]
        mono = mono and all(Ls[i] <= Ls[i + 1] for i in range(7)) and all(ep[i] >= ep[i + 1] for i in range(7))
    check("E4 L_k no decrece y eps_k no crece con k", mono)
    ks = est["k_estrella"]
    inf = est["estimacion_informada"]
    check("E5 la estimacion informada (k*) existe en las dos tasas y cumple el objetivo: error garantizado de P_sil "
          "<= 10 %",
          all(ks[str(e)] is not None
              and fr(inf[str(e)]["certificado"]["P_sil"]["error_relativo_garantizado"]) <= F(1, 10) for e in RATES),
          ks)
    check("N2 en k=8 el estimador encuentra los mismos polimeros y familias que la unidad 2",
          last["polimeros"] == u2res["counters"]["polimeros"]
          and last["familias"] == u2res["counters"]["familias_compatibles"],
          {"polimeros": last["polimeros"], "familias": last["familias"]})

    # ------------------------------------------------ 3. tablas de salida
    out["contraste_por_k"] = rows
    informed = {}
    for e in RATES:
        cert = inf[str(e)]["certificado"]
        informed[str(e)] = {"k": ks[str(e)], "polimeros": inf[str(e)]["polimeros"],
                            "familias": inf[str(e)]["familias"],
                            "recursos_acumulados": inf[str(e)]["recursos_acumulados"],
                            "eps_k": cert["eps_k"]["decimal"]}
        for nm in ("P_0", "P_sil", "Q_sil"):
            val, gar = fr(cert[nm]["valor_informado"]), fr(cert[nm]["error_relativo_garantizado"])
            ex = exact[str(e)][nm]
            informed[str(e)][nm] = {"L": cert[nm]["L"]["decimal"], "U": cert[nm]["U"]["decimal"],
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
        per_k.append({"k": k, "ejecuciones_locales": row["acumulado"]["ejecuciones_locales"],
                      "sindromes_distintos": row["acumulado"]["sindromes_distintos"],
                      "scores": row["acumulado"]["scores"], "pasos_decoder": row["acumulado"]["pasos_decoder"],
                      "pruebas_conexidad": row["acumulado"]["pruebas_conexidad"], "polimeros": row["polimeros"],
                      "familias": row["familias"], "pruebas_compatibilidad": row["pruebas_compatibilidad"],
                      "asignaciones_retenidas": row["asignaciones_retenidas"],
                      "bits_max": max(c["bits_max"] for c in row["tasas"].values()),
                      "segundos_no_deterministas": round(cum + tm.get("gas_k_%d" % k, 0.0), 4)})
    out["recursos_para_G4"] = {
        "estimador_por_k": per_k,
        "estimador_preparacion": est["preparacion"],
        "estimador_pico_memoria_bytes_no_determinista": rec_est.get("pico_memoria_bytes"),
        "referencia_global": ref_cost,
        "referencia_por_clases_equivalente": "el estimador con k=8 ejecuta 10934 historias por clase: es la "
                                             "enumeracion exacta agrupada por clases de respuesta",
        "nota": "tiempos y memoria de la referencia en el archivo de recursos del contraste",
    }


if __name__ == "__main__":
    main()
