#!/usr/bin/env python3
"""H6-U4. Estimador con certificado finito por componentes para PATH25 con
ruido fenomenologico mixto: fallo de dato con probabilidad p y lectura
erronea de check con probabilidad q (PLAN_UNIDAD.txt).

Autonomo: no importa codigo de otras unidades, no lee resultados ni tablas
y no calcula ninguna referencia global. Solo biblioteca estandar y
aritmetica racional exacta.

Uso: python -I -B estimador_h6_u4.py --output <nuevo.json> --sello <nuevo.json>
            --recursos <nuevo.json> [--kmax K]
"""
import argparse
import datetime
import hashlib
import json
import os
import sys
import time
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import combinations

T0 = time.perf_counter()
WALL_CAP = 570.0
A_COLS, SEED = 4, (3, 6, 12)
B_ROWS = len(SEED)
N = A_COLS * A_COLS + B_ROWS * B_ROWS
M = A_COLS * B_ROWS
NS = N + M
U0 = 4369
PAIRS = ((F(1, 50), F(1, 100)), (F(1, 100), F(1, 200)), (F(1, 100), F(0)))
PAIR_NAMES = ("principal", "lectura", "control_q0")
TARGET = F(1, 10)
KMAX_DEFAULT = 5
PRED_RUNS = (37, 703, 8473, 74518, 510415)

C = {"scores": 0, "sindromes_distintos": 0, "ejecuciones_locales": 0, "pruebas_conexidad": 0,
     "decoraciones": 0}
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


def over_time():
    return time.perf_counter() - T0 > WALL_CAP


def bits(v):
    while v:
        b = v & -v
        yield b.bit_length() - 1
        v ^= b


def lin(cols, v):
    r = 0
    while v:
        b = v & -v
        r ^= cols[b.bit_length() - 1]
        v ^= b
    return r


def key(pair):
    return "%s|%s" % (pair[0], pair[1])


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


def show(x):
    with localcontext() as ctx:
        ctx.prec = 30
        return {"exacto": "%d/%d" % (x.numerator, x.denominator),
                "decimal": format(Decimal(x.numerator) / Decimal(x.denominator), ".12E")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--sello", required=True)
    ap.add_argument("--recursos", required=True)
    ap.add_argument("--kmax", type=int, default=KMAX_DEFAULT)
    args = ap.parse_args()
    for p in (args.output, args.sello, args.recursos):
        if os.path.exists(p):
            raise SystemExit("La salida ya existe; use un nombre nuevo: " + p)
    kmax = args.kmax
    assert 1 <= kmax <= KMAX_DEFAULT
    out = {"unit": "H6-U4", "programa": "estimador", "kmax": kmax,
           "forma_V3": "certificado finito por componentes; no es el algoritmo general del teorema",
           "lee_referencia_o_resultados_de_otras_unidades": False}
    status = "PASS"
    try:
        run(kmax, out)
    except Stop as ex:
        status = "FAIL"
        out["stopped_at"] = str(ex)
    out["status"] = status
    out["checks"] = CHECKS
    with open(args.output, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(out, ensure_ascii=False, indent=1, sort_keys=True) + "\n")
    digest = hashlib.sha256(open(args.output, "rb").read()).hexdigest()
    now = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=-5))).isoformat()
    with open(args.sello, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps({"archivo": os.path.basename(args.output), "sha256": digest, "sellado": now,
                            "estimador_sha256": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
                            "status": status, "kmax": kmax, "k_completados": out.get("k_completados")},
                           ensure_ascii=False, indent=1) + "\n")
    TIMES["total"] = round(time.perf_counter() - T0, 4)
    with open(args.recursos, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps({"nota": "tiempos de reloj en segundos y memoria; no deterministas; no son resultado",
                            "tiempos": TIMES, "pico_memoria_bytes": peak_memory_bytes()},
                           ensure_ascii=False, indent=1) + "\n")
    print("%s estimador kmax=%d completados=%s controles=%d t=%.2fs sha=%s"
          % (status, kmax, out.get("k_completados"), len(CHECKS), TIMES["total"], digest[:16]))
    if status != "PASS":
        print("DETENIDO EN:", out.get("stopped_at"))
        sys.exit(1)


def run(kmax, out):
    # ------------------------------------------------ paso 1: preparacion
    t = time.perf_counter()
    a, b = A_COLS, B_ROWS
    h = [[(SEED[r] >> c) & 1 for c in range(a)] for r in range(b)]
    hx, hz = [0] * (a * b), [0] * (b * a)
    for i in range(a):
        for j in range(b):
            row = 0
            for k in range(a):
                if h[j][k]:
                    row |= 1 << (i * a + k)
            for l in range(b):
                if h[l][i]:
                    row |= 1 << (a * a + l * b + j)
            hx[i * b + j] = row
    for i in range(b):
        for j in range(a):
            row = 0
            for k in range(a):
                if h[i][k]:
                    row |= 1 << (k * a + j)
            for l in range(b):
                if h[l][j]:
                    row |= 1 << (a * a + i * b + l)
            hz[i * a + j] = row
    cols = [sum(1 << c for c in range(M) if hx[c] >> v & 1) for v in range(N)]
    check("I0 conmutacion de H_X y H_Z", all((x & z).bit_count() % 2 == 0 for x in hx for z in hz))

    def span(vectors):
        piv = {}
        for v in vectors:
            while v:
                k = v.bit_length() - 1
                if k in piv:
                    v ^= piv[k]
                else:
                    piv[k] = v
                    break
        return piv

    def inside(v, piv):
        while v:
            k = v.bit_length() - 1
            if k not in piv:
                return False
            v ^= piv[k]
        return True

    bx, bz = span(hx), span(hz)
    check("I1 k logico = 1 y u_0 protegido: conmuta con H_Z y no esta en el espacio de filas de H_X",
          N - len(bx) - len(bz) == 1 and all((U0 & g).bit_count() % 2 == 0 for g in hz) and not inside(U0, bx))

    cand = set()
    for g in hz:
        f = g
        while f:
            if lin(cols, f):
                cand.add(f)
            f = (f - 1) & g
    ACT = sorted(cand, key=lambda v: (v.bit_count(), v))
    FOOT = [lin(cols, f) for f in ACT]
    WID = [p.bit_count() for p in FOOT]
    INV = [[] for _ in range(M)]
    for j, ft in enumerate(FOOT):
        for c in bits(ft):
            INV[c].append(j)
    psi = cols + [1 << c for c in range(M)]
    adjK = [0] * NS
    for f, phi in zip(ACT, FOOT):
        rm = 0
        for x in range(NS):
            if psi[x] & phi:
                rm |= 1 << x
        for x in bits(rm):
            adjK[x] |= rm | f
        for y in bits(f):
            adjK[y] |= rm
    for x in range(NS):
        adjK[x] &= ~(1 << x)
    prep = {"candidatos": len(ACT), "entradas_listas_inversas": sum(len(x) for x in INV),
            "aristas_host_K": sum(x.bit_count() for x in adjK) // 2, "sitios": NS}
    check("N2a host K con 388 aristas (cifra publicada)", prep["aristas_host_K"] == 388, prep["aristas_host_K"])
    TIMES["preparacion"] = round(time.perf_counter() - t, 4)

    # ------------------------------------------------ decoder al vuelo
    RES = {}

    def decode(s0):
        got = RES.get(s0)
        if got is not None:
            return got
        s = s0
        cor = touched = 0
        while s:
            readers = set()
            for c in bits(s):
                readers.update(INV[c])
            best, num, den = -1, 0, 1
            for j in sorted(readers):
                d = WID[j]
                g = 2 * (s & FOOT[j]).bit_count() - d
                if best < 0 or g * den > num * d:
                    best, num, den = j, g, d
            C["scores"] += len(readers)
            if best < 0 or 2 * num < den:
                break
            cor ^= ACT[best]
            touched |= ACT[best]
            s ^= FOOT[best]
        C["sindromes_distintos"] += 1
        got = (cor, touched)
        RES[s0] = got
        return got

    DATA_MASK = (1 << N) - 1
    ALL = (1 << NS) - 1
    nP = len(PAIRS)
    zp = [p / (1 - p) for p, q in PAIRS]
    zq = [q / (1 - q) for p, q in PAIRS]
    B0 = [(1 - p) ** N * (1 - q) ** M for p, q in PAIRS]

    POLY = {}      # Y -> {(a,b): [n, nT, nTu]}
    table = []
    done_k = 0
    for size in range(1, kmax + 1):
        if over_time():
            break
        # ---------------------------------------- paso 2: catalogo de tamano `size`
        t = time.perf_counter()
        for S in combinations(range(NS), size):
            sm = s = 0
            for i in S:
                sm |= 1 << i
                s ^= psi[i]
            C["ejecuciones_locales"] += 1
            cor, touched = decode(s)
            Y = sm | touched
            C["pruebas_conexidad"] += 1
            seen = 0
            todo = Y & -Y
            while todo:
                bb = todo & -todo
                todo ^= bb
                seen |= bb
                todo |= adjK[bb.bit_length() - 1] & Y & ~seen
            if seen != Y:
                continue
            C["decoraciones"] += 1
            E = sm & DATA_MASK
            r = E ^ cor
            na = E.bit_count()
            slot = POLY.get(Y)
            if slot is None:
                slot = {}
                POLY[Y] = slot
            ent = slot.get((na, size - na))
            if ent is None:
                ent = [0, 0, 0]
                slot[(na, size - na)] = ent
            ent[0] += 1
            if lin(cols, r) == 0:
                ent[1] += 1
                ent[2] += -1 if (U0 & r).bit_count() & 1 else 1
        TIMES["catalogo_tamano_%d" % size] = round(time.perf_counter() - t, 4)
        if over_time():
            break

        # ---------------------------------------- paso 3: gas truncado por recurrencia
        t = time.perf_counter()
        mono = {}

        def monomial(ri, na, nb_):
            kk = (ri, na, nb_)
            v = mono.get(kk)
            if v is None:
                v = zp[ri] ** na * zq[ri] ** nb_
                mono[kk] = v
            return v

        WD = {}
        BYMIN = [[] for _ in range(NS)]
        for Y, slot in POLY.items():
            vals = []
            for ri in range(nP):
                w0 = wt = wu = F(0)
                for (na, nb_), (n0, nT, nTu) in slot.items():
                    mm = monomial(ri, na, nb_)
                    w0 += n0 * mm
                    if nT:
                        wt += nT * mm
                    if nTu:
                        wu += nTu * mm
                vals += [w0, wt, wu]
            nb = Y
            for x in bits(Y):
                nb |= adjK[x]
            ndec = sum(v[0] for v in slot.values())
            rec = (tuple(vals), nb, ndec, Y)
            WD[Y] = rec
            BYMIN[(Y & -Y).bit_length() - 1].append(rec)
        NV = 3 * nP
        ONE = tuple([F(1)] * NV + [1, 1])
        memo = {0: ONE}
        gas = {"terminos": 0, "consultas_por_subconjunto": 0}

        def xi(J):
            got = memo.get(J)
            if got is not None:
                return got
            x = J & -J
            rest = J ^ x
            acc = list(xi(rest))
            lst = BYMIN[x.bit_length() - 1]
            m_ = rest.bit_count()
            if m_ < 22 and (1 << m_) <= len(lst):
                sub = rest
                while True:
                    recd = WD.get(sub | x)
                    gas["consultas_por_subconjunto"] += 1
                    if recd is not None:
                        tt = xi(J & ~recd[1])
                        w = recd[0]
                        for i in range(NV):
                            acc[i] += w[i] * tt[i]
                        acc[NV] += tt[NV]
                        acc[NV + 1] += recd[2] * tt[NV + 1]
                        gas["terminos"] += 1
                    if sub == 0:
                        break
                    sub = (sub - 1) & rest
            else:
                for recd in lst:
                    if recd[3] & ~J == 0:
                        tt = xi(J & ~recd[1])
                        w = recd[0]
                        for i in range(NV):
                            acc[i] += w[i] * tt[i]
                        acc[NV] += tt[NV]
                        acc[NV + 1] += recd[2] * tt[NV + 1]
                        gas["terminos"] += 1
            got = tuple(acc)
            memo[J] = got
            return got

        # la cadena de sufijos se evalua de abajo arriba para no profundizar la pila
        for i in range(NS - 1, -1, -1):
            xi(ALL & ~((1 << i) - 1))
        top = memo[ALL]
        TIMES["gas_k_%d" % size] = round(time.perf_counter() - t, 4)

        # ---------------------------------------- paso 4: certificado
        row = {"k": size, "polimeros": len(POLY), "familias": top[NV], "historias_retenidas": top[NV + 1],
               "historias_totales": 1 << NS, "estados_del_gas": len(memo), "terminos_del_gas": gas["terminos"],
               "consultas_por_subconjunto": gas["consultas_por_subconjunto"], "acumulado": dict(C), "tasas": {}}
        for ri, (p, q) in enumerate(PAIRS):
            X0, XT, XTu = top[3 * ri], top[3 * ri + 1], top[3 * ri + 2]
            ret = B0[ri] * X0
            eps = 1 - ret
            P0k = B0[ri] * XT
            Psk = B0[ri] * (XT - XTu) / 2
            dist = [F(1)]
            for qq in [p] * N + [q] * M:
                nd = [F(0)] * (len(dist) + 1)
                for i_, v_ in enumerate(dist):
                    nd[i_] += v_ * (1 - qq)
                    nd[i_ + 1] += v_ * qq
                dist = nd
            tail = sum(dist[size + 1:], F(0))
            check("I2 k=%d %s: 0 <= eps_k <= Pr[mas de k fallos]" % (size, key((p, q))), 0 <= eps <= tail)
            ivals = {"P_0": (P0k, P0k + eps), "P_sil": (Psk, Psk + eps),
                     "Q_sil": (Psk / (P0k + eps), (Psk + eps) / (P0k + eps))}
            cell = {"par": PAIR_NAMES[ri], "B_0": show(B0[ri]), "masa_retenida": show(ret), "eps_k": show(eps),
                    "cola_binomial": show(tail),
                    "bits_max": max(z.numerator.bit_length() + z.denominator.bit_length() for z in (X0, XT, XTu))}
            for nm, (lo, hi) in ivals.items():
                cc = {"L": show(lo), "U": show(hi)}
                if lo > 0:
                    cc["valor_informado"] = show(2 * lo * hi / (lo + hi))
                    cc["error_relativo_garantizado"] = show((hi - lo) / (hi + lo))
                    cc["anchura_relativa"] = show((hi - lo) / lo)
                else:
                    cc["valor_informado"] = None
                    cc["error_relativo_garantizado"] = None
                    cc["anchura_relativa"] = None
                cell[nm] = cc
            cell["_Psil_L"] = Psk
            cell["_eps"] = eps
            row["tasas"][key((p, q))] = cell
        table.append(row)
        done_k = size
        del memo, WD, BYMIN

    out["k_completados"] = done_k
    out["truncado_por_techo_de_tiempo"] = done_k < kmax
    if done_k == 0:
        raise Stop("ningun k completado dentro del techo de tiempo")

    # ------------------------------------------------ controles internos y parada
    runs = [row["acumulado"]["ejecuciones_locales"] for row in table]
    check("N1 historias locales acumuladas por k iguales a la prediccion del plan",
          runs == list(PRED_RUNS[:done_k]), runs)
    if done_k >= 3:
        check("N2b hasta k=3 hay 4705 historias con soporte conexo (cifra publicada)",
              table[2]["acumulado"]["decoraciones"] == 4705, table[2]["acumulado"]["decoraciones"])
    for pair in PAIRS:
        kk = key(pair)
        Ls = [row["tasas"][kk]["_Psil_L"] for row in table]
        eps_ = [row["tasas"][kk]["_eps"] for row in table]
        check("I3 %s: L_k no decrece y eps_k no crece" % kk,
              all(Ls[i] <= Ls[i + 1] for i in range(len(Ls) - 1))
              and all(eps_[i] >= eps_[i + 1] for i in range(len(eps_) - 1)))
        check("N3 %s: P_sil^(1) = 0" % kk, Ls[0] == 0)
    kstar = {}
    for pair in PAIRS:
        kk = key(pair)
        ks = None
        for row in table:
            cell = row["tasas"][kk]
            lo = cell["_Psil_L"]
            if lo > 0 and cell["_eps"] / (2 * lo + cell["_eps"]) <= TARGET:
                ks = row["k"]
                break
        kstar[kk] = ks
    for row in table:
        for cell in row["tasas"].values():
            del cell["_Psil_L"]
            del cell["_eps"]

    out["instancia"] = {"codigo": "PATH25", "n": N, "checks": M, "sitios": NS, "historias": 1 << NS,
                        "pares_de_tasas": [{"nombre": nm, "p": str(p), "q": str(q)}
                                           for nm, (p, q) in zip(PAIR_NAMES, PAIRS)], "u_0": U0}
    out["preparacion"] = prep
    out["objetivo"] = "error relativo garantizado de P_sil <= 1/10"
    out["k_estrella"] = kstar
    out["tabla"] = table
    out["estimacion_informada"] = {
        kk: ({"k": kstar[kk], "certificado": table[kstar[kk] - 1]["tasas"][kk],
              "recursos_acumulados": table[kstar[kk] - 1]["acumulado"],
              "polimeros": table[kstar[kk] - 1]["polimeros"], "familias": table[kstar[kk] - 1]["familias"],
              "estados_del_gas": table[kstar[kk] - 1]["estados_del_gas"],
              "terminos_del_gas": table[kstar[kk] - 1]["terminos_del_gas"]}
             if kstar[kk] else None)
        for kk in kstar
    }


if __name__ == "__main__":
    main()
