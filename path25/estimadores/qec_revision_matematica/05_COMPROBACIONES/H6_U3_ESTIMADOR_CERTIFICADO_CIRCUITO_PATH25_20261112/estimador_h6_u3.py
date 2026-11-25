#!/usr/bin/env python3
"""H6-U3. Estimador con certificado finito por componentes para el circuito
de una ronda sobre PATH25 con ocho ubicaciones ruidosas (PLAN_UNIDAD.txt).

Autonomo: no importa codigo de la unidad 2, no lee sus resultados ni la
tabla congelada del piloto, y no calcula ninguna referencia global. Solo
biblioteca estandar y aritmetica racional exacta.

Uso: python -I -B estimador_h6_u3.py --output <nuevo.json> --sello <nuevo.json>
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
from itertools import combinations, product

T0 = time.perf_counter()
WALL_CAP = 290.0
A_COLS, SEED = 4, (3, 6, 12)
B_ROWS = len(SEED)
N = A_COLS * A_COLS + B_ROWS * B_ROWS
M = A_COLS * B_ROWS
U0 = 4369
RATES = (F(1, 100), F(1, 1000))
TARGET = F(1, 10)
LOCS = (
    ("dwait", 1, 0), ("dwait", 2, 0), ("dwait", 3, 0), ("dwait", 16, 0),
    ("dwait", 9, 0), ("dwait", 10, 0), ("cnot", 10, 8), ("read", 1),
)
NAMES = ("W1", "W2", "W3", "W16", "W9", "W10", "C10_8", "R1")
NL = len(LOCS)
PRED_RUNS = (18, 158, 774, 2454, 5366, 8502, 10422, 10934)

C = {"scores": 0, "sindromes_distintos": 0, "ejecuciones_locales": 0, "pasos_decoder": 0,
     "pruebas_conexidad": 0}
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


def guard():
    if time.perf_counter() - T0 > WALL_CAP:
        raise Stop("techo de tiempo de reloj alcanzado")


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
    ap.add_argument("--kmax", type=int, default=NL)
    args = ap.parse_args()
    for p in (args.output, args.sello, args.recursos):
        if os.path.exists(p):
            raise SystemExit("La salida ya existe; use un nombre nuevo: " + p)
    kmax = args.kmax
    assert 1 <= kmax <= NL
    out = {"unit": "H6-U3", "programa": "estimador", "kmax": kmax,
           "forma_V3": "certificado finito por componentes; no es el algoritmo general del teorema",
           "lee_referencia_o_resultados_de_la_unidad_2": False}
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
                            "status": status, "kmax": kmax}, ensure_ascii=False, indent=1) + "\n")
    TIMES["total"] = round(time.perf_counter() - T0, 4)
    with open(args.recursos, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps({"nota": "tiempos de reloj en segundos y memoria; no deterministas; no son resultado",
                            "tiempos": TIMES, "pico_memoria_bytes": peak_memory_bytes()},
                           ensure_ascii=False, indent=1) + "\n")
    print("%s estimador kmax=%d controles=%d t=%.2fs sha=%s" % (status, kmax, len(CHECKS), TIMES["total"], digest[:16]))
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
    k_log = N - len(bx) - len(bz)
    check("I1 k logico = 1 y u_0 protegido: conmuta con H_Z y no esta en el espacio de filas de H_X",
          k_log == 1 and all((U0 & g).bit_count() % 2 == 0 for g in hz) and not inside(U0, bx))

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
    prep = {"candidatos": len(ACT), "entradas_listas_inversas": sum(len(x) for x in INV)}

    psi = cols + [1 << c for c in range(M)]
    adjK = [0] * (N + M)
    for f, phi in zip(ACT, FOOT):
        rm = 0
        for x in range(N + M):
            if psi[x] & phi:
                rm |= 1 << x
        for x in bits(rm):
            adjK[x] |= rm | f
        for y in bits(f):
            adjK[y] |= rm
    for x in range(N + M):
        adjK[x] &= ~(1 << x)
    prep["aristas_host_K"] = sum(x.bit_count() for x in adjK) // 2

    used_c = [set() for _ in range(M)]
    used_v = [set() for _ in range(N)]
    layer = {}
    for c in range(M):
        for v in bits(hx[c]):
            tt = 1
            while tt in used_c[c] or tt in used_v[v]:
                tt += 1
            layer[(c, v)] = tt
            used_c[c].add(tt)
            used_v[v].add(tt)
    L = max(x.bit_count() for x in hx) + max(x.bit_count() for x in cols) - 1
    prep["capas"] = L
    prep["ubicaciones"] = (L + 2) * (N + M) - sum(x.bit_count() for x in hx)

    def Tmask(v, tt):
        return sum(1 << c for c in bits(cols[v]) if layer[(c, v)] > tt)

    def x_of(p):
        return 1 if p in (1, 2) else 0

    # clases de respuesta por ubicacion: estado 0 = categoria base
    NST, CNT, ER, SR, base, kinds = [], [], [], [], [], []
    CW = [[] for _ in RATES]
    BC = [F(1) for _ in RATES]

    def law(kind, cat, e):
        if kind == "dwait":
            return 1 - e if cat == 0 else e / 3
        if kind == "cnot":
            return 1 - e if cat == 0 else e / 15
        return (1 - e if cat & 3 == 0 else e / 3) * (e if cat >> 2 else 1 - e)

    for loc in LOCS:
        kind = loc[0]
        kinds.append(kind)
        if kind == "dwait":
            v, tt = loc[1], loc[2]
            resp = [(0, 0), (0, 0), (1 << v, Tmask(v, tt))]
            st_of = [0, 2, 2, 1]
            base.append(v)
        elif kind == "cnot":
            v, c = loc[1], loc[2]
            tt = layer[(c, v)]
            resp = [(0, 0)]
            for b_ in (0, 1):
                for a_ in (0, 1):
                    resp.append(((1 << v) if a_ else 0, (Tmask(v, tt) if a_ else 0) ^ ((1 << c) if b_ else 0)))
            st_of = [0 if q == 0 else 1 + x_of(q >> 2) + 2 * x_of(q & 3) for q in range(16)]
            base.append(v)
        else:
            c = loc[1]
            resp = [(0, 0), (0, 0), (0, 1 << c)]
            st_of = [0 if r == 0 else 1 + (x_of(r & 3) ^ (r >> 2)) for r in range(8)]
            base.append(N + c)
        ns = len(resp)
        NST.append(ns)
        ER.append([x[0] for x in resp])
        SR.append([x[1] for x in resp])
        cnt = [0] * ns
        for cat in range(len(st_of)):
            cnt[st_of[cat]] += 1
        CNT.append(cnt)
        for ri, e in enumerate(RATES):
            ww = [F(0)] * ns
            for cat in range(len(st_of)):
                ww[st_of[cat]] += law(kind, cat, e)
            CW[ri].append([x / ww[0] for x in ww])
            BC[ri] *= law(kind, 0, e)
    total_assign = 1
    for cnt in CNT:
        total_assign *= sum(cnt)

    NV = NL + N
    vbase = base + list(range(N))
    ADJH = [0] * NV
    for x in range(NV):
        for y in range(NV):
            if x != y and (vbase[x] == vbase[y] or adjK[vbase[x]] >> vbase[y] & 1):
                ADJH[x] |= 1 << y
    prep["vertices_host_primitivo_usados"] = NV
    TIMES["preparacion"] = round(time.perf_counter() - t, 4)

    # ------------------------------------------------ decoder al vuelo
    TRANS = {}

    def step(s):
        got = TRANS.get(s)
        if got is not None:
            return got
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
        C["sindromes_distintos"] += 1
        got = (-1, s) if (best < 0 or 2 * num < den) else (best, s ^ FOOT[best])
        TRANS[s] = got
        return got

    def decode(s):
        C["ejecuciones_locales"] += 1
        cor = touched = 0
        steps = 0
        while True:
            j, nxt = step(s)
            if j < 0:
                break
            cor ^= ACT[j]
            touched |= ACT[j]
            s = nxt
            steps += 1
        C["pasos_decoder"] += steps
        return cor, touched

    def connected(Y):
        C["pruebas_conexidad"] += 1
        seen = 0
        todo = Y & -Y
        while todo:
            bb = todo & -todo
            todo ^= bb
            seen |= bb
            todo |= ADJH[bb.bit_length() - 1] & Y & ~seen
        return seen == Y

    # ------------------------------------------------ pasos 2-4 por tamano
    POLY = {}                 # Y -> [fuentes, N, [[W0, WT, WTu] por tasa]]
    table = []
    snapshots = []
    nR = len(RATES)
    for size in range(1, kmax + 1):
        t = time.perf_counter()
        for S in combinations(range(NL), size):
            src = 0
            for k in S:
                src |= 1 << k
            for st in product(*[range(1, NST[k]) for k in S]):
                mult = 1
                E = s = 0
                for k, sk in zip(S, st):
                    mult *= CNT[k][sk]
                    E ^= ER[k][sk]
                    s ^= SR[k][sk]
                cor, touched = decode(s)
                Y = src | (touched << NL)
                if not connected(Y):
                    continue
                r = E ^ cor
                zero = lin(cols, r) == 0
                mark = -1 if (U0 & r).bit_count() & 1 else 1
                slot = POLY.get(Y)
                if slot is None:
                    slot = [size, 0, [[F(0), F(0), F(0)] for _ in range(nR)]]
                    POLY[Y] = slot
                slot[1] += mult
                for ri in range(nR):
                    w = F(1)
                    for k, sk in zip(S, st):
                        w *= CW[ri][k][sk]
                    slot[2][ri][0] += w
                    if zero:
                        slot[2][ri][1] += w
                        slot[2][ri][2] += mark * w
            guard()
        TIMES["catalogo_tamano_%d" % size] = round(time.perf_counter() - t, 4)

        # paso 3: gas truncado con los polimeros de hasta `size` fuentes
        t = time.perf_counter()
        P = sorted(POLY)
        NB = []
        for Y in P:
            m_ = Y
            for x in bits(Y):
                m_ |= ADJH[x]
            NB.append(m_)
        PW = [POLY[Y] for Y in P]
        acc = {"fam": 0, "N": 0, "compat": 0, "maxfam": 0, "X": [[F(0), F(0), F(0)] for _ in range(nR)]}

        def rec(start, blocked, cnt, ws, depth):
            acc["fam"] += 1
            acc["N"] += cnt
            if depth > acc["maxfam"]:
                acc["maxfam"] = depth
            for ri in range(nR):
                for q in range(3):
                    acc["X"][ri][q] += ws[ri][q]
            acc["compat"] += len(P) - start
            for i in range(start, len(P)):
                if P[i] & blocked == 0:
                    wi = PW[i][2]
                    rec(i + 1, blocked | NB[i], cnt * PW[i][1],
                        [[ws[ri][q] * wi[ri][q] for q in range(3)] for ri in range(nR)], depth + 1)

        rec(0, 0, 1, [[F(1), F(1), F(1)] for _ in range(nR)], 0)
        TIMES["gas_k_%d" % size] = round(time.perf_counter() - t, 4)

        # paso 4: certificado
        row = {"k": size, "polimeros": len(P), "familias": acc["fam"], "tamano_maximo_de_familia": acc["maxfam"],
               "pruebas_compatibilidad": acc["compat"], "asignaciones_retenidas": acc["N"],
               "asignaciones_totales": total_assign, "acumulado": dict(C), "tasas": {}}
        for ri, e in enumerate(RATES):
            X0, XT, XTu = acc["X"][ri]
            ret = BC[ri] * X0
            eps = 1 - ret
            P0k = BC[ri] * XT
            Psk = BC[ri] * (XT - XTu) / 2
            q = []
            for loc in LOCS:
                q.append(1 - law(loc[0], 0, e))
            dist = [F(1)]
            for qq in q:
                nd = [F(0)] * (len(dist) + 1)
                for i_, v_ in enumerate(dist):
                    nd[i_] += v_ * (1 - qq)
                    nd[i_ + 1] += v_ * qq
                dist = nd
            tail = sum(dist[size + 1:], F(0))
            check("I2 k=%d e=%s: 0 <= eps_k <= Pr[mas de k ubicaciones no base]" % (size, e), 0 <= eps <= tail)
            ivals = {"P_0": (P0k, P0k + eps), "P_sil": (Psk, Psk + eps),
                     "Q_sil": (Psk / (P0k + eps), (Psk + eps) / (P0k + eps))}
            cell = {"B_c": show(BC[ri]), "masa_retenida": show(ret), "eps_k": show(eps), "cola_binomial": show(tail),
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
            row["tasas"][str(e)] = cell
        table.append(row)
        snapshots.append(dict(C))

    # ------------------------------------------------ controles internos y parada
    check("N1 ejecuciones locales acumuladas por k iguales a la prediccion del plan",
          [s_["ejecuciones_locales"] for s_ in snapshots] == list(PRED_RUNS[:kmax]),
          [s_["ejecuciones_locales"] for s_ in snapshots])
    for ri, e in enumerate(RATES):
        Ls = [row["tasas"][str(e)]["_Psil_L"] for row in table]
        eps_ = [row["tasas"][str(e)]["_eps"] for row in table]
        check("I3 e=%s: L_k no decrece y eps_k no crece" % e,
              all(Ls[i] <= Ls[i + 1] for i in range(len(Ls) - 1))
              and all(eps_[i] >= eps_[i + 1] for i in range(len(eps_) - 1)))
        check("N3 e=%s: P_sil^(1) = 0" % e, Ls[0] == 0)
        if kmax == NL:
            check("N2a e=%s: eps_8 = 0 exacto" % e, eps_[-1] == 0)
    if kmax == NL:
        check("N2b en k=8: 2887 polimeros y 4374 familias (cifras previstas)",
              table[-1]["polimeros"] == 2887 and table[-1]["familias"] == 4374,
              {"polimeros": table[-1]["polimeros"], "familias": table[-1]["familias"]})
        check("I4 en k=8 se retienen todas las asignaciones", table[-1]["asignaciones_retenidas"] == total_assign)
    kstar = {}
    for e in RATES:
        ks = None
        for row in table:
            cell = row["tasas"][str(e)]
            lo = cell["_Psil_L"]
            if lo > 0 and cell["_eps"] / (2 * lo + cell["_eps"]) <= TARGET:
                ks = row["k"]
                break
        kstar[str(e)] = ks
    for row in table:
        for cell in row["tasas"].values():
            del cell["_Psil_L"]
            del cell["_eps"]

    out["instancia"] = {"codigo": "PATH25", "n": N, "checks": M, "ubicaciones_ruidosas": list(NAMES),
                        "asignaciones": total_assign, "tasas": [str(e) for e in RATES], "u_0": U0}
    out["preparacion"] = prep
    out["objetivo"] = "error relativo garantizado de P_sil <= 1/10"
    out["k_estrella"] = kstar
    out["tabla"] = table
    out["estimacion_informada"] = {
        e: ({"k": kstar[e], "certificado": table[kstar[e] - 1]["tasas"][e],
             "recursos_acumulados": table[kstar[e] - 1]["acumulado"],
             "polimeros": table[kstar[e] - 1]["polimeros"], "familias": table[kstar[e] - 1]["familias"]}
            if kstar[e] else None)
        for e in kstar
    }


if __name__ == "__main__":
    main()
