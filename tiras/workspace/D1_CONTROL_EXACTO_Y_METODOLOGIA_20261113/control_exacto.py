#!/usr/bin/env python3
"""D1-CE: control exacto de eps_k (k <= 3) en las celdas de D1, segun
PLAN_CONTROL.txt (sellado en DEC-036) y METODOLOGIA_FIJACION_PARAMETROS.txt.

Para cada instancia y cada celda (p, q) calcula en aritmetica exacta:
- certificado por componentes (el de U4): eps_comp(k) = 1 - B_0 Xi_0^(k) y la
  cota inferior L_comp(k) de P_sil, con el gas truncado por la recurrencia de
  borrado del vertice menor;
- truncacion ingenua: eps_ing(k) = Pr[mas de k fallos] y L_ing(k), masa de las
  historias silenciosas con a lo sumo k fallos.
Enteros escalados: con p = 1/(D+1), z_p = 1/D; el estado J se multiplica por
el producto de los D de sus sitios y la masa retenida es Xi~/prod(D+1).

Construccion de instancias copiada de d1_viabilidad.py (D1). Solo biblioteca
estandar. Lee FASE_A.json, FASE_B.json y FASE_C.json de D1 y
ESTIMACION_H6_U4.json de U4 solo para validar y para el sesgo de la
prediccion, con su SHA-256 comprobado.

Uso (desde D:\\paper ia cuantica\\qec_revision_matematica\\05_COMPROBACIONES):
  python -I -B <dir>/control_exacto.py --output <nuevo.json> --recursos <nuevo.json> [--solo PATH25]
"""
import argparse
import hashlib
import json
import math
import os
import random
import sys
import time
from array import array
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import combinations

T0 = time.perf_counter()
LLAMADAS_MAX = 1_200_000
SEG_MAX = 7200.0
SEMILLA = 20261112
KMAX = 3
D1 = "D1_VIABILIDAD_ALCANCE_CERTIFICADO_20261112/"
U4 = "H6_U4_ESTIMADOR_CERTIFICADO_FENOMENOLOGICO_PATH25_20261112/"
SHA_FASE_A = "e68435581f661be720cc03be80aae62261b8810fdc0aa09ea5b3af934b87de2c"
SHA_FASE_B = "31c89e888ee86989ca772f87352ec851cc9de2e15dc40b71510159f2bd6b711b"
SHA_U4 = "654367777c23c83cdc6cccbbdb756544f85cca57ad0e86a680d6deec473304e9"
GRID_P = [F(1, 10 ** 3), F(1, 10 ** 4), F(1, 10 ** 5), F(1, 10 ** 6)]
U4_PAIRS = [(F(1, 50), F(1, 100)), (F(1, 100), F(1, 200)), (F(1, 100), F(0))]
TARGETS = (F(1, 20), F(1, 10), F(1, 4))
PROPOSED = ("S1_L6", F(1, 10 ** 4), F(1, 2 * 10 ** 4))

CHECKS = []
TIMES = {}
CALLS = [0]


class Stop(Exception):
    pass


def check(name, cond, detail=None):
    row = {"name": name, "pass": bool(cond)}
    if detail is not None:
        row["detail"] = detail
    CHECKS.append(row)
    if not cond:
        raise Stop(name)


def sha_file(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


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


def dec(x):
    with localcontext() as ctx:
        ctx.prec = 30
        return float(Decimal(x.numerator) / Decimal(x.denominator))


# ------------------------------------------------------------------ copiado de d1_viabilidad.py
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


def classical_distance(rows, ncols):
    best = None
    for x in range(1, 1 << ncols):
        if all((r & x).bit_count() % 2 == 0 for r in rows):
            w = x.bit_count()
            if best is None or w < best:
                best = w
    return best


def path_seed(L):
    return [3 << r for r in range(L - 1)], L


def regular_seed_34():
    rng = random.Random(SEMILLA)
    ncols, nrows = 8, 6
    valid = []
    while len(valid) < 2000:
        sockets = [c for c in range(ncols) for _ in range(3)]
        rng.shuffle(sockets)
        rows, ok = [], True
        for r in range(nrows):
            cs = sockets[4 * r:4 * r + 4]
            if len(set(cs)) != 4:
                ok = False
                break
            rows.append(sum(1 << c for c in cs))
        if not ok or len(span(rows)) != nrows:
            continue
        valid.append((classical_distance(rows, ncols), len(valid), rows))
    best = max(valid, key=lambda t: (t[0], -t[1]))
    return best[2], ncols


def build(name, rows, a):
    b = len(rows)
    h = [[(rows[r] >> c) & 1 for c in range(a)] for r in range(b)]
    N, M = a * a + b * b, a * b
    NS = N + M
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
    assert all((x & z).bit_count() % 2 == 0 for x in hx for z in hz), "H_X y H_Z no conmutan"
    bx, bz = span(hx), span(hz)
    klog = N - len(bx) - len(bz)
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
    coord = [None] * NS
    for i in range(a):
        for k in range(a):
            coord[i * a + k] = (2 * i, 2 * k)
    for l in range(b):
        for j in range(b):
            coord[a * a + l * b + j] = (2 * l + 1, 2 * j + 1)
    for i in range(a):
        for j in range(b):
            coord[N + i * b + j] = (2 * i, 2 * j + 1)
    return {"name": name, "rows": rows, "a": a, "b": b, "N": N, "M": M, "NS": NS, "hx": hx, "hz": hz,
            "cols": cols, "bx": bx, "bz": bz, "klog": klog, "ACT": ACT, "FOOT": FOOT, "WID": WID, "INV": INV,
            "psi": psi, "adjK": adjK, "coord": coord}


def frontier(order, adj):
    placed = front = width = 0
    for v in order:
        placed |= 1 << v
        front = (front | adj[v]) & ~placed
        width = max(width, front.bit_count())
    return width


def greedy_order(adj, NS):
    placed = front = 0
    order = []
    remaining = set(range(NS))
    while remaining:
        best = None
        for v in sorted(remaining):
            nf = ((front | adj[v]) & ~(placed | (1 << v))).bit_count()
            k = (nf, 0 if front >> v & 1 else 1, v)
            if best is None or k < best[0]:
                best = (k, v)
        v = best[1]
        order.append(v)
        remaining.discard(v)
        placed |= 1 << v
        front = (front | adj[v]) & ~placed
    return order


def phase_a_order(inst):
    adj, NS = inst["adjK"], inst["NS"]
    prod = sorted(range(NS), key=lambda v: (inst["coord"][v][0], inst["coord"][v][1], v))
    greedy = greedy_order(adj, NS)
    return prod if frontier(prod, adj) <= frontier(greedy, adj) else greedy


# ------------------------------------------------------------------ nuevo: logicos y caracteres
def nullspace(rows, n):
    R, pc = [], []
    for r in rows:
        for c, pr in zip(pc, R):
            if r >> c & 1:
                r ^= pr
        if r:
            c = r.bit_length() - 1
            for i in range(len(R)):
                if R[i] >> c & 1:
                    R[i] ^= r
            R.append(r)
            pc.append(c)
    pcs = set(pc)
    out = []
    for f in range(n):
        if f in pcs:
            continue
        v = 1 << f
        for c, pr in zip(pc, R):
            if pr >> f & 1:
                v |= 1 << c
        out.append(v)
    return out


def logical_reps(inst):
    N, hz, hx = inst["N"], inst["hz"], inst["hx"]
    ker = nullspace(hz, N)
    piv = dict(span(hx))
    reps = []
    for v in ker:
        if not inside(v, piv):
            reps.append(v)
            w = v
            while w:
                k = w.bit_length() - 1
                if k in piv:
                    w ^= piv[k]
                else:
                    piv[k] = w
                    break
    return reps


# ------------------------------------------------------------------ catalogo
def make_decoder(inst):
    INV, FOOT, WID, ACT = inst["INV"], inst["FOOT"], inst["WID"], inst["ACT"]
    RES = {}

    def decode(s0):
        CALLS[0] += 1
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
            if best < 0 or 2 * num < den:
                break
            cor ^= ACT[best]
            touched |= ACT[best]
            s ^= FOOT[best]
        got = (cor, touched)
        RES[s0] = got
        return got

    return decode, RES


def catalogue(inst, reps):
    NS, N, psi, adj, cols, bz = inst["NS"], inst["N"], inst["psi"], inst["adjK"], inst["cols"], inst["bz"]
    DATA = (1 << N) - 1
    klog = len(reps)
    group = []
    for m in range(1, 1 << klog):
        u = 0
        for i in range(klog):
            if m >> i & 1:
                u ^= reps[i]
        group.append(u)
    decode, RES = make_decoder(inst)
    POLY = {}          # Y -> {(na, nb): [n0, nT, chi_u1, ..., chi_u(2^k-1)]}
    naive_sil = {}     # (na, nb) -> historias silenciosas (todas, conexas o no)
    stats = {}
    i4_bad = 0
    for f in range(1, KMAX + 1):
        t = time.perf_counter()
        st = {"historias": 0, "decoraciones_validas": 0, "silenciosas_conexas_por_na_nb": {},
              "silenciosas_todas_por_na_nb": {}}
        for S in combinations(range(NS), f):
            sm = s = 0
            for i in S:
                sm |= 1 << i
                s ^= psi[i]
            st["historias"] += 1
            cor, touched = decode(s)
            E = sm & DATA
            na = E.bit_count()
            key = (na, f - na)
            r = E ^ cor
            passes = lin(cols, r) == 0
            silent = False
            chis = None
            if passes:
                silent = not inside(r, bz)
                chis = [-1 if (u & r).bit_count() & 1 else 1 for u in group]
                if silent != any(c < 0 for c in chis):
                    i4_bad += 1
                if silent:
                    ks = "%d,%d" % key
                    st["silenciosas_todas_por_na_nb"][ks] = st["silenciosas_todas_por_na_nb"].get(ks, 0) + 1
                    naive_sil[key] = naive_sil.get(key, 0) + 1
            Y = sm | touched
            seen = 0
            todo = Y & -Y
            while todo:
                bb = todo & -todo
                todo ^= bb
                seen |= bb
                todo |= adj[bb.bit_length() - 1] & Y & ~seen
            if seen != Y:
                continue
            st["decoraciones_validas"] += 1
            if silent:
                ks = "%d,%d" % key
                st["silenciosas_conexas_por_na_nb"][ks] = st["silenciosas_conexas_por_na_nb"].get(ks, 0) + 1
            slot = POLY.get(Y)
            if slot is None:
                slot = {}
                POLY[Y] = slot
            ent = slot.get(key)
            if ent is None:
                ent = [0] * (1 + (1 << klog))
                slot[key] = ent
            ent[0] += 1
            if passes:
                ent[1] += 1
                for i, c in enumerate(chis):
                    ent[2 + i] += c
            if CALLS[0] > LLAMADAS_MAX or time.perf_counter() - T0 > SEG_MAX:
                raise Stop("techo de llamadas o de tiempo alcanzado en el catalogo de " + inst["name"])
        st["segundos"] = round(time.perf_counter() - t, 3)
        stats[str(f)] = st
    RES.clear()
    return POLY, naive_sil, stats, i4_bad


# ------------------------------------------------------------------ gas truncado (estructura)
def gas_structure(inst, POLY, order, k):
    NS, N, adj = inst["NS"], inst["N"], inst["adjK"]
    pos = [0] * NS
    for i, v in enumerate(order):
        pos[v] = i

    def relabel(m):
        r = 0
        for x in bits(m):
            r |= 1 << pos[x]
        return r

    DATA = (1 << N) - 1
    DATAr = relabel(DATA)
    polys = []
    BYMIN = [[] for _ in range(NS)]
    for Y, slot in POLY.items():
        ents = {key: v for key, v in slot.items() if key[0] + key[1] <= k}
        if not ents:
            continue
        nb = Y
        for x in bits(Y):
            nb |= adj[x]
        Yr, nbr = relabel(Y), relabel(nb)
        pi = len(polys)
        polys.append(((Y & DATA).bit_count(), (Y & ~DATA).bit_count(), ents))
        BYMIN[(Yr & -Yr).bit_length() - 1].append((Yr, nbr, pi))
    ALL = (1 << NS) - 1
    states = {0}
    stack = [ALL & ~((1 << i) - 1) for i in range(NS)]
    while stack:
        J = stack.pop()
        if J in states:
            continue
        states.add(J)
        x = J & -J
        rest = J ^ x
        if rest not in states:
            stack.append(rest)
        for Yr, nbr, pi in BYMIN[x.bit_length() - 1]:
            if Yr & ~J == 0:
                Jn = J & ~nbr
                if Jn not in states:
                    stack.append(Jn)
    sts = sorted(states, key=lambda m: (m.bit_count(), m))
    index = {J: i for i, J in enumerate(sts)}
    rest_i = array("i", [0] * len(sts))
    xdata = array("b", [0] * len(sts))
    off = array("i", [0] * (len(sts) + 1))
    t_p, t_t, t_d, t_c = array("i"), array("i"), array("i"), array("i")
    for i, J in enumerate(sts):
        off[i] = len(t_p)
        if J == 0:
            continue
        x = J & -J
        rest_i[i] = index[J ^ x]
        xdata[i] = 1 if DATAr & x else 0
        for Yr, nbr, pi in BYMIN[x.bit_length() - 1]:
            if Yr & ~J == 0:
                extra = nbr & J & ~Yr
                t_p.append(pi)
                t_t.append(index[J & ~nbr])
                t_d.append((extra & DATAr).bit_count())
                t_c.append((extra & ~DATAr).bit_count())
    off[len(sts)] = len(t_p)
    return {"polys": polys, "n_states": len(sts), "rest": rest_i, "xdata": xdata, "off": off,
            "tp": t_p, "tt": t_t, "td": t_d, "tc": t_c, "top": index[ALL], "terms": len(t_p)}


def gas_eval(g, inst, p, q, NV):
    """Devuelve (Xi~_0, Xi~_T, [Xi~_T,u]) en el estado total y T = prod(D+1)."""
    N, M = inst["N"], inst["M"]
    Dp = F(1) / p - 1
    assert Dp.denominator == 1
    Dp = int(Dp)
    if q == 0:
        Dq, Tq = 1, 1
    else:
        Dq = F(1) / q - 1
        assert Dq.denominator == 1
        Dq = int(Dq)
        Tq = Dq + 1
    pp = [Dp ** i for i in range(N + 1)]
    pq = [Dq ** j for j in range(M + 1)]
    W0, WT = [], []
    for nd, nc, ents in g["polys"]:
        w = [0] * NV
        for (na, nb), v in ents.items():
            if q == 0 and nb > 0:
                continue
            mult = pp[nd - na] * pq[nc - nb]
            for i in range(NV):
                if v[i]:
                    w[i] += v[i] * mult
        W0.append(w[0])
        WT.append(tuple(w[1:]) if any(w[1:]) else None)
    NT = NV - 1
    V0 = [0] * g["n_states"]
    VT = [None] * g["n_states"]
    V0[0] = 1
    VT[0] = [1] * NT
    rest, xdata, off = g["rest"], g["xdata"], g["off"]
    tp, tt, td, tc = g["tp"], g["tt"], g["td"], g["tc"]
    pw = {}
    for i in range(1, g["n_states"]):
        Dx = Dp if xdata[i] else Dq
        ri = rest[i]
        a0 = Dx * V0[ri]
        aT = [Dx * v for v in VT[ri]]
        for j in range(off[i], off[i + 1]):
            pi, ti = tp[j], tt[j]
            kk = (td[j], tc[j])
            fct = pw.get(kk)
            if fct is None:
                fct = pp[kk[0]] * pq[kk[1]]
                pw[kk] = fct
            a0 += W0[pi] * V0[ti] * fct
            wt = WT[pi]
            if wt is not None:
                vt = VT[ti]
                for r in range(NT):
                    if wt[r]:
                        aT[r] += wt[r] * vt[r] * fct
        V0[i] = a0
        VT[i] = aT
    top = g["top"]
    T = (Dp + 1) ** N * Tq ** M
    return V0[top], VT[top], T


# ------------------------------------------------------------------ truncacion ingenua
def naive(inst, p, q, naive_sil, k):
    n, m = inst["N"], inst["M"]
    L = F(0)
    for (na, nb), c in naive_sil.items():
        if na + nb > k or (q == 0 and nb > 0):
            continue
        L += c * p ** na * (1 - p) ** (n - na) * q ** nb * (1 - q) ** (m - nb)
    le = F(0)
    for a in range(k + 1):
        for b in range(k + 1 - a):
            if q == 0 and b > 0:
                continue
            le += math.comb(n, a) * p ** a * (1 - p) ** (n - a) * math.comb(m, b) * q ** b * (1 - q) ** (m - b)
    return L, 1 - le


def err(L, eps):
    return eps / (2 * L + eps) if L > 0 else None


# ------------------------------------------------------------------ instancia completa
def run_instance(name, inst, order, cells, out_inst, expected_b=None):
    t = time.perf_counter()
    reps = logical_reps(inst)
    check("%s logicos: %d representantes = dimension logica" % (name, len(reps)), len(reps) == inst["klog"])
    check("%s logicos en ker H_Z" % name,
          all((u & g).bit_count() % 2 == 0 for u in reps for g in inst["hz"]))
    klog = len(reps)
    NV = 1 + (1 << klog)
    c0 = CALLS[0]
    POLY, naive_sil, stats, i4_bad = catalogue(inst, reps)
    TIMES[name + "_catalogo"] = round(time.perf_counter() - t, 3)
    check("I4 %s: silenciosa <=> algun caracter logico vale -1 (discrepancias)" % name, i4_bad == 0, i4_bad)
    out_inst["catalogo"] = {"por_f": stats, "polimeros_distintos": len(POLY), "llamadas": CALLS[0] - c0,
                            "silenciosas_todas": {"%d,%d" % kk: v for kk, v in sorted(naive_sil.items())},
                            "qubits_logicos": klog}
    if expected_b is not None:
        pf = expected_b["B1"]["por_f"]
        for f in ("1", "2", "3"):
            check("%s reproduce FASE_B f=%s: historias, decoraciones y silenciosas conexas" % (name, f),
                  stats[f]["historias"] == pf[f]["historias"]
                  and stats[f]["decoraciones_validas"] == pf[f]["decoraciones_validas"]
                  and stats[f]["silenciosas_conexas_por_na_nb"] == pf[f]["silenciosas_por_na_nb"],
                  {"obtenido": [stats[f]["historias"], stats[f]["decoraciones_validas"],
                                stats[f]["silenciosas_conexas_por_na_nb"]]})
        check("%s reproduce FASE_B: polimeros distintos" % name,
              len(POLY) == expected_b["B1"]["polimeros_distintos"], len(POLY))
    table = {}
    for k in range(1, KMAX + 1):
        t = time.perf_counter()
        g = gas_structure(inst, POLY, order, k)
        TIMES["%s_estructura_k%d" % (name, k)] = round(time.perf_counter() - t, 3)
        row = {"estados": g["n_states"], "terminos": g["terms"], "polimeros": len(g["polys"]), "celdas": {}}
        t = time.perf_counter()
        for (p, q) in cells:
            X0, XT, T = gas_eval(g, inst, p, q, NV)
            eps_c = 1 - F(X0, T)
            sumu = sum(XT[1:])
            Lc = F(XT[0] * ((1 << klog) - 1) - sumu, (1 << klog) * T)
            P0 = F(XT[0], T)
            Li, eps_i = naive(inst, p, q, naive_sil, k)
            check("I2 %s k=%d p=%s q=%s: 0 <= eps_comp <= eps_ing" % (name, k, p, q), 0 <= eps_c <= eps_i)
            check("I5 %s k=%d p=%s q=%s: L_ing <= L_comp" % (name, k, p, q), Li <= Lc)
            row["celdas"]["%s|%s" % (p, q)] = {"eps_comp": eps_c, "L_comp": Lc, "P_0_inf": P0,
                                              "eps_ing": eps_i, "L_ing": Li}
        TIMES["%s_gas_k%d" % (name, k)] = round(time.perf_counter() - t, 3)
        table[k] = row
        if time.perf_counter() - T0 > SEG_MAX:
            raise Stop("techo de tiempo alcanzado en el gas de " + name)
    for (p, q) in cells:
        kk = "%s|%s" % (p, q)
        Ls = [table[k]["celdas"][kk]["L_comp"] for k in (1, 2, 3)]
        Es = [table[k]["celdas"][kk]["eps_comp"] for k in (1, 2, 3)]
        Li = [table[k]["celdas"][kk]["L_ing"] for k in (1, 2, 3)]
        check("I3 %s %s: L_comp y L_ing no decrecen; eps_comp no crece" % (name, kk),
              Ls[0] <= Ls[1] <= Ls[2] and Li[0] <= Li[1] <= Li[2] and Es[0] >= Es[1] >= Es[2])
    if expected_b is not None and "B3" in expected_b and "k_le_3" in expected_b["B3"]:
        b3 = expected_b["B3"]["k_le_3"]
        check("%s reproduce FASE_B: estados y terminos de la DP con k<=3" % name,
              table[3]["estados"] == b3["estados"] and table[3]["terminos"] == b3["terminos"],
              [table[3]["estados"], table[3]["terminos"]])
    return table


def summarize(name, inst, table, cells):
    NS = inst["NS"]
    calls_k = {k: sum(math.comb(NS, f) for f in range(1, k + 1)) for k in (1, 2, 3, 4)}
    res = {}
    for (p, q) in cells:
        kk = "%s|%s" % (p, q)
        rows = {}
        for k in (1, 2, 3):
            c = table[k]["celdas"][kk]
            ec, ei = err(c["L_comp"], c["eps_comp"]), err(c["L_ing"], c["eps_ing"])
            R = c["eps_ing"] / c["eps_comp"] if c["eps_comp"] > 0 else None
            rows[str(k)] = {"eps_comp": show(c["eps_comp"]), "eps_ing": show(c["eps_ing"]),
                            "L_comp": show(c["L_comp"]), "L_ing": show(c["L_ing"]),
                            "P_0_cota_inferior": show(c["P_0_inf"]),
                            "e_comp": show(ec) if ec is not None else None,
                            "e_ing": show(ei) if ei is not None else None,
                            "R": show(R) if R is not None else None,
                            "llamadas": calls_k[k],
                            "_e_comp": ec, "_e_ing": ei, "_R": R}
        _, eps4 = naive(inst, p, q, {}, 4)
        L3i = table[3]["celdas"][kk]["L_ing"]
        bound4 = eps4 / (2 * L3i + eps4) if L3i > 0 else None
        kst = {}
        for tg in TARGETS:
            kc = next((k for k in (1, 2, 3) if rows[str(k)]["_e_comp"] is not None
                       and rows[str(k)]["_e_comp"] <= tg), None)
            ki = next((k for k in (1, 2, 3) if rows[str(k)]["_e_ing"] is not None
                       and rows[str(k)]["_e_ing"] <= tg), None)
            if ki is None:
                ki = 4 if (bound4 is not None and bound4 <= tg) else ">4"
            kst[str(tg)] = {"componentes": kc if kc is not None else ">3", "ingenua": ki,
                            "llamadas_componentes": calls_k[kc] if kc is not None else None,
                            "llamadas_ingenua": calls_k[ki] if isinstance(ki, int) else None}
        res[kk] = {"p": str(p), "q": str(q), "por_k": rows, "eps_ing_4": show(eps4),
                   "cota_e_ing_4": show(bound4) if bound4 is not None else None, "k_estrella": kst}
    return res


def strip_private(res):
    for cell in res.values():
        for row in cell["por_k"].values():
            for kk in [x for x in row if x.startswith("_")]:
                del row[kk]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--recursos", required=True)
    ap.add_argument("--solo", choices=["PATH25"], default=None)
    args = ap.parse_args()
    for pth in (args.output, args.recursos):
        if os.path.exists(pth):
            raise SystemExit("La salida ya existe; use un nombre nuevo: " + pth)
    me = sha_file(os.path.abspath(__file__))
    out = {"unit": "D1-CE", "plan": "PLAN_CONTROL.txt", "programa_sha256": me, "solo": args.solo,
           "instancias": {}}
    status = "PASS"
    try:
        run(args, out)
    except Stop as ex:
        status = "FAIL"
        out["detenido_en"] = str(ex)
    out["status"] = status
    out["checks"] = CHECKS
    out["llamadas_totales"] = CALLS[0]
    with open(args.output, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(out, ensure_ascii=False, indent=1, sort_keys=True, default=str) + "\n")
    TIMES["total"] = round(time.perf_counter() - T0, 3)
    with open(args.recursos, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps({"nota": "tiempos de reloj en segundos y memoria; no deterministas; no son resultado",
                            "tiempos": TIMES, "pico_memoria_bytes": peak_memory_bytes(),
                            "llamadas_totales": CALLS[0]}, ensure_ascii=False, indent=1) + "\n")
    print("%s controles=%d llamadas=%d t=%.1fs" % (status, len(CHECKS), CALLS[0], TIMES["total"]))
    if status != "PASS":
        print("DETENIDO EN:", out.get("detenido_en"))
        sys.exit(1)


def run(args, out):
    check("FASE_A.json de D1 con su SHA-256", sha_file(D1 + "FASE_A.json") == SHA_FASE_A)
    check("FASE_B.json de D1 con su SHA-256", sha_file(D1 + "FASE_B.json") == SHA_FASE_B)
    check("ESTIMACION_H6_U4.json con su SHA-256", sha_file(U4 + "ESTIMACION_H6_U4.json") == SHA_U4)
    A = json.load(open(D1 + "FASE_A.json", encoding="utf-8"))["instancias"]
    B = json.load(open(D1 + "FASE_B.json", encoding="utf-8"))["instancias"]
    U = json.load(open(U4 + "ESTIMACION_H6_U4.json", encoding="utf-8"))

    insts = [("PATH25",) + path_seed(4)]
    if args.solo is None:
        insts += [("S1_L6",) + path_seed(6), ("S1_L7",) + path_seed(7), ("S2_34_6x8",) + regular_seed_34()]
    cells_grid = [(p, q) for p in GRID_P for q in (p, p / 2, F(0))]
    decisions_input = {}
    for name, rows, a in insts:
        t = time.perf_counter()
        inst = build(name, rows, a)
        mats = json.dumps({"hx": inst["hx"], "hz": inst["hz"]}, sort_keys=True).encode()
        check("%s matrices iguales a FASE_A de D1" % name,
              hashlib.sha256(mats).hexdigest() == A[name]["matrices_sha256"])
        order = list(range(inst["NS"])) if name == "PATH25" else phase_a_order(inst)
        cells = U4_PAIRS if name == "PATH25" else cells_grid
        oi = {"N": inst["N"], "M": inst["M"], "sitios": inst["NS"], "qubits_logicos": inst["klog"],
              "orden": "original" if name == "PATH25" else "fase A de D1"}
        exp_b = B[name] if name != "PATH25" else None
        table = run_instance(name, inst, order, cells, oi, exp_b)
        oi["gas_por_k"] = {str(k): {"estados": table[k]["estados"], "terminos": table[k]["terminos"],
                                    "polimeros": table[k]["polimeros"]} for k in (1, 2, 3)}
        if name == "PATH25":
            pf = B["PATH25"]["B1"]["por_f"]
            st = oi["catalogo"]["por_f"]
            check("PATH25 reproduce FASE_B: historias, decoraciones y silenciosas conexas por f",
                  all(st[f]["historias"] == pf[f]["historias"]
                      and st[f]["decoraciones_validas"] == pf[f]["decoraciones_validas"]
                      and st[f]["silenciosas_conexas_por_na_nb"] == pf[f]["silenciosas_por_na_nb"]
                      for f in ("1", "2", "3")))
            check("PATH25 polimeros, estados y terminos del gas por k = 37/387/3918, 234/238/240 y 207/1168/7179",
                  [table[k]["polimeros"] for k in (1, 2, 3)] == [37, 387, 3918]
                  and [table[k]["estados"] for k in (1, 2, 3)] == [234, 238, 240]
                  and [table[k]["terminos"] for k in (1, 2, 3)] == [207, 1168, 7179])
            for k in (1, 2, 3):
                urow = U["tabla"][k - 1]
                assert urow["k"] == k
                for (p, q) in U4_PAIRS:
                    kk = "%s|%s" % (p, q)
                    uc = urow["tasas"][kk]
                    a1, b1 = uc["eps_k"]["exacto"].split("/")
                    a2, b2 = uc["P_sil"]["L"]["exacto"].split("/")
                    mine = table[k]["celdas"][kk]
                    check("PATH25 k=%d %s: eps_k y L de P_sil identicos a U4" % (k, kk),
                          mine["eps_comp"] == F(int(a1), int(b1)) and mine["L_comp"] == F(int(a2), int(b2)))
        res = summarize(name, inst, table, cells)
        decisions_input[name] = res
        oi["celdas"] = res
        oi["segundos"] = round(time.perf_counter() - t, 3)
        out["instancias"][name] = oi
        print("%-10s ok %.1fs llamadas=%d" % (name, oi["segundos"], CALLS[0]), flush=True)
    if args.solo is None:
        out["decisiones"] = decide(decisions_input)
        out["sesgo_prediccion_D1"] = bias(decisions_input)
    for res in decisions_input.values():
        strip_private(res)


def decide(R):
    def r3(name, kk):
        return R[name][kk]["por_k"]["3"]["_R"]

    keys = list(R["S1_L6"].keys())
    out = {}
    name, p, q = PROPOSED
    kk = "%s|%s" % (p, q)
    cell = R[name][kk]
    e3 = cell["por_k"]["3"]["_e_comp"]
    ok = e3 is not None and e3 <= F(1, 10)
    ki = cell["k_estrella"][str(F(1, 10))]["ingenua"]
    out["D_a"] = {"celda": "%s %s" % (name, kk), "e_comp_3": show(e3) if e3 is not None else None,
                  "e_ing_3": show(cell["por_k"]["3"]["_e_ing"]),
                  "componentes_alcanzan_el_10_por_ciento_con_k3": ok, "k_estrella_ingenua_10_por_ciento": ki,
                  "ahorro_llamadas": ("125671 frente a 2798341" if ok and ki == 4 else None)}
    out["D_b"] = {}
    for name in ("S1_L6", "S1_L7"):
        vals = [r3(name, kk) for kk in keys]
        n_rel = sum(1 for v in vals if v >= 10)
        n_nul = sum(1 for v in vals if v < 2)
        cat = "relevante" if n_rel == len(vals) else ("nula" if n_nul > 0 else "zona gris")
        out["D_b"][name] = {"categoria": cat, "celdas_R3_mayor_o_igual_10": n_rel,
                            "celdas_R3_menor_que_2": n_nul, "celdas_R3_entre_2_y_10": len(vals) - n_rel - n_nul,
                            "R3_min": dec(min(vals)), "R3_max": dec(max(vals))}
    out["D_b"]["S1_conjunto"] = ("relevante" if all(out["D_b"][n]["categoria"] == "relevante"
                                                   for n in ("S1_L6", "S1_L7")) else "no relevante en conjunto")
    grow = [r3("S1_L7", kk) > r3("S1_L6", kk) for kk in keys]
    out["D_c"] = {"crece_con_N": all(grow), "celdas_que_crecen": sum(grow), "celdas": len(grow)}
    vals = [r3("S2_34_6x8", kk) for kk in R["S2_34_6x8"]]
    neg_ok = all(v < 2 for v in vals)
    out["control_negativo_S2"] = {"R3_menor_que_2_en_todas": neg_ok, "R3_min": dec(min(vals)),
                                  "R3_max": dec(max(vals))}
    if not neg_ok:
        out["control_negativo_S2"]["accion"] = "detener la lectura y revisar antes de informar (plan)"
    return out


def bias(R):
    """Sesgo de la prediccion de D1: e_comp(3) predicho (FASE_C.json) / exacto."""
    path = D1 + "FASE_C.json"
    C = json.load(open(path, encoding="utf-8"))["instancias"]
    out = {"FASE_C_sha256": sha_file(path), "celdas": {}}
    for name in ("S1_L6", "S1_L7", "S2_34_6x8"):
        for cell in C[name]["tasas"]:
            kk = "%s|%s" % (cell["p"], cell["q"])
            ex = R[name].get(kk)
            if ex is None:
                continue
            pred = cell["error_predicho_por_k"].get("3")
            exact = ex["por_k"]["3"]["_e_comp"]
            out["celdas"]["%s %s" % (name, kk)] = {
                "e_comp_3_predicho": pred, "e_comp_3_exacto": dec(exact) if exact is not None else None,
                "cociente_predicho_sobre_exacto": (pred / dec(exact)) if (pred and exact) else None}
    return out


if __name__ == "__main__":
    main()
