#!/usr/bin/env python3
"""B1 (PLAN_B1.txt, sellado en DEC-040): certificado con k=4, exacto, en S1 con
L=6, 7 y 8, con PATH25 como control positivo.

Importa control_exacto.py de D1-CE tras comprobar su SHA-256 y usa sin cambios
su constructor, decoder, logicos, truncacion ingenua y DP de referencia. Piezas
nuevas: (1) catalogo en paralelo, con las mismas historias que el de D1-CE,
repartidas por primer sitio; (2) DP del gas sin guardar terminos, con mascaras
de exclusion por sitio. Ambas se contrastan con las de D1-CE.

Uso (desde D:\\paper ia cuantica\\qec_revision_matematica\\05_COMPROBACIONES):
  python -I -B <dir>/b1_exacto_k4.py --output <nuevo.json> --recursos <nuevo.json> --tmp <carpeta temporal>
"""
import argparse
import hashlib
import importlib.util
import json
import math
import multiprocessing as mp
import os
import pickle
import sys
import time
from fractions import Fraction as F
from itertools import combinations

T0 = time.perf_counter()
CE_DIR = "D1_CONTROL_EXACTO_Y_METODOLOGIA_20261113/"
B0_DIR = "B0_ETAPA0_MAPA_REGULAR_Y_S1_L8_20261113/"
SHA_CE = "c003050183a14364750a33784c7c815c0c07b5814f552a8b52a14d142db3306b"
SHA_RES_CE = "334d6485968c914c299bbc7e5362cf49f4b975b740d47d78bccfaec47255024d"
SHA_E1 = "1d444fead5056eccdeb036d1b6a8bd8403753a90a051af224296fa714f31bb59"
SHA_PLAN = "69daf6a5061ab430b437f6d2febea343256bb8fc81fae598cf74b36e2f2c5be8"
GIB = 1 << 30
LLAMADAS_MAX = 55_000_000
SEG_MAX = 6 * 3600.0
MEM_TOTAL = 7.5 * GIB
MEM_MAIN = 2 * GIB
MEM_CAT = 0.75 * GIB
W_CAT = 6
W_GAS_MAX = 3
GATE_STATES = 5_000_000
GATE_TERMS = 2_500_000_000
CAL_TERMS = 20_000_000
KMAX = 4
TARGETS = (F(1, 20), F(1, 10), F(1, 4))
TIMES = {}
WORKER_PEAKS = {}
CALLS_W = [0]
CE = [None]


def sha_file(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def load_ce():
    path = CE_DIR + "control_exacto.py"
    if sha_file(path) != SHA_CE:
        raise SystemExit("control_exacto.py no es el sellado en D1-CE")
    spec = importlib.util.spec_from_file_location("control_exacto", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def mem_now(peak=False):
    """Working set actual (o maximo) del proceso en bytes (Windows); None si no se puede medir."""
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
            return int(pmc.PeakWorkingSetSize if peak else pmc.WorkingSetSize)
    except Exception:
        pass
    return None


def mem_peak():
    return mem_now(peak=True)


def guard(where):
    ce = CE[0]
    if time.perf_counter() - T0 > SEG_MAX:
        raise ce.Stop("techo de tiempo alcanzado en " + where)
    m = mem_now()
    if m is not None and m > MEM_MAIN:
        raise ce.Stop("techo de memoria del proceso principal alcanzado en " + where)


def frac(s):
    a, b = s.split("/")
    return F(int(a), int(b))


# ------------------------------------------------------------------ catalogo en paralelo
_W = {}


def cat_init(name, rows, a, reps):
    ce = load_ce()
    _W["ce"] = ce
    _W["inst"] = ce.build(name, rows, a)
    _W["reps"] = reps


def cat_task(task):
    f, firsts = task
    ce, inst, reps = _W["ce"], _W["inst"], _W["reps"]
    NS, N, psi, adj, cols, bz = inst["NS"], inst["N"], inst["psi"], inst["adjK"], inst["cols"], inst["bz"]
    lin, inside = ce.lin, ce.inside
    DATA = (1 << N) - 1
    klog = len(reps)
    group = []
    for m in range(1, 1 << klog):
        u = 0
        for i in range(klog):
            if m >> i & 1:
                u ^= reps[i]
        group.append(u)
    decode, RES = ce.make_decoder(inst)
    c0 = ce.CALLS[0]
    POLY, naive_sil = {}, {}
    st = {"historias": 0, "decoraciones_validas": 0, "silenciosas_conexas_por_na_nb": {},
          "silenciosas_todas_por_na_nb": {}}
    i4_bad = 0
    for first in firsts:
        for rest in combinations(range(first + 1, NS), f - 1):
            sm = 1 << first
            s = psi[first]
            for i in rest:
                sm |= 1 << i
                s ^= psi[i]
            st["historias"] += 1
            if len(RES) > 300000:
                RES.clear()
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
        m = mem_now()
        if m is not None and m > MEM_CAT:
            raise RuntimeError("techo de memoria de un proceso del catalogo (%d bytes)" % m)
    RES.clear()
    return {"POLY": POLY, "naive_sil": naive_sil, "st": st, "i4_bad": i4_bad, "calls": ce.CALLS[0] - c0,
            "peak": ce.peak_memory_bytes(), "pid": os.getpid()}


def split_firsts(NS, f, nbins):
    items = sorted(((math.comb(NS - 1 - i, f - 1), i) for i in range(NS - f + 1)), reverse=True)
    bins = [[0, []] for _ in range(nbins)]
    for w, i in items:
        b = min(bins, key=lambda x: x[0])
        b[0] += w
        b[1].append(i)
    return [(f, sorted(b[1])) for b in bins if b[1]]


def merge_poly(POLY, part):
    for Y, slot in part.items():
        g = POLY.get(Y)
        if g is None:
            POLY[Y] = slot
            continue
        for key, ent in slot.items():
            e = g.get(key)
            if e is None:
                g[key] = ent
            else:
                for j in range(len(e)):
                    e[j] += ent[j]


def catalogue_parallel(name, rows, a, inst, reps, kmax, ctx):
    ce = CE[0]
    NS = inst["NS"]
    POLY, naive_sil, stats = {}, {}, {}
    i4_bad = 0
    calls = 0
    with ctx.Pool(W_CAT, initializer=cat_init, initargs=(name, rows, a, reps)) as pool:
        for f in range(1, kmax + 1):
            t = time.perf_counter()
            if ce.CALLS[0] + CALLS_W[0] + math.comb(NS, f) > LLAMADAS_MAX:
                raise ce.Stop("techo de llamadas: el catalogo de %s con f=%d no cabe" % (name, f))
            st = {"historias": 0, "decoraciones_validas": 0, "silenciosas_conexas_por_na_nb": {},
                  "silenciosas_todas_por_na_nb": {}}
            try:
                for res in pool.imap_unordered(cat_task, split_firsts(NS, f, 4 * W_CAT)):
                    merge_poly(POLY, res["POLY"])
                    for kk, v in res["naive_sil"].items():
                        naive_sil[kk] = naive_sil.get(kk, 0) + v
                    st["historias"] += res["st"]["historias"]
                    st["decoraciones_validas"] += res["st"]["decoraciones_validas"]
                    for fld in ("silenciosas_conexas_por_na_nb", "silenciosas_todas_por_na_nb"):
                        for kk, v in res["st"][fld].items():
                            st[fld][kk] = st[fld].get(kk, 0) + v
                    i4_bad += res["i4_bad"]
                    calls += res["calls"]
                    CALLS_W[0] += res["calls"]
                    pk = WORKER_PEAKS.setdefault("catalogo", {})
                    pk[res["pid"]] = max(pk.get(res["pid"], 0), res["peak"] or 0)
                    guard("catalogo de " + name)
            except Exception as ex:
                raise ce.Stop("catalogo de %s: %s" % (name, ex))
            st["silenciosas_conexas_por_na_nb"] = dict(sorted(st["silenciosas_conexas_por_na_nb"].items()))
            st["silenciosas_todas_por_na_nb"] = dict(sorted(st["silenciosas_todas_por_na_nb"].items()))
            stats[str(f)] = st
            TIMES["%s_catalogo_f%d" % (name, f)] = round(time.perf_counter() - t, 3)
            print("  %s catalogo f=%d: %d historias, %d decoraciones validas (%.1fs)"
                  % (name, f, st["historias"], st["decoraciones_validas"], time.perf_counter() - t), flush=True)
    return POLY, naive_sil, stats, i4_bad, calls


# ------------------------------------------------------------------ DP del gas sin guardar terminos
def prep_gas(inst, POLY, order, k):
    """Misma seleccion de polimeros y mismo reetiquetado que gas_structure de D1-CE."""
    ce = CE[0]
    NS, N, adj = inst["NS"], inst["N"], inst["adjK"]
    pos = [0] * NS
    for i, v in enumerate(order):
        pos[v] = i

    def relabel(m):
        r = 0
        for x in ce.bits(m):
            r |= 1 << pos[x]
        return r

    DATA = (1 << N) - 1
    polys, YR, NB = [], [], []
    BY = [[] for _ in range(NS)]
    for Y, slot in POLY.items():
        ents = {key: v for key, v in slot.items() if key[0] + key[1] <= k}
        if not ents:
            continue
        nb = Y
        for x in ce.bits(Y):
            nb |= adj[x]
        Yr, nbr = relabel(Y), relabel(nb)
        pi = len(polys)
        polys.append(((Y & DATA).bit_count(), (Y & ~DATA).bit_count(), ents))
        YR.append(Yr)
        NB.append(nbr)
        BY[(Yr & -Yr).bit_length() - 1].append(pi)
    return {"polys": polys, "YR": YR, "NB": NB, "BY": BY, "DATAr": relabel(DATA), "NS": NS, "N": N,
            "M": inst["M"]}


def anchor_masks(G):
    """Para cada sitio ancla x: ventana, mascara total y, por sitio y de la ventana, la mascara de los
    polimeros anclados en x que contienen y. Un polimero anclado en el minimo x de J esta contenido en J
    si y solo si ninguno de sus sitios es un hueco de J en la ventana."""
    NS, YR, BY = G["NS"], G["YR"], G["BY"]
    WIN, FULL, CM = [0] * NS, [0] * NS, [None] * NS
    for x in range(NS):
        lst = BY[x]
        cm = {}
        win = 0
        for j, pi in enumerate(lst):
            m = YR[pi] & ~(1 << x)
            win |= m
            bj = 1 << j
            while m:
                b = m & -m
                m ^= b
                cm[b] = cm.get(b, 0) | bj
        WIN[x], FULL[x], CM[x] = win, (1 << len(lst)) - 1, cm
    return WIN, FULL, CM


def contained(J, x, WIN, FULL, CM):
    holes = WIN[x] & ~J
    excl = 0
    cm = CM[x]
    while holes:
        b = holes & -holes
        holes ^= b
        excl |= cm[b]
    return FULL[x] & ~excl


def explore(G, A, limits=None):
    """Mismo conjunto de estados que gas_structure de D1-CE. Devuelve (estados ordenados, terminos) o
    (None, info) si se supera la compuerta."""
    NS, NB, BY = G["NS"], G["NB"], G["BY"]
    WIN, FULL, CM = A
    ALL = (1 << NS) - 1
    states = {0}
    stack = [ALL & ~((1 << i) - 1) for i in range(NS)]
    terms = 0
    n = 0
    while stack:
        J = stack.pop()
        if J in states:
            continue
        states.add(J)
        x = J & -J
        xi = x.bit_length() - 1
        rest = J ^ x
        if rest not in states:
            stack.append(rest)
        c = contained(J, xi, WIN, FULL, CM)
        lst = BY[xi]
        while c:
            b = c & -c
            c ^= b
            terms += 1
            Jn = J & ~NB[lst[b.bit_length() - 1]]
            if Jn not in states:
                stack.append(Jn)
        n += 1
        if limits is not None and n % 20000 == 0:
            if len(states) > limits["states"] or terms > limits["terms"]:
                return None, {"estados_al_menos": len(states), "terminos_al_menos": terms,
                              "motivo": "estados o terminos por encima de la compuerta"}
            guard("exploracion de estados")
    if limits is not None and (len(states) > limits["states"] or terms > limits["terms"]):
        return None, {"estados_al_menos": len(states), "terminos_al_menos": terms,
                      "motivo": "estados o terminos por encima de la compuerta"}
    return sorted(states, key=lambda m: (m.bit_count(), m)), terms


def poly_weights(G, p, q, NV, pp, pq):
    W0, WT = [], []
    for nd, nc, ents in G["polys"]:
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
    return W0, WT


def eval_fly(G, A, sts, index, p, q, NV, cal_terms=None, mem_limit=None):
    """Misma recurrencia y mismos enteros escalados que gas_eval de D1-CE, recalculando los terminos."""
    N, M, NB, YR, BY, DATAr = G["N"], G["M"], G["NB"], G["YR"], G["BY"], G["DATAr"]
    WIN, FULL, CM = A
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
    W0, WT = poly_weights(G, p, q, NV, pp, pq)
    NBO = [nb & ~y for nb, y in zip(NB, YR)]
    NT = NV - 1
    ns = len(sts)
    V0 = [0] * ns
    VT = [None] * ns
    V0[0] = 1
    VT[0] = [1] * NT
    pw = {}
    count = 0
    t = time.perf_counter()
    for i in range(1, ns):
        J = sts[i]
        x = J & -J
        xi = x.bit_length() - 1
        Dx = Dp if DATAr & x else Dq
        ri = index[J ^ x]
        a0 = Dx * V0[ri]
        aT = [Dx * v for v in VT[ri]]
        c = contained(J, xi, WIN, FULL, CM)
        lst = BY[xi]
        while c:
            b = c & -c
            c ^= b
            pi = lst[b.bit_length() - 1]
            ti = index[J & ~NB[pi]]
            extra = NBO[pi] & J
            td = (extra & DATAr).bit_count()
            kk = (td, extra.bit_count() - td)
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
            count += 1
        V0[i] = a0
        VT[i] = aT
        if cal_terms is not None and count >= cal_terms:
            return {"terminos": count, "estados": i, "segundos": time.perf_counter() - t}
        if mem_limit is not None and i % 20000 == 0:
            m = mem_now()
            if m is not None and m > mem_limit:
                raise RuntimeError("techo de memoria de un proceso de la DP (%d bytes)" % m)
    if cal_terms is not None:
        return {"terminos": count, "estados": ns - 1, "segundos": time.perf_counter() - t}
    top = ns - 1
    assert sts[top] == (1 << G["NS"]) - 1
    T = (Dp + 1) ** N * Tq ** M
    return V0[top], VT[top], T


def gas_init(path):
    with open(path, "rb") as fh:
        G, sts = pickle.load(fh)
    _W["G"] = G
    _W["sts"] = sts
    _W["index"] = {J: i for i, J in enumerate(sts)}
    _W["A"] = anchor_masks(G)


def gas_task(task):
    p, q, NV, mem_limit = task
    t = time.perf_counter()
    X0, XT, T = eval_fly(_W["G"], _W["A"], _W["sts"], _W["index"], p, q, NV, mem_limit=mem_limit)
    return {"p": p, "q": q, "X0": X0, "XT": XT, "T": T, "segundos": round(time.perf_counter() - t, 3),
            "peak": mem_peak(), "pid": os.getpid()}


def finish_cell(p, q, X0, XT, T, klog):
    eps_c = 1 - F(X0, T)
    sumu = sum(XT[1:])
    Lc = F(XT[0] * ((1 << klog) - 1) - sumu, (1 << klog) * T)
    P0 = F(XT[0], T)
    return eps_c, Lc, P0


def run_dp(name, k, G, sts, cells, NV, klog, nproc, ctx, tmpdir, mem_limit=None):
    """Evalua las celdas con la DP nueva; devuelve {celda: (eps_comp, L_comp, P0)} y segundos por celda."""
    out, secs = {}, {}
    if nproc <= 1:
        A = anchor_masks(G)
        index = {J: i for i, J in enumerate(sts)}
        for (p, q) in cells:
            t = time.perf_counter()
            X0, XT, T = eval_fly(G, A, sts, index, p, q, NV)
            out["%s|%s" % (p, q)] = finish_cell(p, q, X0, XT, T, klog)
            secs["%s|%s" % (p, q)] = round(time.perf_counter() - t, 3)
            guard("DP de %s k=%d" % (name, k))
        return out, secs
    path = os.path.join(tmpdir, "b1_gas_%s_k%d.pkl" % (name, k))
    with open(path, "wb") as fh:
        pickle.dump((G, sts), fh, protocol=pickle.HIGHEST_PROTOCOL)
    try:
        with ctx.Pool(nproc, initializer=gas_init, initargs=(path,)) as pool:
            try:
                for res in pool.imap_unordered(gas_task, [(p, q, NV, mem_limit) for (p, q) in cells]):
                    kk = "%s|%s" % (res["p"], res["q"])
                    out[kk] = finish_cell(res["p"], res["q"], res["X0"], res["XT"], res["T"], klog)
                    secs[kk] = res["segundos"]
                    pk = WORKER_PEAKS.setdefault("dp_%s_k%d" % (name, k), {})
                    pk[res["pid"]] = max(pk.get(res["pid"], 0), res["peak"] or 0)
                    print("    %s k=%d celda %s lista (%.0fs)" % (name, k, kk, res["segundos"]), flush=True)
                    guard("DP de %s k=%d" % (name, k))
            except Exception as ex:
                raise CE[0].Stop("DP de %s k=%d: %s" % (name, k, ex))
    finally:
        os.remove(path)
    return out, secs


# ------------------------------------------------------------------ resumen y decisiones
def summarize4(name, inst, table, cells, kmax, evaluated):
    ce = CE[0]
    NS = inst["NS"]
    calls_k = {k: sum(math.comb(NS, f) for f in range(1, k + 1)) for k in range(1, kmax + 2)}
    res = {}
    for (p, q) in cells:
        kk = "%s|%s" % (p, q)
        rows = {}
        for k in range(1, kmax + 1):
            c = table[k]["celdas"][kk]
            comp_ok = evaluated[k]
            ec = ce.err(c["L_comp"], c["eps_comp"]) if comp_ok else None
            ei = ce.err(c["L_ing"], c["eps_ing"])
            R = (c["eps_ing"] / c["eps_comp"] if comp_ok and c["eps_comp"] > 0 else None)
            rows[str(k)] = {"eps_comp": ce.show(c["eps_comp"]) if comp_ok else None,
                            "eps_ing": ce.show(c["eps_ing"]),
                            "L_comp": ce.show(c["L_comp"]) if comp_ok else None,
                            "L_ing": ce.show(c["L_ing"]),
                            "P_0_cota_inferior": ce.show(c["P_0_inf"]) if comp_ok else None,
                            "e_comp": ce.show(ec) if ec is not None else None,
                            "e_ing": ce.show(ei) if ei is not None else None,
                            "R": ce.show(R) if R is not None else None,
                            "componentes_evaluadas": comp_ok, "llamadas": calls_k[k],
                            "_e_comp": ec, "_e_ing": ei, "_R": R}
        _, epsn = ce.naive(inst, p, q, {}, kmax + 1)
        Lk = table[kmax]["celdas"][kk]["L_ing"]
        boundn = epsn / (2 * Lk + epsn) if Lk > 0 else None
        kst = {}
        for tg in TARGETS:
            kc = next((k for k in range(1, kmax + 1) if rows[str(k)]["_e_comp"] is not None
                       and rows[str(k)]["_e_comp"] <= tg), None)
            ki = next((k for k in range(1, kmax + 1) if rows[str(k)]["_e_ing"] is not None
                       and rows[str(k)]["_e_ing"] <= tg), None)
            if ki is None:
                ki = (kmax + 1) if (boundn is not None and boundn <= tg) else ">%d" % (kmax + 1)
            kst[str(tg)] = {"componentes": kc if kc is not None else ">%d" % kmax, "ingenua": ki,
                            "llamadas_componentes": calls_k[kc] if kc is not None else None,
                            "llamadas_ingenua": calls_k[ki] if isinstance(ki, int) else None}
        res[kk] = {"p": str(p), "q": str(q), "por_k": rows, "eps_ing_%d" % (kmax + 1): ce.show(epsn),
                   "cota_e_ing_%d_informativa" % (kmax + 1): ce.show(boundn) if boundn is not None else None,
                   "k_estrella": kst}
    return res


def strip_private(res):
    for cell in res.values():
        for row in cell["por_k"].values():
            for kk in [x for x in row if x.startswith("_")]:
                del row[kk]


def m8_category(vals):
    n_rel = sum(1 for v in vals if v >= 10)
    n_nul = sum(1 for v in vals if v < 2)
    cat = "relevante" if n_rel == len(vals) else ("nula" if n_nul > 0 else "zona gris")
    return {"categoria": cat, "celdas_mayor_o_igual_10": n_rel, "celdas_menor_que_2": n_nul,
            "celdas_entre_2_y_10": len(vals) - n_rel - n_nul}


def decide(R, keys, k4_done, nsil):
    ce = CE[0]
    S1 = ("S1_L6", "S1_L7", "S1_L8")
    out = {"D_f1": {}, "D_f2": {}}
    any_adv = False
    for name in S1:
        if not k4_done[name]:
            out["D_f1"][name] = {"evaluada": False}
            continue
        row = {"evaluada": True}
        for tg in TARGETS:
            adv = [kk for kk in keys if R[name][kk]["por_k"]["4"]["_e_comp"] is not None
                   and R[name][kk]["por_k"]["4"]["_e_comp"] <= tg
                   and (R[name][kk]["por_k"]["4"]["_e_ing"] is None or R[name][kk]["por_k"]["4"]["_e_ing"] > tg)]
            row["celdas_de_ventaja_%s" % tg] = len(adv)
            row["lista_%s" % tg] = adv
        any_adv = any_adv or row["celdas_de_ventaja_1/10"] > 0
        out["D_f1"][name] = row
    out["D_f1"]["lectura"] = "ventaja a igual coste" if any_adv else "sin ventaja a igual coste"
    for name in S1:
        if not k4_done[name]:
            out["D_f2"][name] = {"evaluada": False}
            continue
        vals = [R[name][kk]["por_k"]["4"]["_R"] for kk in keys]
        row = m8_category(vals)
        row["R4_min"] = ce.dec(min(vals))
        row["R4_max"] = ce.dec(max(vals))
        grow_k = [R[name][kk]["por_k"]["4"]["_R"] > R[name][kk]["por_k"]["3"]["_R"] for kk in keys]
        row["celdas_R4_mayor_que_R3"] = sum(grow_k)
        out["D_f2"][name] = row
    done = [n for n in S1 if k4_done[n]]
    if len(done) >= 2:
        grow = [all(R[done[j]][kk]["por_k"]["4"]["_R"] < R[done[j + 1]][kk]["por_k"]["4"]["_R"]
                    for j in range(len(done) - 1)) for kk in keys]
        out["D_f2"]["crecimiento_con_N"] = {"instancias": done, "celdas_que_crecen": sum(grow), "celdas": len(keys)}
    else:
        out["D_f2"]["crecimiento_con_N"] = {"instancias": done, "celdas_que_crecen": None, "celdas": len(keys)}
    n4 = nsil["S1_L8"]
    out["D_f3"] = {"testigo_L8_con_k4": n4["total_f4"] > 0, "n_4": n4["total_f4"], "n_4_por_na_nb": n4["f4"],
                   "celdas_con_L_ing_4_positivo": sum(1 for kk in keys
                                                      if R["S1_L8"][kk]["por_k"]["4"]["L_ing"]["exacto"] != "0/1")}
    return out


def subset_sampling(name, inst, naive_sil, kmax):
    NS = inst["NS"]
    byf = {}
    for (na, nb), c in naive_sil.items():
        byf[na + nb] = byf.get(na + nb, 0) + c
    ws = sorted(f for f, c in byf.items() if c > 0)
    if not ws:
        return {"w_estrella": None, "nota": "sin fallos silenciosos con <= %d fallos" % kmax}
    w = ws[0]
    n = byf[w]
    C = math.comb(NS, w)
    pi = F(n, C)
    nss = math.ceil(F(38416, 10000) * (1 - pi) / (F(1, 100) * pi))
    return {"w_estrella": w, "n_w_estrella": n, "C_N_w": C, "pi": float(pi), "muestras_10_por_ciento_95": nss,
            "llamadas_catalogo_k_igual_w": sum(math.comb(NS, f) for f in range(1, w + 1)),
            "nota": "solo q = p; el muestreo no da certificado ni acota la cola"}


# ------------------------------------------------------------------ programa principal
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--recursos", required=True)
    ap.add_argument("--tmp", required=True)
    ap.add_argument("--solo", choices=["PATH25"], default=None,
                    help="prueba de desarrollo: solo PATH25, con la DP en paralelo y la compuerta en k=4")
    args = ap.parse_args()
    for pth in (args.output, args.recursos):
        if os.path.exists(pth):
            raise SystemExit("La salida ya existe; use un nombre nuevo: " + pth)
    if not os.path.isdir(args.tmp):
        raise SystemExit("No existe la carpeta temporal: " + args.tmp)
    ce = load_ce()
    CE[0] = ce
    ce.LLAMADAS_MAX = LLAMADAS_MAX
    ce.SEG_MAX = SEG_MAX
    out = {"unit": "B1", "plan": "PLAN_B1.txt", "plan_sha256": SHA_PLAN,
           "programa_sha256": sha_file(os.path.abspath(__file__)), "control_exacto_sha256": SHA_CE,
           "solo": args.solo, "instancias": {}}
    status = "PASS"
    try:
        run(ce, out, args)
    except ce.Stop as ex:
        status = "FAIL"
        out["detenido_en"] = str(ex)
    out["status"] = status
    out["checks"] = ce.CHECKS
    out["llamadas_totales"] = ce.CALLS[0] + CALLS_W[0]
    with open(args.output, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(out, ensure_ascii=False, indent=1, sort_keys=True, default=str) + "\n")
    TIMES["total"] = round(time.perf_counter() - T0, 3)
    with open(args.recursos, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps({"nota": "tiempos de reloj en segundos y memoria; no deterministas; no son resultado",
                            "tiempos": TIMES, "pico_memoria_principal_bytes": ce.peak_memory_bytes(),
                            "pico_memoria_procesos_bytes": WORKER_PEAKS,
                            "llamadas_totales": ce.CALLS[0] + CALLS_W[0],
                            "llamadas_proceso_principal": ce.CALLS[0], "llamadas_procesos": CALLS_W[0]},
                           ensure_ascii=False, indent=1, default=str) + "\n")
    print("%s controles=%d llamadas=%d t=%.1fs" % (status, len(ce.CHECKS), ce.CALLS[0] + CALLS_W[0],
                                                   TIMES["total"]))
    if status != "PASS":
        print("DETENIDO EN:", out.get("detenido_en"))
        sys.exit(1)


def run(ce, out, args):
    check = ce.check
    ctx = mp.get_context("spawn")
    check("PLAN_B1.txt con su SHA-256", sha_file("B1_K4_EXACTO_S1_20261113/PLAN_B1.txt") == SHA_PLAN)
    check("FASE_A.json de D1 con su SHA-256", sha_file(ce.D1 + "FASE_A.json") == ce.SHA_FASE_A)
    check("FASE_B.json de D1 con su SHA-256", sha_file(ce.D1 + "FASE_B.json") == ce.SHA_FASE_B)
    check("ESTIMACION_H6_U4.json con su SHA-256", sha_file(ce.U4 + "ESTIMACION_H6_U4.json") == ce.SHA_U4)
    check("RESULTADO_CONTROL.json de D1-CE con su SHA-256", sha_file(CE_DIR + "RESULTADO_CONTROL.json") == SHA_RES_CE)
    check("RESULTADO_E1.json de B0 con su SHA-256", sha_file(B0_DIR + "RESULTADO_E1.json") == SHA_E1)
    A = json.load(open(ce.D1 + "FASE_A.json", encoding="utf-8"))["instancias"]
    B = json.load(open(ce.D1 + "FASE_B.json", encoding="utf-8"))["instancias"]
    U = json.load(open(ce.U4 + "ESTIMACION_H6_U4.json", encoding="utf-8"))
    RCE = json.load(open(CE_DIR + "RESULTADO_CONTROL.json", encoding="utf-8"))["instancias"]
    RE1 = json.load(open(B0_DIR + "RESULTADO_E1.json", encoding="utf-8"))["instancias"]
    sealed = {"PATH25": None, "S1_L6": RCE["S1_L6"], "S1_L7": RCE["S1_L7"], "S1_L8": RE1["S1_L8"]}
    grid = [(p, q) for p in ce.GRID_P for q in (p, p / 2, F(0))]
    keys = ["%s|%s" % (p, q) for (p, q) in grid]
    results, k4_done, nsil = {}, {}, {}
    todo = (("PATH25", 4), ("S1_L6", 6), ("S1_L7", 7), ("S1_L8", 8))
    for name, L in (todo[:1] if args.solo else todo):
        t_inst = time.perf_counter()
        rows, a = ce.path_seed(L)
        inst = ce.build(name, rows, a)
        mats = hashlib.sha256(json.dumps({"hx": inst["hx"], "hz": inst["hz"]}, sort_keys=True).encode()).hexdigest()
        if name in A:
            check("%s matrices iguales a FASE_A de D1" % name, mats == A[name]["matrices_sha256"])
        else:
            check("%s matrices iguales a E1 de B0" % name, mats == RE1[name]["matrices_sha256"])
        order = list(range(inst["NS"])) if name == "PATH25" else ce.phase_a_order(inst)
        cells = ce.U4_PAIRS if name == "PATH25" else grid
        reps = ce.logical_reps(inst)
        check("%s logicos: %d representantes = dimension logica" % (name, len(reps)), len(reps) == inst["klog"])
        klog = len(reps)
        NV = 1 + (1 << klog)
        oi = {"N": inst["N"], "M": inst["M"], "sitios": inst["NS"], "qubits_logicos": klog,
              "matrices_sha256": mats, "orden": "original" if name == "PATH25" else "fase A de D1"}
        print("%s: %d sitios" % (name, inst["NS"]), flush=True)

        # ---- catalogo
        t = time.perf_counter()
        c_main0 = ce.CALLS[0]
        POLY, naive_sil, stats, i4_bad, calls = catalogue_parallel(name, rows, a, inst, reps, KMAX, ctx)
        TIMES[name + "_catalogo"] = round(time.perf_counter() - t, 3)
        check("I4 %s: silenciosa <=> algun caracter logico vale -1 (discrepancias)" % name, i4_bad == 0, i4_bad)
        check("%s catalogo: historias = C(N,f) para f = 1..4" % name,
              all(stats[str(f)]["historias"] == math.comb(inst["NS"], f) for f in range(1, KMAX + 1)))
        n3 = sum(1 for slot in POLY.values() if any(kk[0] + kk[1] <= 3 for kk in slot))
        if name == "PATH25":
            pf = B["PATH25"]["B1"]["por_f"]
            check("PATH25 catalogo f<=3 reproduce FASE_B: historias, decoraciones y silenciosas conexas por f",
                  all(stats[f]["historias"] == pf[f]["historias"]
                      and stats[f]["decoraciones_validas"] == pf[f]["decoraciones_validas"]
                      and stats[f]["silenciosas_conexas_por_na_nb"] == pf[f]["silenciosas_por_na_nb"]
                      for f in ("1", "2", "3")))
        else:
            sc = sealed[name]["catalogo"]
            check("%s catalogo f<=3 identico al sellado (historias, decoraciones, silenciosas conexas y todas)"
                  % name, all(stats[f]["historias"] == sc["por_f"][f]["historias"]
                              and stats[f]["decoraciones_validas"] == sc["por_f"][f]["decoraciones_validas"]
                              and stats[f]["silenciosas_conexas_por_na_nb"]
                              == sc["por_f"][f]["silenciosas_conexas_por_na_nb"]
                              and stats[f]["silenciosas_todas_por_na_nb"]
                              == sc["por_f"][f]["silenciosas_todas_por_na_nb"] for f in ("1", "2", "3")))
            check("%s polimeros distintos con f<=3 = %d (sellado)" % (name, sc["polimeros_distintos"]),
                  n3 == sc["polimeros_distintos"], n3)
        if name in ("PATH25", "S1_L6"):
            t = time.perf_counter()
            ce.KMAX = KMAX
            ce.LLAMADAS_MAX = LLAMADAS_MAX - CALLS_W[0]
            P2, ns2, st2, i42 = ce.catalogue(inst, reps)
            ce.KMAX = 3
            same_st = all({kx: v for kx, v in st2[f].items() if kx != "segundos"} == stats[f] for f in stats)
            check("%s control cruzado k=4: catalogo paralelo = catalogo de D1-CE (polimeros, entradas, silenciosas "
                  "y recuentos)" % name, P2 == POLY and ns2 == naive_sil and same_st and i42 == i4_bad)
            del P2
            TIMES[name + "_catalogo_referencia"] = round(time.perf_counter() - t, 3)
        sil_f = {}
        for (na, nb), c in sorted(naive_sil.items()):
            sil_f.setdefault(str(na + nb), {})["%d,%d" % (na, nb)] = c
        nsil[name] = {"f4": sil_f.get("4", {}), "total_f4": sum(sil_f.get("4", {}).values())}
        oi["catalogo"] = {"por_f": stats, "polimeros_distintos": len(POLY), "polimeros_distintos_f_le_3": n3,
                          "llamadas": calls + (ce.CALLS[0] - c_main0),
                          "silenciosas_todas": {"%d,%d" % kk: v for kk, v in sorted(naive_sil.items())},
                          "silenciosas_por_f": sil_f}
        print("  %s catalogo listo (%.1fs); polimeros %d; silenciosas %s" %
              (name, TIMES[name + "_catalogo"], len(POLY), oi["catalogo"]["silenciosas_todas"]), flush=True)

        # ---- DP por k
        table = {}
        evaluated = {}
        gas_k = {}
        polyholder = [POLY]
        for k in range(1, KMAX + 1):
            t = time.perf_counter()
            G = prep_gas(inst, polyholder[0], order, k)
            if k == KMAX and name == "S1_L8":
                polyholder[0] = None
                POLY = None
            gate = (name == "S1_L8" or args.solo is not None) and k == KMAX
            sts, terms = explore(G, anchor_masks(G),
                                 {"states": GATE_STATES, "terms": GATE_TERMS} if gate else None)
            TIMES["%s_estructura_k%d" % (name, k)] = round(time.perf_counter() - t, 3)
            row = {"polimeros": len(G["polys"]), "celdas": {}}
            if sts is None:
                row.update({"estados": None, "terminos": None, "compuerta": dict(terms, pasa=False)})
                comp = None
            else:
                row.update({"estados": len(sts), "terminos": terms})
                print("  %s k=%d: %d polimeros, %d estados, %d terminos (%.1fs)" %
                      (name, k, len(G["polys"]), len(sts), terms, time.perf_counter() - t), flush=True)
                if k <= 3:
                    if name == "PATH25":
                        exp = {1: (37, 234, 207), 2: (387, 238, 1168), 3: (3918, 240, 7179)}[k]
                    else:
                        gk = sealed[name]["gas_por_k"][str(k)]
                        exp = (gk["polimeros"], gk["estados"], gk["terminos"])
                    check("%s k=%d DP nueva: polimeros, estados y terminos = %s (sellados)" % (name, k, list(exp)),
                          (row["polimeros"], row["estados"], row["terminos"]) == exp,
                          [row["polimeros"], row["estados"], row["terminos"]])
                nproc = 1 if name in ("PATH25", "S1_L6") and args.solo is None else W_GAS_MAX
                mem_limit = None
                if gate:
                    nproc, mem_limit, info = gate_check(G, sts, terms, NV, ctx)
                    row["compuerta"] = info
                if nproc == 0:
                    comp = None
                else:
                    t2 = time.perf_counter()
                    comp, secs = run_dp(name, k, G, sts, cells, NV, klog, nproc, ctx, args.tmp, mem_limit)
                    TIMES["%s_gas_k%d" % (name, k)] = round(time.perf_counter() - t2, 3)
                    TIMES["%s_gas_k%d_por_celda" % (name, k)] = secs
            evaluated[k] = comp is not None
            for (p, q) in cells:
                kk = "%s|%s" % (p, q)
                Li, eps_i = ce.naive(inst, p, q, naive_sil, k)
                cell = {"eps_ing": eps_i, "L_ing": Li}
                if comp is not None:
                    eps_c, Lc, P0 = comp[kk]
                    cell.update({"eps_comp": eps_c, "L_comp": Lc, "P_0_inf": P0})
                    check("I2 %s k=%d p=%s q=%s: 0 <= eps_comp <= eps_ing" % (name, k, p, q), 0 <= eps_c <= eps_i)
                    check("I5 %s k=%d p=%s q=%s: L_ing <= L_comp" % (name, k, p, q), Li <= Lc)
                    ec, ei = ce.err(Lc, eps_c), ce.err(Li, eps_i)
                    if ec is not None and ei is not None:
                        check("I6 %s k=%d p=%s q=%s: e_comp <= e_ing" % (name, k, p, q), ec <= ei)
                else:
                    cell.update({"eps_comp": None, "L_comp": None, "P_0_inf": None})
                row["celdas"][kk] = cell
            if k <= 3 and comp is not None:
                if name == "PATH25":
                    urow = U["tabla"][k - 1]
                    for (p, q) in cells:
                        kk = "%s|%s" % (p, q)
                        uc = urow["tasas"][kk]
                        check("PATH25 k=%d %s: eps_k y L de P_sil de la DP nueva identicos a U4" % (k, kk),
                              row["celdas"][kk]["eps_comp"] == frac(uc["eps_k"]["exacto"])
                              and row["celdas"][kk]["L_comp"] == frac(uc["P_sil"]["L"]["exacto"]))
                else:
                    ref = sealed[name]["celdas"]
                    same = all(row["celdas"][kk][fld] == frac(ref[kk]["por_k"][str(k)][fld]["exacto"])
                               for kk in keys for fld in ("eps_comp", "L_comp", "eps_ing", "L_ing"))
                    check("%s k=%d DP nueva y truncacion ingenua identicas a las selladas (12 celdas)" % (name, k),
                          same)
            if k == KMAX and name in ("PATH25", "S1_L6", "S1_L7") and comp is not None:
                t2 = time.perf_counter()
                POLYref = polyholder[0]
                g = ce.gas_structure(inst, POLYref, order, k)
                ok = (g["n_states"], g["terms"], len(g["polys"])) == (row["estados"], row["terminos"],
                                                                       row["polimeros"])
                for (p, q) in cells:
                    X0, XT, T = ce.gas_eval(g, inst, p, q, NV)
                    eps_c, Lc, P0 = finish_cell(p, q, X0, XT, T, klog)
                    c = row["celdas"]["%s|%s" % (p, q)]
                    ok = ok and eps_c == c["eps_comp"] and Lc == c["L_comp"] and P0 == c["P_0_inf"]
                    guard("DP de referencia de " + name)
                check("%s control cruzado k=4: DP nueva = DP de D1-CE (estados, terminos, eps_comp, L_comp, P_0; "
                      "todas las celdas)" % name, ok, [g["n_states"], g["terms"]])
                del g
                TIMES["%s_gas_k%d_referencia" % (name, k)] = round(time.perf_counter() - t2, 3)
            table[k] = row
            gas_k[str(k)] = {x: row[x] for x in ("polimeros", "estados", "terminos") if x in row}
            if "compuerta" in row:
                gas_k[str(k)]["compuerta"] = row["compuerta"]
            del G, sts
            guard("fin de k=%d en %s" % (k, name))
        for (p, q) in cells:
            kk = "%s|%s" % (p, q)
            ks = [k for k in range(1, KMAX + 1) if evaluated[k]]
            Ls = [table[k]["celdas"][kk]["L_comp"] for k in ks]
            Es = [table[k]["celdas"][kk]["eps_comp"] for k in ks]
            Li = [table[k]["celdas"][kk]["L_ing"] for k in range(1, KMAX + 1)]
            Ei = [table[k]["celdas"][kk]["eps_ing"] for k in range(1, KMAX + 1)]
            check("I3 %s %s: L_comp y L_ing no decrecen; eps_comp y eps_ing no crecen (k = 1..4)" % (name, kk),
                  all(Ls[j] <= Ls[j + 1] for j in range(len(Ls) - 1))
                  and all(Es[j] >= Es[j + 1] for j in range(len(Es) - 1))
                  and all(Li[j] <= Li[j + 1] for j in range(KMAX - 1))
                  and all(Ei[j] >= Ei[j + 1] for j in range(KMAX - 1)))
        oi["gas_por_k"] = gas_k
        res = summarize4(name, inst, table, cells, KMAX, evaluated)
        oi["celdas"] = res
        oi["muestreo_por_subconjuntos"] = subset_sampling(name, inst, naive_sil, KMAX) if name != "PATH25" else None
        results[name] = res
        k4_done[name] = evaluated[KMAX]
        oi["segundos"] = round(time.perf_counter() - t_inst, 3)
        out["instancias"][name] = oi
        del POLY, polyholder, table
        print("%-8s ok %.1fs llamadas=%d" % (name, oi["segundos"], ce.CALLS[0] + CALLS_W[0]), flush=True)
    if args.solo is not None:
        for res in results.values():
            strip_private(res)
        return
    out["decisiones"] = decide(results, keys, k4_done, nsil)
    d = out["decisiones"]
    f1 = d["D_f1"]

    def adv(n):
        return "%d/12" % f1[n]["celdas_de_ventaja_1/10"] if f1[n].get("evaluada") else "NO_EVALUADA"

    def cat(n):
        return d["D_f2"][n]["categoria"].upper().replace(" ", "_") if d["D_f2"][n].get("evaluada") else "NO_EVALUADA"

    out["lectura_global"] = ("D-f1 %s (L6 %s, L7 %s, L8 %s); D-f2 L6 %s, L7 %s, L8 %s, CRECE_CON_N %d/%d; "
                             "D-f3 TESTIGO_L8 %s (n_4 = %d)" % (
                                 f1["lectura"].upper().replace(" ", "_"), adv("S1_L6"), adv("S1_L7"), adv("S1_L8"),
                                 cat("S1_L6"), cat("S1_L7"), cat("S1_L8"),
                                 d["D_f2"]["crecimiento_con_N"]["celdas_que_crecen"] or 0,
                                 d["D_f2"]["crecimiento_con_N"]["celdas"],
                                 "SI" if d["D_f3"]["testigo_L8_con_k4"] else "NO", d["D_f3"]["n_4"]))
    for res in results.values():
        strip_private(res)


def gate_check(G, sts, terms, NV, ctx):
    """Compuerta fijada en PLAN_B1.txt para L=8 con k=4."""
    ce = CE[0]
    ns = len(sts)
    info = {"estados": ns, "terminos": terms, "limite_estados": GATE_STATES, "limite_terminos": GATE_TERMS}
    if ns > GATE_STATES or terms > GATE_TERMS:
        info.update(pasa=False, motivo="estados o terminos por encima de la compuerta")
        return 0, None, info
    A = anchor_masks(G)
    index = {J: i for i, J in enumerate(sts)}
    p = q = F(1, 10 ** 6)
    cal = eval_fly(G, A, sts, index, p, q, NV, cal_terms=min(CAL_TERMS, terms))
    t_term = cal["segundos"] / max(cal["terminos"], 1)
    del A, index
    # memoria prevista por proceso: valores de la DP con p = 1e-6 y estructura de polimeros
    bits_T = (G["N"] + G["M"]) * math.log2(10 ** 6)
    per_int = 28 + 4 * math.ceil(bits_T / 30)
    per_state = 2 * per_int + 8 + 8 + 64 + 160
    per_poly = 650
    per_proc = int(1.25 * (ns * per_state + len(G["polys"]) * per_poly) + 150 * (1 << 20))
    main_now = mem_now() or 0
    avail = MEM_TOTAL - main_now
    nproc = int(min(W_GAS_MAX, avail // per_proc)) if per_proc > 0 else W_GAS_MAX
    remaining = SEG_MAX - (time.perf_counter() - T0)
    projected = terms * t_term * 12 / nproc if nproc > 0 else None
    info.update({"calibracion": {"celda": "1/1000000|1/1000000", "terminos": cal["terminos"],
                                 "segundos": round(cal["segundos"], 3),
                                 "segundos_por_termino": t_term},
                 "memoria_prevista_por_proceso_bytes": per_proc, "memoria_principal_bytes": main_now,
                 "procesos": nproc, "tiempo_previsto_s": round(projected, 1) if projected is not None else None,
                 "tiempo_restante_s": round(remaining, 1)})
    if nproc < 1:
        info.update(pasa=False, motivo="la memoria prevista no permite ningun proceso dentro del techo")
        return 0, None, info
    if projected > 0.9 * remaining:
        info.update(pasa=False, motivo="el tiempo previsto no cabe en el techo de tiempo restante")
        return 0, None, info
    info["pasa"] = True
    mem_limit = int(avail // nproc)
    info["limite_memoria_por_proceso_bytes"] = mem_limit
    print("  compuerta L8 k=4: pasa; %d procesos; previsto %.0f s" % (nproc, projected), flush=True)
    return nproc, mem_limit, info


if __name__ == "__main__":
    main()
