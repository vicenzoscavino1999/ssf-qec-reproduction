#!/usr/bin/env python3
"""T2-0 (PLAN_T2_0.txt, sellado en SELLO_PLAN.json): calibracion del certificado por componentes con el ruido
de circuito del apendice C (una ronda, todas las ubicaciones) en las tiras HGP(camino_m, camino_6).

Reutiliza sin cambios, tras comprobar su SHA-256: b4_tiras.py de B4 (construccion de las tiras, utilidades y
techos) y, a traves de el, control_exacto.py (decoder SSF) y b1_exacto_k4.py (prep_gas, anchor_masks, explore,
eval_fly). Piezas nuevas: calendario y respuestas del circuito, propagacion independiente del marco de Pauli,
leyes de fibra (lema F), catalogo con clases, exhaustivo con clases, control de recomposicion, regla de suma con
clases (lema S) y las medidas y decisiones del plan.

Uso (desde D:\\paper ia cuantica\\qec_revision_matematica\\05_COMPROBACIONES):
  python -I -B T2_0_CALIBRACION_CIRCUITO_TIRAS_20261113/t2_calibracion.py --output <nuevo.json>
      --recursos <nuevo.json> [--prueba]
  (--prueba: prueba de desarrollo solo en m = 6 con k <= 3; no es la corrida de la unidad)
"""
import argparse
import gc
import hashlib
import json
import math
import multiprocessing as mp
import os
import random
import sys
import time
from fractions import Fraction as F
from itertools import combinations, product

sys.set_int_max_str_digits(0)
T0 = time.perf_counter()
B4DIR = "B4_TIRAS_ESCALADO_N_20261113"
UDIR = "T2_0_CALIBRACION_CIRCUITO_TIRAS_20261113/"
PLAN = UDIR + "PLAN_T2_0.txt"
SHA = {
    PLAN: "29c7c4949f314958db4f5a61d9b5e94d9a05e3c5face9359683be900312a0aa3",
    B4DIR + "/b4_tiras.py": "b1a2cb8f93affcc3713eadd4ff860a2d74756b91752011cd6b16da44def2dd15",
    B4DIR + "/RESULTADO_B4.json": "640db376441656cfb72e16aa00e0e99e2fc2eb6047905dfe031f5d994634b5f7",
}


def sha_file(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


for _p, _s in SHA.items():
    if sha_file(_p) != _s:
        raise SystemExit("SHA-256 distinto del sellado: " + _p)
sys.path.insert(0, os.path.abspath(B4DIR))
import b4_tiras as b4  # noqa: E402  (tras comprobar su SHA-256)

D = 10 ** 14
PX = F(1, D + 1)
GIB = 1 << 30
SEG_MAX = 3 * 3600.0
CALLS_MAX = 1_200_000_000
MEM_TOTAL = 18 * GIB
MEM_MAIN = 12 * GIB
MEM_W = GIB
PROCS = 14
GATE_STATES = 5_000_000
GATE_TERMS = 2_500_000_000
REC_TARGET = 10000
REC_TRIES = 2_000_000
RESP_RANDOM = 2000
SEED = 20261113
LIMITS = []
b4.T0, b4.SEG_MAX, b4.MEM_MAIN, b4.LLAMADAS_MAX = T0, SEG_MAX, MEM_MAIN, CALLS_MAX
bits = b4.bits
check = b4.check


def guard(where):
    if time.perf_counter() - T0 > SEG_MAX:
        raise b4.Stop("techo de tiempo alcanzado en " + where)
    m = b4.mem_now()
    if m is not None and m > MEM_MAIN:
        raise b4.Stop("techo de memoria del proceso principal alcanzado en " + where)
    if b4.CALLS[0] > CALLS_MAX:
        raise b4.Stop("techo de llamadas alcanzado en " + where)


# ------------------------------------------------------------------ circuito del apendice C
def schedule(inst):
    """Calendario voraz: incidencias (c, v) en orden lexicografico, menor color positivo libre."""
    hx, N, M, cols = inst["hx"], inst["N"], inst["M"], inst["cols"]
    used_c = [set() for _ in range(M)]
    used_v = [set() for _ in range(N)]
    layer = {}
    for c in range(M):
        for v in bits(hx[c]):
            t = 1
            while t in used_c[c] or t in used_v[v]:
                t += 1
            layer[(c, v)] = t
            used_c[c].add(t)
            used_v[v].add(t)
    L = max(x.bit_count() for x in hx) + max(cols[v].bit_count() for v in range(N)) - 1
    assert max(layer.values()) <= L
    return layer, L


def locations(inst, layer, L):
    N, M = inst["N"], inst["M"]
    by_v = [dict() for _ in range(N)]
    by_c = [dict() for _ in range(M)]
    for (c, v), t in layer.items():
        by_v[v][t] = c
        by_c[c][t] = v
    locs = []
    for v in range(N):
        locs.append(("dwait", v, 0, None))
        for t in range(1, L + 1):
            locs.append(("cnot", v, t, by_v[v][t]) if t in by_v[v] else ("dwait", v, t, None))
        locs.append(("dwait", v, L + 1, None))
    for c in range(M):
        locs.append(("prep", c, 0, None))
        for t in range(1, L + 1):
            if t not in by_c[c]:
                locs.append(("await", c, t, None))
        locs.append(("read", c, L + 1, None))
    return locs


VISIBLE = {"dwait": [(1, 0)], "prep": [(0, 1)], "await": [(0, 1)], "read": [(0, 1)],
           "cnot": [(1, 0), (0, 1), (1, 1)]}


def tmask(inst, layer, v, t):
    out = 0
    for c in bits(inst["cols"][v]):
        if layer[(c, v)] > t:
            out |= 1 << c
    return out


def resp_formula(inst, layer, loc, cat):
    """eq:circ-response: (e, s) de la categoria visible cat = (a, b) en la ubicacion loc."""
    kind, base, t, c = loc
    a, b = cat
    if kind == "dwait":
        return (1 << base, tmask(inst, layer, base, t))
    if kind == "cnot":
        return ((1 << base) if a else 0, (tmask(inst, layer, base, t) if a else 0) ^ ((1 << c) if b else 0))
    return (0, 1 << base)


def propagate(inst, layer, L, faults):
    """Propagacion independiente del marco X: CNOT dato -> ancilla capa a capa; fallos tras su operacion."""
    by_t = {}
    for (c, v), t in layer.items():
        by_t.setdefault(t, []).append((c, v))
    inj = {}
    for loc, cat in faults:
        inj.setdefault(loc[2], []).append((loc, cat))
    st = [0, 0, 0]          # X de datos, X de ancillas, flips clasicos de lectura

    def inject(t):
        for (kind, base, tt, c), (a, b) in inj.get(t, []):
            if kind == "dwait":
                st[0] ^= 1 << base
            elif kind == "cnot":
                if a:
                    st[0] ^= 1 << base
                if b:
                    st[1] ^= 1 << c
            elif kind in ("prep", "await"):
                st[1] ^= 1 << base
            elif kind == "read":
                st[2] ^= 1 << base

    inject(0)
    for t in range(1, L + 1):
        for (c, v) in by_t.get(t, []):
            if st[0] >> v & 1:
                st[1] ^= 1 << c
        inject(t)
    inject(L + 1)
    return st[0], st[1] ^ st[2]


def fibres(inst, locs):
    N = inst["N"]
    fd, fc = [[] for _ in range(N)], [[] for _ in range(inst["M"])]
    for loc in locs:
        (fd if loc[0] in ("dwait", "cnot") else fc)[loc[1]].append(loc)
    return fd, fc


def site_classes(inst, layer, L, locs, restricted):
    """Lema F: clases no nulas alcanzables de cada sitio (respuesta global y clave canonica) y tipo de ley."""
    N, M, cols = inst["N"], inst["M"], inst["cols"]
    fd, fc = fibres(inst, locs)
    CL, CK, LTK = [], [], []
    for v in range(N):
        chk = sorted(bits(cols[v]), key=lambda c: layer[(c, v)])
        lays = tuple(layer[(c, v)] for c in chk)
        LTK.append(("d", lays))
        canon = (lambda e, s, chk=chk: (1 if e else 0,) + tuple(s >> c & 1 for c in chk))
        if restricted:
            CL.append([(1 << v, cols[v])])
            CK.append([canon(1, cols[v])])
            continue
        reach = {(0, 0)}
        for loc in fd[v]:
            rs = [resp_formula(inst, layer, loc, cat) for cat in VISIBLE[loc[0]]]
            reach = reach | {(e ^ e2, s ^ s2) for (e, s) in reach for (e2, s2) in rs}
        cl = sorted(reach - {(0, 0)}, key=lambda es: canon(*es))
        CL.append(cl)
        CK.append([canon(*es) for es in cl])
    for c in range(M):
        LTK.append(("c", sum(1 for loc in fc[c] if loc[0] == "await")))
        CL.append([(0, 1 << c)])
        CK.append([(1,)])
    return CL, CK, LTK


# ------------------------------------------------------------------ leyes exactas (subfamilia concreta)
def conv(dist, opts):
    new = {}
    for k1, p1 in dist.items():
        for k2, p2 in opts:
            k = tuple(x ^ y for x, y in zip(k1, k2))
            new[k] = new.get(k, 0) + p1 * p2
    return new


def law_of(lt, L, rc, rm):
    w = 2 * rc / 3
    if lt[0] == "c":
        dist = {(0,): F(1)}
        for _ in range(1 + lt[1]):
            dist = conv(dist, [((0,), 1 - w), ((1,), w)])
        fr = w * (1 - rm) + (1 - w) * rm
        return conv(dist, [((0,), 1 - fr), ((1,), fr)])
    lays = lt[1]
    d = len(lays)
    zero = (0,) * (d + 1)

    def T(t):
        return tuple(1 if lays[i] > t else 0 for i in range(d))

    dist = {zero: F(1)}
    for t in range(L + 2):
        if t in lays:
            j = lays.index(t)
            a10 = (1,) + T(t)
            b01 = (0,) + tuple(1 if i == j else 0 for i in range(d))
            a11 = tuple(x ^ y for x, y in zip(a10, b01))
            dist = conv(dist, [(zero, 1 - 4 * rc / 5), (a10, 4 * rc / 15), (b01, 4 * rc / 15),
                               (a11, 4 * rc / 15)])
        else:
            dist = conv(dist, [(zero, 1 - w), ((1,) + T(t), w)])
    return dist


def location_law_sums(rc, rm):
    """Las leyes de ubicacion completas (con categorias invisibles) suman 1."""
    one_q = (1 - rc) + 3 * (rc / 3)
    cnot = (1 - rc) + 15 * (rc / 15)
    cnot_x = (1 - rc + 3 * rc / 15) + 3 * (4 * rc / 15)
    read = one_q * ((1 - rm) + rm)
    return [one_q, cnot, cnot_x, read]


# ------------------------------------------------------------------ procesos de trabajo
_W = {}


def w_init(m, restricted):
    sys.set_int_max_str_digits(0)
    ce, rect, _ = b4.load_all(with_b1=False)
    inst = b4.make_inst(ce, rect, "tira", m)
    dec, RES = ce.make_decoder(inst)
    layer, L = schedule(inst)
    locs = locations(inst, layer, L)
    CL, CK, LTK = site_classes(inst, layer, L, locs, restricted)
    lts = sorted(set(LTK))
    _W.update(ce=ce, inst=inst, dec=dec, RES=RES, group=b4.group_of(ce.logical_reps(inst)), CL=CL,
              LT=[lts.index(x) for x in LTK], CONN={})


def connected(Y, adj, CONN):
    ok = CONN.get(Y)
    if ok is None:
        seen = 0
        todo = Y & -Y
        while todo:
            bb = todo & -todo
            todo ^= bb
            seen |= bb
            todo |= adj[bb.bit_length() - 1] & Y & ~seen
        ok = seen == Y
        if len(CONN) > 400000:
            CONN.clear()
        CONN[Y] = ok
    return ok


def classify(sm, sites, idxs):
    W = _W
    ce, inst, CL = W["ce"], W["inst"], W["CL"]
    E = s = 0
    for b, i in zip(sites, idxs):
        e, x = CL[b][i]
        E ^= e
        s ^= x
    if len(W["RES"]) > 300000:
        W["RES"].clear()
    cor, touched = W["dec"](s)
    Y = sm | touched
    r = E ^ cor
    passes = ce.lin(inst["cols"], r) == 0
    silent, chis, bad = False, None, 0
    if passes:
        silent = not ce.inside(r, inst["bz"])
        chis = tuple(-1 if (u & r).bit_count() & 1 else 1 for u in W["group"])
        bad = int(silent != any(c < 0 for c in chis))
    return Y, passes, silent, chis, bad


def encode(sm, idxs, NS):
    """Decoracion (sitios, clases) como un entero: mascara de sitios y 3 bits por clase (clase <= 6)."""
    out, sh = sm, NS
    for i in idxs:
        out |= i << sh
        sh += 3
    return out


def w_grow(task):
    f, chunk, want_rec, keep_codes, nchar = task
    W = _W
    inst, CL, LT = W["inst"], W["CL"], W["LT"]
    adj, DATA = inst["adjK"], (1 << inst["N"]) - 1
    poly, rec, sil, codes = {}, {}, {}, []
    calls = bad = nval = 0
    for sm in chunk:
        sites = list(bits(sm))
        na = (sm & DATA).bit_count()
        key = (na, len(sites) - na)
        yu = 0
        for idxs in product(*[range(len(CL[b])) for b in sites]):
            Y, passes, silent, chis, b_ = classify(sm, sites, idxs)
            calls += 1
            bad += b_
            if not connected(Y, adj, W["CONN"]):
                continue
            nval += 1
            yu |= Y
            slot = poly.get(Y)
            if slot is None:
                slot = {}
                poly[Y] = slot
            ent = slot.get(key)
            if ent is None:
                ent = [0] * nchar
                slot[key] = ent
            ent[0] += 1
            if passes:
                ent[1] += 1
                for j, c in enumerate(chis):
                    ent[2 + j] += c
            if silent:
                sig = tuple(sorted((LT[b], i) for b, i in zip(sites, idxs)))
                sil[sig] = sil.get(sig, 0) + 1
            if keep_codes:
                codes.append(encode(sm, idxs, inst["NS"]))
        if yu and want_rec:
            rec[sm] = yu
    W["RES"].clear()
    m = b4.mem_now()
    if m is not None and m > MEM_W:
        raise RuntimeError("techo de memoria de un proceso del catalogo (%d bytes)" % m)
    return {"poly": poly, "rec": rec, "sil": sil, "codes": codes, "calls": calls, "bad": bad, "nval": nval,
            "peak": b4.mem_now(peak=True), "pid": os.getpid()}


def w_exh(task):
    f, firsts = task
    W = _W
    inst, CL = W["inst"], W["CL"]
    NS, adj = inst["NS"], inst["adjK"]
    codes, hist, silh, bad = [], 0, 0, 0
    for first in firsts:
        for rest in combinations(range(first + 1, NS), f - 1):
            sites = [first] + list(rest)
            sm = 0
            for i in sites:
                sm |= 1 << i
            for idxs in product(*[range(len(CL[b])) for b in sites]):
                Y, passes, silent, chis, b_ = classify(sm, sites, idxs)
                hist += 1
                bad += b_
                silh += int(silent)
                if connected(Y, adj, W["CONN"]):
                    codes.append(encode(sm, idxs, NS))
        m = b4.mem_now()
        if m is not None and m > MEM_W:
            raise RuntimeError("techo de memoria de un proceso del exhaustivo (%d bytes)" % m)
    W["RES"].clear()
    return {"codes": codes, "hist": hist, "silh": silh, "bad": bad, "peak": b4.mem_now(peak=True),
            "pid": os.getpid()}


# ------------------------------------------------------------------ catalogo con clases (regla de B2)
def grow(pool, inst, name, nclass, kmax, keep_codes_upto, nchar):
    NS, adj = inst["NS"], inst["adjK"]
    rec_prev, POLY, sil, info, codes = None, {}, {}, {}, {}
    bad = 0
    for f in range(1, kmax + 1):
        t = time.perf_counter()
        if f == 1:
            cand = [1 << i for i in range(NS)]
        else:
            cs = set()
            for A, yu in rec_prev.items():
                nb = yu
                for v in bits(yu):
                    nb |= adj[v]
                nb &= ~A
                for x in bits(nb):
                    cs.add(A | (1 << x))
            cand = sorted(cs)
            del cs
        wts = [math.prod(nclass[b] for b in bits(sm)) for sm in cand]
        need = sum(wts)
        if b4.CALLS[0] + need > CALLS_MAX:
            raise b4.Stop("techo de llamadas: catalogo de %s con f=%d (%d asignaciones) no cabe" % (name, f, need))
        target = max(1, need // (8 * PROCS))
        chunks, cur, acc = [], [], 0
        for sm, w in zip(cand, wts):
            cur.append(sm)
            acc += w
            if acc >= target:
                chunks.append(cur)
                cur, acc = [], 0
        if cur:
            chunks.append(cur)
        want_rec = f < kmax
        keep = f <= keep_codes_upto
        rec_new, sf, cf = {}, {}, []
        nval = calls = 0
        try:
            for res in pool.imap_unordered(w_grow, [(f, ch, want_rec, keep, nchar) for ch in chunks]):
                calls += res["calls"]
                bad += res["bad"]
                nval += res["nval"]
                b4.note_peak("catalogo_" + name, res)
                for Y, slot in res["poly"].items():
                    tg = POLY.get(Y)
                    if tg is None:
                        POLY[Y] = slot
                        continue
                    for kk, ent in slot.items():
                        e0 = tg.get(kk)
                        if e0 is None:
                            tg[kk] = ent
                        else:
                            for j in range(nchar):
                                e0[j] += ent[j]
                for A, yu in res["rec"].items():
                    rec_new[A] = rec_new.get(A, 0) | yu
                for sg, c in res["sil"].items():
                    sf[sg] = sf.get(sg, 0) + c
                if keep:
                    cf.extend(res["codes"])
                res.clear()
                guard("catalogo de %s f=%d" % (name, f))
        except (b4.Stop, b4.CheckFail):
            raise
        except Exception as ex:
            raise b4.Stop("catalogo de %s f=%d: %s" % (name, f, ex))
        b4.CALLS[0] += calls
        check("%s f=%d: cada asignacion de cada candidato se decodifica una vez (%d)" % (name, f, need),
              calls == need, verbose=False)
        sil[f] = sf
        if keep:
            codes[f] = cf
        info[str(f)] = {"candidatos": len(cand), "asignaciones": need, "validas": nval,
                        "conjuntos_validos": len(rec_new) if want_rec else None,
                        "segundos": round(time.perf_counter() - t, 3)}
        print("  %s f=%d: %d candidatos, %d asignaciones, %d validas (%.1fs)" % (
            name, f, len(cand), need, nval, time.perf_counter() - t), flush=True)
        rec_prev = rec_new
        del cand, chunks, wts
    check("I4 %s: silenciosa <=> algun caracter logico vale -1 (discrepancias = 0)" % name, bad == 0, bad)
    return POLY, sil, info, codes


def exhaustive(pool, inst, name, kmax):
    NS = inst["NS"]
    out = {}
    bad = 0
    for f in range(1, kmax + 1):
        t = time.perf_counter()
        codes, hist, silh = [], 0, 0
        try:
            for res in pool.imap_unordered(w_exh, b4.split_firsts(NS, f, 4 * PROCS)):
                codes.extend(res["codes"])
                hist += res["hist"]
                silh += res["silh"]
                bad += res["bad"]
                b4.note_peak("exhaustivo_" + name, res)
                guard("exhaustivo de %s f=%d" % (name, f))
        except (b4.Stop, b4.CheckFail):
            raise
        except Exception as ex:
            raise b4.Stop("exhaustivo de %s f=%d: %s" % (name, f, ex))
        b4.CALLS[0] += hist
        out[f] = {"codes": codes, "historias": hist, "silenciosas": silh}
        print("  %s exhaustivo f=%d: %d historias, %d validas, %d silenciosas (%.1fs)" % (
            name, f, hist, len(codes), silh, time.perf_counter() - t), flush=True)
    check("I4 %s exhaustivo: silenciosa <=> algun caracter logico vale -1" % name, bad == 0, bad)
    return out


# ------------------------------------------------------------------ DP, regla de suma y negativo
def structure(b1, inst, POLY, kmax, label, keep_k):
    coord, NS = inst["coord"], inst["NS"]
    order = sorted(range(NS), key=lambda v: (coord[v][0], coord[v][1], v))
    rows, keep = {}, None
    for k in range(1, kmax + 1):
        t = time.perf_counter()
        G = b1.prep_gas(inst, POLY, order, k)
        A = b1.anchor_masks(G)
        sts, terms = b1.explore(G, A)
        rows[str(k)] = {"polimeros": len(G["polys"]), "estados": len(sts), "terminos": terms,
                        "segundos": round(time.perf_counter() - t, 3)}
        print("  %s k=%d: %d polimeros, %d estados, %d terminos (%.1fs)" % (
            label, k, len(G["polys"]), len(sts), terms, time.perf_counter() - t), flush=True)
        if k == keep_k:
            keep = (G, A, sts)
        else:
            del G, A, sts
        guard("estructura de %s k=%d" % (label, k))
    return rows, keep


def esym(nlist, kmax):
    e = [1] + [0] * kmax
    for nb in nlist:
        for j in range(kmax, 0, -1):
            e[j] += nb * e[j - 1]
    return e


def criterion(nlist, k, X0):
    NS = len(nlist)
    e = esym(nlist, k + 1)
    T = math.prod(D + nb for nb in nlist)
    thr = D ** (NS - k)
    tail = T - sum(e[j] * D ** (NS - j) for j in range(k + 1))
    pre = tail < thr
    S = T - X0
    check("regla de suma NS=%d k=%d: 0 <= T_cat - X0" % (NS, k), S >= 0, verbose=False)
    below = S < thr
    ver = ("COMPLETO" if below else "INCOMPLETO") if pre else "NO_DECIDIBLE"
    return {"k": k, "condicion_previa": pre, "cola_sobre_umbral": float(F(tail, thr)), "veredicto": ver,
            "margen_log10": (math.log10(thr) - math.log10(S)) if S > 0 else None,
            "peso_orden_siguiente_aprox": float(F(S, D ** (NS - k - 1)))}


def sum_rule(b1, keep, nlist, k, label):
    t = time.perf_counter()
    G, A, sts = keep
    index = {J: i for i, J in enumerate(sts)}
    X0, _, T = b1.eval_fly(G, A, sts, index, PX, PX, 1, mem_limit=MEM_MAIN)
    check("regla de suma %s k=%d: normalizacion de eval_fly T = (D+1)^NS, luego X0 = D^NS Xi(P)" % (label, k),
          T == (D + 1) ** len(nlist), verbose=False)
    cr = criterion(nlist, k, X0)
    cr["segundos"] = round(time.perf_counter() - t, 3)
    print("  regla de suma %s k=%d: %s (margen %s ordenes; %.1fs)" % (
        label, k, cr["veredicto"], "%.2f" % cr["margen_log10"] if cr["margen_log10"] is not None else "-",
        time.perf_counter() - t), flush=True)
    guard("regla de suma de " + label)
    return cr


def remove_one(POLY, d):
    Ys = [Y for Y, slot in POLY.items() if any(a + b == d for (a, b) in slot)]
    Y = min(Ys)
    key = min(kk for kk in POLY[Y] if kk[0] + kk[1] == d)
    new = dict(POLY)
    slot = {kk: list(v) for kk, v in POLY[Y].items()}
    before = slot[key][0]
    slot[key][0] -= 1
    if slot[key][0] == 0:
        del slot[key]
    if slot:
        new[Y] = slot
    else:
        del new[Y]
    return new, {"d": d, "soporte_hex": hex(Y), "clave": list(key), "recuento_antes": before}


# ------------------------------------------------------------------ controles del circuito y de las leyes
def control_circuit(inst, name, grid_rates):
    layer, L = schedule(inst)
    locs = locations(inst, layer, L)
    N, M = inst["N"], inst["M"]
    hX = sum(x.bit_count() for x in inst["hx"])
    Nloc = (L + 2) * (N + M) - hX
    check("C-resp %s: N_loc = (L+2)(n+m_X) - h_X = %d" % (name, Nloc), len(locs) == Nloc, len(locs))
    mism = 0
    nsingle = 0
    for loc in locs:
        for cat in VISIBLE[loc[0]]:
            nsingle += 1
            if resp_formula(inst, layer, loc, cat) != propagate(inst, layer, L, [(loc, cat)]):
                mism += 1
    check("C-resp %s: respuesta por formula = propagacion en las %d ubicaciones y categorias visibles"
          % (name, nsingle), mism == 0, mism)
    rng = random.Random(SEED + inst["NS"])
    mism = 0
    for _ in range(RESP_RANDOM):
        ls = rng.sample(locs, rng.randint(2, 20))
        faults = [(loc, rng.choice(VISIBLE[loc[0]])) for loc in ls]
        e = s = 0
        for loc, cat in faults:
            e2, s2 = resp_formula(inst, layer, loc, cat)
            e ^= e2
            s ^= s2
        if (e, s) != propagate(inst, layer, L, faults):
            mism += 1
    check("C-resp %s: XOR de respuestas = propagacion conjunta en %d asignaciones aleatorias"
          % (name, RESP_RANDOM), mism == 0, mism)
    CL, CK, LTK = site_classes(inst, layer, L, locs, False)
    cnt = {}
    for v in range(N):
        d = inst["cols"][v].bit_count()
        cnt.setdefault(d, set()).add(len(CL[v]))
    check("C-ley %s: clases no nulas por sitio de dato 7 (d = 2) y 3 (d = 1); check 1" % name,
          cnt.get(2, {7}) == {7} and cnt.get(1, {3}) == {3} and all(len(CL[N + c]) == 1 for c in range(M)),
          {str(d): sorted(v) for d, v in cnt.items()})
    lts = sorted(set(LTK))
    rep = {lt: LTK.index(lt) for lt in lts}
    ok_sum = ok_sup = True
    for (rc, rm) in grid_rates:
        if any(x != 1 for x in location_law_sums(rc, rm)):
            ok_sum = False
        for lt in lts:
            law = law_of(lt, L, rc, rm)
            if sum(law.values()) != 1:
                ok_sum = False
            sup = {k for k, p in law.items() if p > 0} - {tuple([0] * len(next(iter(law))))}
            if sup != set(CK[rep[lt]]):
                ok_sup = False
    check("C-ley %s: leyes de ubicacion y de fibra suman 1 en las %d celdas" % (name, len(grid_rates)), ok_sum)
    check("C-ley %s: soporte de cada ley de fibra = clases alcanzables (%d tipos de ley)" % (name, len(lts)),
          ok_sup)
    deg = max(x.bit_count() for x in inst["adjK"])
    return {"capas": L, "colores_usados": max(layer.values()), "h_X": hX, "ubicaciones": Nloc,
            "host_primitivo": Nloc + N, "tipos_de_ley": len(lts), "grado_max_K": deg,
            "clases_por_sitio": {"dato_d2": 7, "dato_d1": 3, "check": 1}}, (layer, L, locs, CL, CK, LTK)


def control_recomp(ce, inst, CL, name, target):
    dec, RES = ce.make_decoder(inst)
    NS, N, adj = inst["NS"], inst["N"], inst["adjK"]
    DATA = (1 << N) - 1
    rng = random.Random(SEED + 7 * NS)
    done = tries = 0
    fails = []
    hist = {}
    c0 = ce.CALLS[0]

    def run(sites, idxs):
        E = s = 0
        for b, i in zip(sites, idxs):
            e, x = CL[b][i]
            E ^= e
            s ^= x
        cor, touched = dec(s)
        return E, cor, touched

    while done < target and tries < REC_TRIES:
        tries += 1
        sites = sorted(rng.sample(range(NS), rng.randint(4, 8)))
        idxs = [rng.randrange(len(CL[b])) for b in sites]
        sm = sum(1 << b for b in sites)
        E, cor, touched = run(sites, idxs)
        Y = sm | touched
        comps = b4_components(Y, adj)
        parts = []
        for cmask in comps:
            sub = [(b, i) for b, i in zip(sites, idxs) if cmask >> b & 1]
            Ec, corc, touc = run([b for b, _ in sub], [i for _, i in sub])
            parts.append((cmask, sum(1 << b for b, _ in sub), Ec, corc, touc))
        if sum(1 for p in parts if p[4]) < 2:
            continue
        done += 1
        hist[str(len(comps))] = hist.get(str(len(comps)), 0) + 1
        ok = all((Sc | touc) == cmask for cmask, Sc, Ec, corc, touc in parts)
        cx = tx = rx = 0
        for cmask, Sc, Ec, corc, touc in parts:
            cx ^= corc
            tx |= touc
            rx ^= Ec ^ corc
        if not ok or cx != cor or tx != touched or rx != (E ^ cor) or (E & ~DATA):
            fails.append({"sitios": sites, "clases": idxs})
            break
        if len(RES) > 300000:
            RES.clear()
    b4.CALLS[0] += ce.CALLS[0] - c0
    check("C-rec %s: %d historias con >= 2 componentes con acciones; ejecucion global = suma de las "
          "componentes" % (name, done), not fails and done >= target, {"intentos": tries, "fallos": fails[:1]})
    return {"historias_multicomponente": done, "intentos": tries, "componentes_por_historia": hist,
            "llamadas": ce.CALLS[0] - c0}


def b4_components(Y, adj):
    out = []
    rest = Y
    while rest:
        seen = 0
        todo = rest & -rest
        while todo:
            bb = todo & -todo
            todo ^= bb
            seen |= bb
            todo |= adj[bb.bit_length() - 1] & Y & ~seen
        out.append(seen)
        rest &= ~seen
    return out


# ------------------------------------------------------------------ metricas informativas
def leading_psil(sil, LTK, L, cells, NS):
    fs = min((f for f in sil if sil[f]), default=None)
    if fs is None:
        return fs, None
    lts = sorted(set(LTK))
    out = {}
    for (rc, rm) in cells:
        laws = [law_of(lt, L, rc, rm) for lt in lts]
        zero = [law[tuple([0] * len(next(iter(law))))] for law in laws]
        classes = [sorted(k for k in law if any(k)) for law in laws]
        odds = [[law[c] / z for c in cl] for law, cl, z in zip(laws, classes, zero)]
        Bz = math.prod(zero[lts.index(lt)] for lt in LTK)
        tot = sum(cnt * math.prod(odds[lt][i] for lt, i in sg) for sg, cnt in sil[fs].items())
        P = Bz * tot
        out["%s|%s" % (rc, rm)] = {"P_sil_dominante": float(P), "montecarlo_10pc_95": float(384 / P) if P else None}
    return fs, out


# ------------------------------------------------------------------ una tira
def run_strip(ctx, ce, rect, b1, m, kmax, R, args, I, cells_rates):
    name = "TIRA_%dx6" % m
    print("== %s (k <= %d)" % (name, kmax), flush=True)
    t_inst = time.perf_counter()
    inst = b4.make_inst(ce, rect, "tira", m)
    check("%s: matrices iguales a las de B4" % name, b4.mats_sha(inst) == R[name]["matrices_sha256"])
    circ, (layer, L, locs, CL, CK, LTK) = control_circuit(inst, name, cells_rates)
    nclass = [len(c) for c in CL]
    klog = len(ce.logical_reps(inst))
    nchar = 1 + (1 << klog)
    o = {"sitios": inst["NS"], "circuito": circ}
    keep_codes = 3 if (m == 6 and not args.no_exh) else 0
    c0 = b4.CALLS[0]
    with ctx.Pool(PROCS, initializer=w_init, initargs=(m, False)) as pool:
        POLY, sil, cinfo, codes = grow(pool, inst, name, nclass, kmax, keep_codes, nchar)
        o["catalogo"] = {"por_f": cinfo, "polimeros_distintos": len(POLY), "llamadas": b4.CALLS[0] - c0}
        if keep_codes:
            ex = exhaustive(pool, inst, name, 3)
            for f in range(1, 4):
                same = set(ex[f]["codes"]) == set(codes[f])
                check("C-exh %s f=%d: decoraciones validas del exhaustivo = catalogo por crecimiento (%d)"
                      % (name, f, len(ex[f]["codes"])), same and len(ex[f]["codes"]) == len(codes[f]))
            o["exhaustivo"] = {str(f): {"historias": ex[f]["historias"], "validas": len(ex[f]["codes"]),
                                        "silenciosas": ex[f]["silenciosas"]} for f in ex}
            fs_exh = min((f for f in ex if ex[f]["silenciosas"]), default=None)
            o["f_s_exhaustivo"] = fs_exh
            del ex, codes
    gc.collect()
    rows, keep = structure(b1, inst, POLY, kmax, name, kmax)
    npoly = len(POLY)
    if m != 6:
        del POLY
        gc.collect()
    o["estructura"] = rows
    o["cocientes_frente_a_B4"] = ratios_b4(o, R[name], kmax)
    o["regla_de_suma"] = sum_rule(b1, keep, nclass, kmax, name)
    check("C-suma %s k=%d: el criterio dice COMPLETO" % (name, kmax), o["regla_de_suma"]["veredicto"] == "COMPLETO",
          o["regla_de_suma"])
    del keep
    gc.collect()
    if m == 6:
        P2, info = remove_one(POLY, kmax)
        _, keep2 = structure(b1, inst, P2, kmax, name + "_menos_%d" % kmax, kmax)
        info["criterio"] = sum_rule(b1, keep2, nclass, kmax, name + "_menos_%d" % kmax)
        check("C-suma negativo %s: sin una decoracion de %d sitios el criterio dice INCOMPLETO" % (name, kmax),
              info["criterio"]["veredicto"] == "INCOMPLETO", info["criterio"])
        o["regla_de_suma_negativo"] = info
        del P2, keep2, POLY
    if m in (6, 12):
        o["recomposicion"] = control_recomp(ce, inst, CL, name, REC_TARGET if not args.prueba else 1000)
    fs, lead = leading_psil(sil, LTK, L, cells_rates, inst["NS"])
    o["f_s_circuito"] = fs
    o["silenciosas_por_f"] = {str(f): sum(v.values()) for f, v in sil.items()}
    o["silenciosas_f_s_por_firma"] = len(sil[fs]) if fs else 0
    o["P_sil_dominante"] = lead
    e = esym(nclass, kmax)
    o["truncacion_ingenua_decodificaciones"] = {str(k): sum(e[1:k + 1]) for k in range(1, kmax + 1)}
    o["segundos"] = round(time.perf_counter() - t_inst, 3)
    I[name] = o
    print("  %s: f_s = %s; %d polimeros; %.1fs" % (name, fs, npoly, o["segundos"]), flush=True)
    del sil
    gc.collect()


def control_fen(ctx, ce, rect, b1, R, I):
    name = "TIRA_6x6"
    print("== C-fen: %s con las clases fenomenologicas" % name, flush=True)
    inst = b4.make_inst(ce, rect, "tira", 6)
    layer, L = schedule(inst)
    locs = locations(inst, layer, L)
    CL, CK, LTK = site_classes(inst, layer, L, locs, True)
    nchar = 1 + (1 << len(ce.logical_reps(inst)))
    with ctx.Pool(PROCS, initializer=w_init, initargs=(6, True)) as pool:
        POLY, sil, cinfo, _ = grow(pool, inst, name + "_fen", [1] * inst["NS"], 4, 0, nchar)
    rc = R[name]["catalogo"]
    check("C-fen: validas por f iguales a B4 %s" % [rc["por_f"][str(f)]["validas"] for f in range(1, 5)],
          all(cinfo[str(f)]["validas"] == rc["por_f"][str(f)]["validas"] for f in range(1, 5)))
    check("C-fen: polimeros distintos = %d (B4)" % rc["polimeros_distintos"], len(POLY) == rc["polimeros_distintos"])
    rows, keep = structure(b1, inst, POLY, 4, name + "_fen", 4)
    gk = R[name]["gas_por_k"]["4"]
    got = [rows["4"]["polimeros"], rows["4"]["estados"], rows["4"]["terminos"]]
    check("C-fen: estructura de la DP con k=4 = %s (B4)" % [gk["polimeros"], gk["estados"], gk["terminos"]],
          got == [gk["polimeros"], gk["estados"], gk["terminos"]], got)
    G, A, sts = keep
    index = {J: i for i, J in enumerate(sts)}
    grid = [(p, q) for p in ce.GRID_P for q in (p, p / 2, F(0))]
    same = 0
    for (p, q) in grid:
        X0, _, T = b1.eval_fly(G, A, sts, index, p, q, 1)
        kk = "%s|%s" % (p, q)
        same += int(1 - F(X0, T) == b4.frac(R[name]["celdas"][kk]["por_k"]["4"]["eps_comp"]["exacto"]))
    check("C-fen: eps_comp(4) identico al sellado en B4 en las 12 celdas", same == 12, same)
    I["C_fen"] = {"catalogo_por_f": {f: cinfo[f]["validas"] for f in cinfo}, "polimeros": len(POLY),
                  "estructura_k4": got, "celdas_identicas": same}
    del POLY, keep, G, A, sts
    gc.collect()


def ratios_b4(o, Rn, kmax):
    out = {}
    for k in range(1, kmax + 1):
        g, mine = Rn["gas_por_k"][str(k)], o["estructura"][str(k)]
        calls = sum(o["catalogo"]["por_f"][str(f)]["asignaciones"] for f in range(1, k + 1))
        out[str(k)] = {"llamadas": calls / Rn["llamadas_por_k"][str(k)]["catalogo_por_crecimiento"],
                       "polimeros": mine["polimeros"] / g["polimeros"], "estados": mine["estados"] / g["estados"],
                       "terminos": mine["terminos"] / g["terminos"]}
    return out


def project(I, out):
    """D-c1 con las reglas fijadas en el plan."""
    def val(m, k, what):
        o = I.get("TIRA_%dx6" % m)
        if not o or str(k) not in o.get("estructura", {}):
            return None
        if what == "llamadas":
            return sum(o["catalogo"]["por_f"][str(f)]["asignaciones"] for f in range(1, k + 1))
        return o["estructura"][str(k)][what]

    res = {}
    for what in ("estados", "terminos", "llamadas"):
        v6, v12 = val(6, 4, what), val(12, 4, what)
        v24_3, v12_3 = val(24, 3, what), val(12, 3, what)
        p24 = (v24_3 * v12 / v12_3) if None not in (v24_3, v12, v12_3) and v12_3 else None
        p48 = (3 * p24 - 2 * v12) if None not in (p24, v12) else None
        res[what] = {"6": v6, "12": v12, "24_proyectado": p24, "48_proyectado": p48}
    dec = {}
    for m, key in (("6", "6"), ("12", "12"), ("24", "24_proyectado"), ("48", "48_proyectado")):
        st, te = res["estados"][key], res["terminos"][key]
        if st is None or te is None:
            dec[m] = "NO_MEDIDA"
        else:
            dec[m] = "CABE" if st <= GATE_STATES and te <= GATE_TERMS else "NO_CABE"
    out["proyecciones_k4"] = res
    out["decision_D_c1"] = dec


def run(out, args):
    ce, rect, b1 = b4.load_all(with_b1=True)
    ctx = mp.get_context("spawn")
    R = json.load(open(B4DIR + "/RESULTADO_B4.json", encoding="utf-8"))["instancias"]
    cells_rates = [(rc, rm) for rc in ce.GRID_P for rm in (rc, F(0))]
    I = out["instancias"]
    control_fen(ctx, ce, rect, b1, R, I)
    if args.prueba:
        run_strip(ctx, ce, rect, b1, 6, 3, R, args, I, cells_rates)
        return
    run_strip(ctx, ce, rect, b1, 6, 4, R, args, I, cells_rates)
    run_strip(ctx, ce, rect, b1, 12, 4, R, args, I, cells_rates)
    try:
        run_strip(ctx, ce, rect, b1, 24, 3, R, args, I, cells_rates)
    except b4.Stop as ex:
        LIMITS.append("TIRA_24x6 NO_MEDIDA: %s" % ex)
        I["TIRA_24x6"] = {"no_medida": str(ex)}
        print("  TIRA_24x6 NO_MEDIDA: %s" % ex, flush=True)
    project(I, out)
    fs = {n: I[n].get("f_s_circuito") for n in ("TIRA_6x6", "TIRA_12x6", "TIRA_24x6") if n in I}
    out["decision_D_c0"] = "PASS"
    out["decision_D_c2"] = {n: (v if v is not None else "NO_MEDIDA") for n, v in fs.items()}
    out["lectura_global"] = "D-c0 PASS; D-c1 %s; D-c2 f_s %s" % (
        ", ".join("m=%s %s" % (m, v) for m, v in out["decision_D_c1"].items()),
        ", ".join("%s=%s" % (n, v) for n, v in out["decision_D_c2"].items()))
    print("LECTURA:", out["lectura_global"], flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--recursos", required=True)
    ap.add_argument("--prueba", action="store_true")
    ap.add_argument("--no-exh", dest="no_exh", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args()
    out = {"unit": "T2-0", "plan": PLAN, "plan_sha256": SHA[PLAN], "referencias_sha256": SHA,
           "programa_sha256": sha_file(os.path.abspath(__file__)), "prueba": bool(args.prueba), "D": str(D),
           "instancias": {}}
    status = "FAIL"
    try:
        run(out, args)
        status = "PASS" if not LIMITS else "PASS_CON_NO_MEDIDAS"
    except b4.CheckFail as ex:
        out["detenido_en"] = "control fallido: " + str(ex)
        out["decision_D_c0"] = "FALLA"
    except b4.Stop as ex:
        out["detenido_en"] = str(ex)
        status = "LIMITE_DE_RECURSOS"
    out["status"] = status
    out["limites_alcanzados"] = LIMITS
    out["checks"] = b4.CHECKS
    out["llamadas_totales"] = b4.CALLS[0]
    with open(args.output, "w", encoding="ascii", newline="\n") as fh:
        fh.write(json.dumps(out, ensure_ascii=True, indent=1, sort_keys=True, default=str) + "\n")
    with open(args.recursos, "w", encoding="ascii", newline="\n") as fh:
        fh.write(json.dumps({"nota": "tiempos de reloj en segundos y memoria; no deterministas; no son resultado",
                             "tiempos": dict(b4.TIMES, total=round(time.perf_counter() - T0, 3)),
                             "pico_memoria_principal_bytes": b4.mem_now(peak=True),
                             "pico_memoria_procesos_bytes": b4.WPEAK, "llamadas_totales": b4.CALLS[0]},
                            ensure_ascii=True, indent=1, sort_keys=True) + "\n")
    print("%s controles=%d llamadas=%d t=%.1fs" % (status, len(b4.CHECKS), b4.CALLS[0], time.perf_counter() - T0),
          flush=True)


if __name__ == "__main__":
    main()
