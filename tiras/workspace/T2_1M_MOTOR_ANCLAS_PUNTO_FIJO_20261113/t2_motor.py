#!/usr/bin/env python3
"""T2-1M (PLAN_T2_1M.txt, sellado en SELLO_PLAN.json): motor por anclas con punto fijo dirigido.

Unidad de ingenieria. Construye el catalogo de las tiras en disco, partido en 256 grupos de anclas, con tres
recuentos (n_all, n_plus, n_minus) por soporte Y y firma de clases. Recorre la DP de B1 hacia adelante, ancla por
ancla (lema A), con la recurrencia normalizada (lema N), los canales E/O (lema EO) y redondeo dirigido a P = 384
bits fraccionarios (lema R). Valida el motor (C-anc, C-B4, C-T20, C-exacto, C-dom, C-anch, I4 y leyes) y mide si
T2-1 cabe (D-m2). No calcula intervalos de P_sil.

Reutiliza sin cambios, tras comprobar su SHA-256: t2_calibracion.py de T2-0 (calendario, leyes, clases, classify,
connected, criterion y controles del circuito) y, a traves de el, b4_tiras.py (tiras y utilidades),
control_exacto.py (decoder SSF) y b1_exacto_k4.py (anchor_masks y contained).

Uso (desde D:\\paper ia cuantica\\qec_revision_matematica\\05_COMPROBACIONES):
  python -I -B T2_1M_MOTOR_ANCLAS_PUNTO_FIJO_20261113/t2_motor.py --output <nuevo.json>
      --recursos <nuevo.json> --trabajo <directorio nuevo> [--prueba]
  (--prueba: prueba de desarrollo con m = 6 (y m = 12 fenomenologico); no es la corrida de la unidad)
"""
import argparse
import bisect
import ctypes
import gc
import hashlib
import json
import math
import multiprocessing as mp
import os
import pickle
import sys
import threading
import time
import zlib
from array import array
from decimal import ROUND_CEILING, ROUND_FLOOR, Decimal, localcontext
from fractions import Fraction as F
from itertools import product
from operator import mul

sys.set_int_max_str_digits(0)
T0 = time.perf_counter()
UDIR = "T2_1M_MOTOR_ANCLAS_PUNTO_FIJO_20261113/"
T20 = "T2_0_CALIBRACION_CIRCUITO_TIRAS_20261113/"
B4D = "B4_TIRAS_ESCALADO_N_20261113/"
PLAN = UDIR + "PLAN_T2_1M.txt"
SHA = {
    PLAN: "81f475fe847d9968239d79fa6196f581f682da7b66013ce3cff83d6485ef27a7",
    T20 + "t2_calibracion.py": "195bf3317f7af114f4cdc1daad714e828534ddd57abab29c81daa7126d17cc61",
    T20 + "RESULTADO_T2_0.json": "fb56c4324034549476813bd1a0d2bacd6f46c27d552fbf71ba2d0813d62df38e",
    B4D + "b4_tiras.py": "b1a2cb8f93affcc3713eadd4ff860a2d74756b91752011cd6b16da44def2dd15",
    B4D + "RESULTADO_B4.json": "640db376441656cfb72e16aa00e0e99e2fc2eb6047905dfe031f5d994634b5f7",
}


def sha_file(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


for _p, _s in SHA.items():
    if sha_file(_p) != _s:
        raise SystemExit("SHA-256 distinto del sellado: " + _p)
sys.path.insert(0, os.path.abspath(T20))
import t2_calibracion as t2  # noqa: E402  (tras comprobar su SHA-256)

b4 = t2.b4
bits = b4.bits
check = b4.check

PREC = 384
NGROUPS = 256
CHUNKS_F = 64 * 14
PROCS = 14
KMAX = 4
GIB = 1 << 30
SEG_MAX = 8 * 3600.0
CALLS_MAX = 7_500_000_000
MEM_TOTAL = 24 * GIB
MEM_MAIN = 12 * GIB
MEM_CAT = 3 * GIB // 2
MEM_DP = 4 * GIB
DISK_MAX = 60 * 10 ** 9
REL_W = F(1, 10 ** 30)
SUM_W = F(1, 10 ** 12)
DOM_LO, DOM_HI = 1 - F(1, 10 ** 12), 1 + F(1, 10 ** 2)
TOL_MARGEN, TOL_PESO = 1e-6, 1e-9
CAL_GROUPS = 16
CELL_CAL = (F(1, 10 ** 4), F(1, 10 ** 4))
CELL_DOM = (F(1, 10 ** 6), F(1, 10 ** 6))
CELLS_EXACT = ((F(1, 10 ** 3), F(1, 10 ** 3)), (F(1, 10 ** 6), F(0)))
T21_CELLS = 8
HOURS_OK = 6.0
FEN_SPLIT = {6: 1, 12: 2, 24: 4, 48: 6}     # trabajos de C-B4 con k = 4 (celdas repartidas)
LIMITS = []
TIMES = {}
MEMPEAK = {}
OVER = [None]
DISK = [0, 0]          # bytes actuales y maximos en trabajo/
b4.T0, b4.SEG_MAX, b4.MEM_MAIN, b4.LLAMADAS_MAX = T0, SEG_MAX, MEM_MAIN, CALLS_MAX


def guard(where):
    if time.perf_counter() - T0 > SEG_MAX:
        raise b4.Stop("techo de tiempo alcanzado en " + where)
    m = b4.mem_now()
    if m is not None and m > MEM_MAIN:
        raise b4.Stop("techo de memoria del proceso principal alcanzado en " + where)
    if OVER[0]:
        raise b4.Stop("techo de memoria total alcanzado (%s) en %s" % (OVER[0], where))
    if DISK[0] > DISK_MAX:
        raise b4.Stop("techo de disco alcanzado en " + where)
    if b4.CALLS[0] > CALLS_MAX:
        raise b4.Stop("techo de llamadas alcanzado en " + where)


def mem_pid(pid):
    """Working set de otro proceso (Windows), en bytes; 0 si no se puede medir."""
    try:
        from ctypes import wintypes

        class PMC(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                        ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]
        k32 = ctypes.windll.kernel32
        k32.OpenProcess.restype = wintypes.HANDLE
        k32.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
        k32.K32GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(PMC), wintypes.DWORD]
        k32.CloseHandle.argtypes = [wintypes.HANDLE]
        h = k32.OpenProcess(0x1000 | 0x0010, False, pid)
        if not h:
            return 0
        try:
            pmc = PMC()
            pmc.cb = ctypes.sizeof(PMC)
            if k32.K32GetProcessMemoryInfo(h, ctypes.byref(pmc), pmc.cb):
                return int(pmc.WorkingSetSize)
        finally:
            k32.CloseHandle(h)
    except Exception:
        pass
    return 0


class Monitor:
    """Memoria total (principal + procesos del pool) cada 5 s; marca OVER si supera el techo."""

    def __init__(self, pool, stage):
        self.pool, self.stage = pool, stage
        self.ev = threading.Event()
        self.th = threading.Thread(target=self.loop, daemon=True)
        self.th.start()

    def loop(self):
        while not self.ev.wait(5.0):
            try:
                tot = (b4.mem_now() or 0) + sum(mem_pid(p.pid) for p in list(self.pool._pool) if p.pid)
            except Exception:
                continue
            if tot > MEMPEAK.get(self.stage, 0):
                MEMPEAK[self.stage] = tot
            if tot > MEM_TOTAL and OVER[0] is None:
                OVER[0] = "%d bytes en %s" % (tot, self.stage)

    def stop(self):
        self.ev.set()
        self.th.join()


# ------------------------------------------------------------------ orden de producto, anclas y firmas
def layout(inst):
    coord, NS = inst["coord"], inst["NS"]
    order = sorted(range(NS), key=lambda v: (coord[v][0], coord[v][1], v))
    pos = [0] * NS
    for i, v in enumerate(order):
        pos[v] = i
    gsize = -(-NS // NGROUPS)
    return order, pos, gsize, -(-NS // gsize)


def relabel(msk, pos):
    r = 0
    for x in bits(msk):
        r |= 1 << pos[x]
    return r


def sig_syms(sig):
    out = []
    while sig:
        out.append((sig & 63) - 1)
        sig >>= 6
    return out


def anchor_of(y):
    return (y & -y).bit_length() - 1


# ------------------------------------------------------------------ procesos del catalogo
INC_P = 1 << 32
INC_M = 1 << 64
M32 = (1 << 32) - 1


def w_init_m(m, restricted):
    t2.w_init(m, restricted)
    W = t2._W
    inst = W["inst"]
    order, pos, gsize, _ = layout(inst)
    adjR = [0] * inst["NS"]
    for v in range(inst["NS"]):
        adjR[pos[v]] = relabel(inst["adjK"][v], pos)
    W.update(order=order, pos=pos, gsize=gsize, adjR=adjR)


def w_cat(task):
    """Mismas pruebas y misma regla de validez que w_grow de T2-0; agrega por (Y, firma) y reparte por ancla."""
    f, chunk, want_rec = task
    W = t2._W
    inst, CL, LT, pos, gsize = W["inst"], W["CL"], W["LT"], W["pos"], W["gsize"]
    adj = inst["adjK"]
    classify, connected, CONN = t2.classify, t2.connected, W["CONN"]
    poly, rec = {}, {}
    calls = bad = nval = 0
    for sm in chunk:
        sites = list(bits(sm))
        yu = 0
        for idxs in product(*[range(len(CL[b])) for b in sites]):
            Y, passes, silent, chis, b_ = classify(sm, sites, idxs)
            calls += 1
            bad += b_
            if not connected(Y, adj, CONN):
                continue
            nval += 1
            yu |= Y
            sig = 0
            for s in sorted(LT[b] * 8 + i for b, i in zip(sites, idxs)):
                sig = (sig << 6) | (s + 1)
            inc = 1 if not passes else (1 + INC_P if chis[0] > 0 else 1 + INC_M)
            slot = poly.get(Y)
            if slot is None:
                poly[Y] = {sig: inc}
            else:
                slot[sig] = slot.get(sig, 0) + inc
        if yu and want_rec:
            rec[sm] = yu
    W["RES"].clear()
    byanc = {}
    nent = 0
    for Y, slot in poly.items():
        YR = relabel(Y, pos)
        sg = sorted(slot)
        vals = [slot[s] for s in sg]
        nent += len(sg)
        byanc.setdefault(anchor_of(YR), []).append(
            (YR, array("I", sg), array("I", [v & M32 for v in vals]), array("I", [(v >> 32) & M32 for v in vals]),
             array("I", [v >> 64 for v in vals])))
    npoly = len(poly)
    del poly
    blobs = {a: zlib.compress(pickle.dumps(lst, protocol=5), 1) for a, lst in byanc.items()}
    del byanc
    m = b4.mem_now()
    if m is not None and m > MEM_CAT:
        raise RuntimeError("techo de memoria de un proceso del catalogo (%d bytes)" % m)
    return {"blobs": blobs, "rec": rec, "calls": calls, "bad": bad, "nval": nval, "npoly": npoly, "nent": nent,
            "peak": b4.mem_now(peak=True), "pid": os.getpid()}


def w_merge(task):
    """Fusiona, ancla por ancla, los trozos de un grupo (soportes repetidos entre trozos) y escribe una seccion en
    columnas por ancla: [ancla (4 bytes)][longitud (8 bytes)][zlib(pickle(seccion))]."""
    g, raw, out, want_neg = task
    W = t2._W
    NS, gsize, adjR, order = W["inst"]["NS"], W["gsize"], W["adjR"], W["order"]
    lo, hi = g * gsize, min(NS, (g + 1) * gsize)
    idx = {}
    dup = wrong = npoly = nbytes = 0
    npk, nek = [0] * (KMAX + 1), [0] * (KMAX + 1)
    sigs, per_anchor = set(), {}
    neg = None
    hsh = hashlib.sha256()
    with open(raw, "rb") as fh, open(out, "wb") as fo:
        while True:
            h = fh.read(12)
            if len(h) < 12:
                break
            a, n = int.from_bytes(h[:4], "little"), int.from_bytes(h[4:], "little")
            idx.setdefault(a, []).append((fh.tell(), n))
            fh.seek(n, 1)
        outside = sum(len(v) for a, v in idx.items() if not lo <= a < hi)
        for a in sorted(idx):
            D = {}
            for p0, n in idx[a]:
                fh.seek(p0)
                lst = pickle.loads(zlib.decompress(fh.read(n)))
                for rec in lst:
                    YR = rec[0]
                    if anchor_of(YR) != a:
                        wrong += 1
                    cur = D.get(YR)
                    if cur is None:
                        D[YR] = rec[1:]
                        continue
                    dup += 1
                    acc = {}
                    for arrs in (cur, rec[1:]):
                        for s, x, y, z in zip(*arrs):
                            e = acc.get(s)
                            acc[s] = (x, y, z) if e is None else (e[0] + x, e[1] + y, e[2] + z)
                    ks = sorted(acc)
                    D[YR] = (array("I", ks), array("I", [acc[s][0] for s in ks]),
                             array("I", [acc[s][1] for s in ks]), array("I", [acc[s][2] for s in ks]))
                del lst
            YRs = sorted(D)
            NBs, off = [], array("Q", [0])
            sig, nA, nP, nM = array("I"), array("I"), array("I"), array("I")
            for YR in YRs:
                sg, x, y, z = D.pop(YR)
                nb = YR
                for i in bits(YR):
                    nb |= adjR[i]
                NBs.append(nb)
                sig.extend(sg)
                nA.extend(x)
                nP.extend(y)
                nM.extend(z)
                off.append(len(sig))
                for k in range(1, KMAX + 1):
                    cap = 1 << (6 * k)
                    if sg[0] < cap:
                        npk[k] += 1
                        nek[k] += bisect.bisect_left(sg, cap)
                if want_neg and sg[-1] >= 1 << 18:
                    Y = 0
                    for i in bits(YR):
                        Y |= 1 << order[i]
                    if neg is None or Y < neg[0]:
                        neg = (Y, YR)
            del D
            sigs.update(sig)
            npoly += len(YRs)
            per_anchor[a] = len(YRs)
            blob = zlib.compress(pickle.dumps({"a": a, "YR": YRs, "NB": NBs, "off": off, "sig": sig, "nA": nA,
                                               "nP": nP, "nM": nM}, protocol=5), 1)
            del YRs, NBs, off, sig, nA, nP, nM
            head = a.to_bytes(4, "little") + len(blob).to_bytes(8, "little")
            fo.write(head)
            fo.write(blob)
            hsh.update(head)
            hsh.update(blob)
            nbytes += len(head) + len(blob)
            del blob
            m = b4.mem_now()
            if m is not None and m > MEM_CAT:
                raise RuntimeError("techo de memoria de un proceso de la fusion (%d bytes)" % m)
    return {"g": g, "polimeros": npoly, "duplicados_fusionados": dup, "fuera_de_grupo": outside + wrong,
            "pol_k": npk, "ent_k": nek, "firmas": sorted(sigs), "bytes": nbytes, "sha256": hsh.hexdigest(),
            "polimeros_por_ancla": per_anchor, "neg": neg, "peak": b4.mem_now(peak=True), "pid": os.getpid()}


# ------------------------------------------------------------------ catalogo en disco (regla de B2, como T2-0)
def build_catalogue(pool, inst, name, nclass, kmax, wdir, want_neg):
    NS, adj = inst["NS"], inst["adjK"]
    _, _, gsize, ng = layout(inst)
    os.makedirs(wdir, exist_ok=False)

    def raw(g):
        return os.path.join(wdir, "crudo_g%03d.bin" % g)

    fhs = {}
    rec_prev, info, bad = None, {}, 0
    t_cat = time.perf_counter()
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
        target = max(1, need // CHUNKS_F)
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
        rec_new = {}
        nval = calls = nent = 0
        try:
            for res in pool.imap_unordered(w_cat, [(f, ch, want_rec) for ch in chunks]):
                calls += res["calls"]
                bad += res["bad"]
                nval += res["nval"]
                nent += res["nent"]
                b4.note_peak("catalogo_" + name, res)
                for a, blob in res["blobs"].items():
                    g = a // gsize
                    fh = fhs.get(g)
                    if fh is None:
                        fh = open(raw(g), "ab")
                        fhs[g] = fh
                    fh.write(a.to_bytes(4, "little") + len(blob).to_bytes(8, "little"))
                    fh.write(blob)
                    DISK[0] += 12 + len(blob)
                DISK[1] = max(DISK[1], DISK[0])
                for A, yu in res["rec"].items():
                    rec_new[A] = rec_new.get(A, 0) | yu
                res.clear()
                guard("catalogo de %s f=%d" % (name, f))
        except (b4.Stop, b4.CheckFail):
            raise
        except Exception as ex:
            raise b4.Stop("catalogo de %s f=%d: %s" % (name, f, ex))
        b4.CALLS[0] += calls
        check("%s f=%d: cada asignacion de cada candidato se decodifica una vez (%d)" % (name, f, need),
              calls == need, verbose=False)
        info[str(f)] = {"candidatos": len(cand), "asignaciones": need, "validas": nval,
                        "conjuntos_validos": len(rec_new) if want_rec else None, "entradas_por_trozo": nent,
                        "trozos": len(chunks), "segundos": round(time.perf_counter() - t, 3)}
        print("  %s f=%d: %d candidatos, %d asignaciones, %d validas (%.1fs)" % (
            name, f, len(cand), need, nval, time.perf_counter() - t), flush=True)
        rec_prev = rec_new
        del cand, chunks, wts
    for fh in fhs.values():
        fh.close()
    check("I4 %s: silenciosa <=> caracter logico -1 (discrepancias = 0)" % name, bad == 0, bad)
    seg_cat = time.perf_counter() - t_cat
    t = time.perf_counter()
    stats = []
    tasks = [(g, raw(g), os.path.join(wdir, "grupo_%03d.bin" % g), want_neg) for g in sorted(fhs)]
    try:
        for res in pool.imap_unordered(w_merge, tasks):
            stats.append(res)
            DISK[0] += res["bytes"]
            DISK[1] = max(DISK[1], DISK[0])
            b4.note_peak("fusion_" + name, res)
            guard("fusion de " + name)
    except (b4.Stop, b4.CheckFail):
        raise
    except Exception as ex:
        raise b4.Stop("fusion de %s: %s" % (name, ex))
    for g in fhs:
        DISK[0] -= os.path.getsize(raw(g))
        os.remove(raw(g))
    stats.sort(key=lambda r: r["g"])
    outside = sum(r["fuera_de_grupo"] for r in stats)
    check("C-anc %s: cada soporte fusionado esta en el grupo de su ancla (fuera de grupo = 0; %d grupos)"
          % (name, len(stats)), outside == 0, outside)
    sigs = sorted(set().union(*[r["firmas"] for r in stats]))
    negs = [r["neg"] for r in stats if r["neg"] is not None]
    cat = {"nombre": name, "dir": wdir, "NS": NS, "gsize": gsize, "ngroups": ng,
           "manifiesto": {str(r["g"]): {"archivo": "grupo_%03d.bin" % r["g"], "sha256": r["sha256"],
                                         "bytes": r["bytes"], "polimeros": r["polimeros"]} for r in stats},
           "firmas": sigs, "por_f": info,
           "polimeros_por_k": {str(k): sum(r["pol_k"][k] for r in stats) for k in range(1, kmax + 1)},
           "entradas_por_k": {str(k): sum(r["ent_k"][k] for r in stats) for k in range(1, kmax + 1)},
           "bytes": sum(r["bytes"] for r in stats), "duplicados_fusionados": sum(r["duplicados_fusionados"]
                                                                                  for r in stats),
           "polimeros_por_ancla_max": max((max(r["polimeros_por_ancla"].values(), default=0) for r in stats),
                                          default=0),
           "neg": min(negs) if negs else None, "segundos_catalogo": round(seg_cat, 3),
           "segundos_fusion": round(time.perf_counter() - t, 3)}
    print("  %s: %d polimeros (k=4), %d entradas, %.2f GB, %d firmas; fusion %.1fs" % (
        name, cat["polimeros_por_k"][str(kmax)], cat["entradas_por_k"][str(kmax)], cat["bytes"] / 1e9, len(sigs),
        time.perf_counter() - t), flush=True)
    return cat


def read_group(path, sha):
    """Lee un grupo, comprueba su SHA-256 y devuelve (bytes, {ancla: (inicio, longitud)})."""
    with open(path, "rb") as fh:
        blob = fh.read()
    if hashlib.sha256(blob).hexdigest() != sha:
        raise RuntimeError("SHA-256 distinto en " + path)
    secs, p = {}, 0
    while p < len(blob):
        a, n = int.from_bytes(blob[p:p + 4], "little"), int.from_bytes(blob[p + 4:p + 12], "little")
        secs[a] = (p + 12, n)
        p += 12 + n
    return blob, secs


def load_section(cat, a):
    g = a // cat["gsize"]
    ent = cat["manifiesto"][str(g)]
    blob, secs = read_group(os.path.join(cat["dir"], ent["archivo"]), ent["sha256"])
    p, n = secs[a]
    return g, pickle.loads(zlib.decompress(memoryview(blob)[p:p + n]))


# ------------------------------------------------------------------ DP hacia adelante por anclas
def dp_job(spec):
    try:
        return _dp(spec)
    except Exception as ex:
        return {"etiqueta": spec["etiqueta"], "error": "%s: %s" % (type(ex).__name__, ex), "pid": os.getpid(),
                "peak": b4.mem_now(peak=True)}


def _dp(spec):
    """Lemas A, N, EO y R. Modos: estructura (estados y terminos), celdas (U, E, O con cotas inferior y superior
    a 2^-P; opcionalmente exactas) y suma (U con odds 1/D)."""
    t_all = time.perf_counter()
    b1 = b4.load_all(with_b1=True)[2]
    anchor_masks, contained = b1.anchor_masks, b1.contained
    P = PREC
    NS, gsize, ng, k, mode = spec["NS"], spec["gsize"], spec["ngroups"], spec["k"], spec["modo"]
    cap = 1 << (6 * k)
    cells, exact = spec["celdas"], spec["exacto"]
    nc = len(cells)
    if mode == "suma" and (nc != 1 or exact):
        raise ValueError("modo suma: una celda, sin exacto")
    tabs = []
    for c in cells:
        WL, DL, WX = {}, {}, {}
        ina = True
        for s in spec["firmas"]:
            w = F(1)
            for sym in sig_syms(s):
                w *= c["odds"][sym]
            q, r = divmod(w.numerator << P, w.denominator)
            WL[s] = q
            DL[s] = 1 if r else 0
            if not r:
                ina = False
            if exact:
                WX[s] = w
        tabs.append((WL, DL, WX, ina))
    nslot = 0 if mode == "estructura" else (2 if mode == "suma" else 6 * nc)
    nx = 3 * nc if exact else 0
    ONE2 = 1 << (2 * P)
    ALL = (1 << NS) - 1
    if mode == "estructura":
        pend = {0: {ALL}}
    else:
        init = [ONE2, ONE2] if mode == "suma" else [ONE2, ONE2, ONE2, ONE2, 0, 0] * nc
        if exact:
            init = init + [F(1), F(1), F(0)] * nc
        pend = {0: {ALL: init}}
    zero = [0] * nslot + [F(0)] * nx
    RS = range(nslot)
    R6 = range(0, nslot, 6)
    RXC = range(nslot, nslot + nx, 3)
    states = terms = badA = polys = nent = removed = 0
    npend = maxpend = 1
    maxpoly = 0
    t_load = t_w = 0.0
    bytes_loaded = 0
    bps = 0.0
    limit = spec["limite_grupos"]
    gmax = ng if limit is None else min(limit, ng)
    man, rmv, deadline = spec["manifiesto"], spec["quitar"], spec["plazo"]
    t_loop = time.perf_counter()
    nanch = 0
    BYe = [()] * NS
    for g in range(gmax):
        lo, hi = g * gsize, min(NS, (g + 1) * gsize)
        if not any(x in pend for x in range(lo, hi)):
            continue
        secs, mv, blob = {}, None, None
        ent = man.get(str(g))
        if ent is not None:
            t1 = time.perf_counter()
            blob, secs = read_group(os.path.join(spec["dir"], ent["archivo"]), ent["sha256"])
            bytes_loaded += len(blob)
            mv = memoryview(blob)
            t_load += time.perf_counter() - t1
        rm = rmv if (rmv is not None and rmv[0] == g) else None
        for x in range(lo, hi):
            Dx = pend.pop(x, None)
            if Dx is None:
                continue
            YRk, NBk, Wk = [], [], []
            sec = secs.get(x)
            if sec is not None:
                t1 = time.perf_counter()
                data = pickle.loads(zlib.decompress(mv[sec[0]:sec[0] + sec[1]]))
                t2_ = time.perf_counter()
                t_load += t2_ - t1
                YR, NB, off, sig = data["YR"], data["NB"], data["off"], data["sig"]
                nA, nP, nM = data["nA"], data["nP"], data["nM"]
                if data["a"] != x:
                    raise RuntimeError("seccion con ancla distinta")
                for i in range(len(YR)):
                    a, b = off[i], off[i + 1]
                    e = bisect.bisect_left(sig, cap, a, b)
                    if e == a:
                        continue
                    y = YR[i]
                    if anchor_of(y) != x:
                        raise RuntimeError("polimero fuera de su ancla")
                    A = nA[a:e]
                    if rm is not None and y == rm[1]:
                        j = bisect.bisect_left(sig, rm[2], a, e)
                        if j == e or sig[j] != rm[2] or A[j - a] == 0:
                            raise RuntimeError("decoracion a quitar no encontrada")
                        A[j - a] -= 1
                        removed += 1
                        if not any(A):
                            continue
                    YRk.append(y)
                    NBk.append(NB[i])
                    nent += e - a
                    if nslot:
                        sl = sig[a:e]
                        w = []
                        wx = []
                        if mode != "suma":
                            Pp, Mm = nP[a:e], nM[a:e]
                            sameP = Pp == A
                            hasM = any(Mm)
                        for WL, DL, WX, ina in tabs:
                            wl = list(map(WL.__getitem__, sl))
                            Wl = sum(map(mul, A, wl))
                            Wh = Wl + (sum(A) if ina else sum(map(mul, A, map(DL.__getitem__, sl))))
                            if mode == "suma":
                                w += (Wl, Wh)
                                continue
                            if sameP:
                                Pl, Ph = Wl, Wh
                            else:
                                Pl = sum(map(mul, Pp, wl))
                                Ph = Pl + (sum(Pp) if ina else sum(map(mul, Pp, map(DL.__getitem__, sl))))
                            if hasM:
                                Ml = sum(map(mul, Mm, wl))
                                Mh = Ml + (sum(Mm) if ina else sum(map(mul, Mm, map(DL.__getitem__, sl))))
                            else:
                                Ml = Mh = 0
                            w += (Wl, Wh, Pl, Ph, Ml, Mh)
                            if exact:
                                ws = [WX[s] for s in sl]
                                wx += (sum(map(mul, A, ws), F(0)), sum(map(mul, Pp, ws), F(0)),
                                       sum(map(mul, Mm, ws), F(0)))
                        Wk.append(w + wx)
                del data, YR, NB, off, sig, nA, nP, nM
                t_w += time.perf_counter() - t2_
            npl = len(YRk)
            polys += npl
            if npl > maxpoly:
                maxpoly = npl
            BYe[x] = range(npl)
            WIN, FULL, CM = anchor_masks({"NS": NS, "YR": YRk, "BY": BYe})
            BYe[x] = ()
            xb = 1 << x
            if mode == "estructura":
                for J in Dx:
                    states += 1
                    npend -= 1
                    ch = J ^ xb
                    a = anchor_of(ch) if ch else NS
                    if a <= x:
                        badA += 1
                    d = pend.get(a)
                    if d is None:
                        pend[a] = {ch}
                        npend += 1
                    elif ch not in d:
                        d.add(ch)
                        npend += 1
                    c = contained(J, x, WIN, FULL, CM)
                    if c:
                        for wi, wd in enumerate(memoryview(c.to_bytes(((c.bit_length() + 63) >> 6) << 3, "little")).cast("Q")):
                            while wd:
                                bb = wd & -wd
                                wd ^= bb
                                ch = J & ~NBk[(wi << 6) | (bb.bit_length() - 1)]
                                terms += 1
                                a = (ch & -ch).bit_length() - 1 if ch else NS
                                if a <= x:
                                    badA += 1
                                d = pend.get(a)
                                if d is None:
                                    pend[a] = {ch}
                                    npend += 1
                                elif ch not in d:
                                    d.add(ch)
                                    npend += 1
                    if npend > maxpend:
                        maxpend = npend
            elif mode == "suma":
                for J, acc in Dx.items():
                    states += 1
                    npend -= 1
                    vl = acc[0] >> P
                    vh = -((-acc[1]) >> P)
                    ch = J ^ xb
                    a = anchor_of(ch) if ch else NS
                    if a <= x:
                        badA += 1
                    d = pend.get(a)
                    if d is None:
                        d = {}
                        pend[a] = d
                    t = d.get(ch)
                    if t is None:
                        t = [0, 0]
                        d[ch] = t
                        npend += 1
                    t[0] += vl << P
                    t[1] += vh << P
                    c = contained(J, x, WIN, FULL, CM)
                    if c:
                        for wi, wd in enumerate(memoryview(c.to_bytes(((c.bit_length() + 63) >> 6) << 3, "little")).cast("Q")):
                            while wd:
                                bb = wd & -wd
                                wd ^= bb
                                pi = (wi << 6) | (bb.bit_length() - 1)
                                ch = J & ~NBk[pi]
                                terms += 1
                                a = (ch & -ch).bit_length() - 1 if ch else NS
                                if a <= x:
                                    badA += 1
                                d = pend.get(a)
                                if d is None:
                                    d = {}
                                    pend[a] = d
                                t = d.get(ch)
                                if t is None:
                                    t = [0, 0]
                                    d[ch] = t
                                    npend += 1
                                w = Wk[pi]
                                t[0] += vl * w[0]
                                t[1] += vh * w[1]
                    if npend > maxpend:
                        maxpend = npend
            else:
                for J, acc in Dx.items():
                    states += 1
                    npend -= 1
                    v = [acc[i] >> P if not i & 1 else -((-acc[i]) >> P) for i in RS]
                    ch = J ^ xb
                    a = anchor_of(ch) if ch else NS
                    if a <= x:
                        badA += 1
                    d = pend.get(a)
                    if d is None:
                        d = {}
                        pend[a] = d
                    t = d.get(ch)
                    if t is None:
                        t = zero[:]
                        d[ch] = t
                        npend += 1
                    for i in RS:
                        t[i] += v[i] << P
                    for i in range(nslot, nslot + nx):
                        t[i] += acc[i]
                    c = contained(J, x, WIN, FULL, CM)
                    if c:
                        for wi, wd in enumerate(memoryview(c.to_bytes(((c.bit_length() + 63) >> 6) << 3, "little")).cast("Q")):
                            while wd:
                                bb = wd & -wd
                                wd ^= bb
                                pi = (wi << 6) | (bb.bit_length() - 1)
                                ch = J & ~NBk[pi]
                                terms += 1
                                a = (ch & -ch).bit_length() - 1 if ch else NS
                                if a <= x:
                                    badA += 1
                                d = pend.get(a)
                                if d is None:
                                    d = {}
                                    pend[a] = d
                                t = d.get(ch)
                                if t is None:
                                    t = zero[:]
                                    d[ch] = t
                                    npend += 1
                                w = Wk[pi]
                                for b in R6:
                                    El, Eh, Ol, Oh = v[b + 2], v[b + 3], v[b + 4], v[b + 5]
                                    Pl, Ph, Mh = w[b + 2], w[b + 3], w[b + 5]
                                    t[b] += v[b] * w[b]
                                    t[b + 1] += v[b + 1] * w[b + 1]
                                    if Mh:
                                        Ml = w[b + 4]
                                        t[b + 2] += El * Pl + Ol * Ml
                                        t[b + 3] += Eh * Ph + Oh * Mh
                                        t[b + 4] += El * Ml + Ol * Pl
                                        t[b + 5] += Eh * Mh + Oh * Ph
                                    else:
                                        t[b + 2] += El * Pl
                                        t[b + 3] += Eh * Ph
                                        if Oh:
                                            t[b + 4] += Ol * Pl
                                            t[b + 5] += Oh * Ph
                                for bx in RXC:
                                    u, e_, o_ = acc[bx], acc[bx + 1], acc[bx + 2]
                                    Wx, Px, Mx = w[bx], w[bx + 1], w[bx + 2]
                                    t[bx] += u * Wx
                                    t[bx + 1] += e_ * Px + o_ * Mx
                                    t[bx + 2] += e_ * Mx + o_ * Px
                    if npend > maxpend:
                        maxpend = npend
            del Dx, YRk, NBk, Wk, WIN, FULL, CM
            nanch += 1
            if nanch % 16 == 0 and pend:
                big = max(pend.values(), key=len)
                smp = 0
                cnt = 0
                if mode == "estructura":
                    for J in big:
                        smp += sys.getsizeof(J) + 60
                        cnt += 1
                        if cnt >= 32:
                            break
                else:
                    for J, acc in big.items():
                        smp += sys.getsizeof(J) + sys.getsizeof(acc) + sum(sys.getsizeof(z) for z in acc) + 100
                        cnt += 1
                        if cnt >= 32:
                            break
                if cnt:
                    bps = max(bps, smp / cnt)
            if time.time() > deadline:
                raise RuntimeError("techo de tiempo de la unidad")
            m_ = b4.mem_now()
            if m_ is not None and m_ > MEM_DP:
                raise RuntimeError("techo de memoria de un proceso de la DP (%d bytes)" % m_)
        if mv is not None:
            mv.release()
        del secs, mv, blob
    t_loop = time.perf_counter() - t_loop
    out = {"etiqueta": spec["etiqueta"], "k": k, "modo": mode, "claves": [c["clave"] for c in cells],
           "completo": limit is None, "grupos": gmax, "terminos": terms, "polimeros": polys,
           "polimeros_por_ancla_max": maxpoly, "pendientes_max": maxpend, "anclas_malas": badA, "entradas": nent,
           "quitadas": removed, "bytes_cargados": bytes_loaded, "bytes_por_estado_pendiente": round(bps, 1),
           "segundos_carga": round(t_load, 3), "segundos_pesos": round(t_w, 3),
           "segundos_dp": round(t_loop - t_load - t_w, 3)}
    if limit is None:
        Df = pend.pop(NS, None)
        if Df is None or pend:
            raise RuntimeError("quedan estados sin procesar (%d anclas)" % len(pend))
        states += 1
        out["estados"] = states
        if nslot:
            acc = Df[0]
            out["valores"] = [[acc[b + i] >> P if not i & 1 else -((-acc[b + i]) >> P)
                               for i in range(6 if mode != "suma" else 2)] for b in range(0, nslot, 6)]
            if nx:
                out["exactos"] = [acc[bx:bx + 3] for bx in RXC]
    else:
        out["estados_procesados"] = states
    out["segundos"] = round(time.perf_counter() - t_all, 3)
    out["peak"] = b4.mem_now(peak=True)
    out["pid"] = os.getpid()
    return out


# ------------------------------------------------------------------ odds por celda y B
def odds_fen(cx, p, q):
    od = {}
    for t, lt in enumerate(cx["lts"]):
        od[t * 8] = p / (1 - p) if lt[0] == "d" else q / (1 - q)
    return od, (1 - p) ** cx["N"] * (1 - q) ** cx["M"]


def odds_circ(cx, rc, rm):
    od, B = {}, F(1)
    for t, lt in enumerate(cx["lts"]):
        law = t2.law_of(lt, cx["L"], rc, rm)
        z = law[tuple([0] * len(next(iter(law))))]
        for i, key in enumerate(cx["CKrep"][t]):
            od[t * 8 + i] = law[key] / z
        B *= z ** cx["cnt"][t]
    return od, B


def odds_sum(cx):
    return {t * 8 + i: F(1, t2.D) for t in range(len(cx["lts"])) for i in range(len(cx["CKrep"][t]))}


def dec(x, up):
    with localcontext() as c:
        c.prec = 25
        c.rounding = ROUND_CEILING if up else ROUND_FLOOR
        return str(Decimal(x.numerator) / Decimal(x.denominator))


def relw(lo, hi):
    if hi == 0 and lo == 0:
        return F(0)
    return abs(hi - lo) / max(abs(hi), abs(lo))


def cell_intervals(B, vals):
    S = 1 << PREC
    Ul, Uh, El, Eh, Ol, Oh = vals
    return {"eps_comp": (1 - B * F(Uh, S), 1 - B * F(Ul, S)), "L_comp": (B * F(Ol, S), B * F(Oh, S)),
            "P_0": (B * F(El + Ol, S), B * F(Eh + Oh, S))}


# ------------------------------------------------------------------ una instancia: leyes, clases y catalogo
def catalogue_instance(ctx, ce, rect, m, restricted, wroot, RB4, RT, grid8, I):
    name = "TIRA_%dx6" % m
    tag = ("fen_" if restricted else "circ_") + name
    print("== catalogo %s (k <= %d)" % (tag, KMAX), flush=True)
    inst = b4.make_inst(ce, rect, "tira", m)
    check("%s: matrices iguales a las de B4" % tag, b4.mats_sha(inst) == RB4[name]["matrices_sha256"])
    if restricted:
        layer, L = t2.schedule(inst)
        locs = t2.locations(inst, layer, L)
        CL, CK, LTK = t2.site_classes(inst, layer, L, locs, True)
        circ = None
    else:
        circ, (layer, L, locs, CL, CK, LTK) = t2.control_circuit(inst, tag, grid8)
    lts = sorted(set(LTK))
    LT = [lts.index(x) for x in LTK]
    NS = inst["NS"]
    check("%s: simbolos (tipo de ley, clase) caben en 6 bits (%d tipos, <= %d clases)" % (
        tag, len(lts), max(len(c) for c in CL)), len(lts) <= 7 and max(len(c) for c in CL) <= 8)
    CKrep = [CK[LT.index(t)] for t in range(len(lts))]
    check("%s: las clases de cada sitio son las de su tipo de ley (odds por simbolo bien definidos)" % tag,
          all(CK[b] == CKrep[LT[b]] for b in range(NS)))
    check("%s: k_log = 1" % tag, len(ce.logical_reps(inst)) == 1)
    nclass = [len(c) for c in CL]
    with ctx.Pool(PROCS, initializer=w_init_m, initargs=(m, restricted)) as pool:
        mon = Monitor(pool, "catalogo_" + tag)
        try:
            cat = build_catalogue(pool, inst, tag, nclass, KMAX, os.path.join(wroot, tag),
                                  want_neg=(not restricted and m == 6))
        finally:
            mon.stop()
    gc.collect()
    cat["cx"] = {"lts": lts, "cnt": [LT.count(t) for t in range(len(lts))], "CKrep": CKrep,
                 "L": L, "N": inst["N"], "M": inst["M"], "nclass": nclass}
    o = {"sitios": NS, "grupos": cat["ngroups"], "anclas_por_grupo": cat["gsize"], "circuito": circ,
         "catalogo": {"por_f": cat["por_f"], "polimeros_por_k": cat["polimeros_por_k"],
                      "entradas_por_k": cat["entradas_por_k"], "bytes": cat["bytes"], "firmas": len(cat["firmas"]),
                      "duplicados_fusionados": cat["duplicados_fusionados"],
                      "segundos_catalogo": cat["segundos_catalogo"], "segundos_fusion": cat["segundos_fusion"]}}
    if restricted:
        rc = RB4[name]["catalogo"]["por_f"]
        got = [cat["por_f"][str(f)]["validas"] for f in range(1, KMAX + 1)]
        exp = [rc[str(f)]["validas"] for f in range(1, KMAX + 1)]
        check("C-B4 %s: validas por f = B4 %s" % (tag, exp), got == exp, got)
        gk = RB4[name]["gas_por_k"]
        got = [cat["polimeros_por_k"][str(k)] for k in (3, 4)]
        exp = [gk[str(k)]["polimeros"] for k in (3, 4)]
        check("C-B4 %s: polimeros con k = 3 y 4 = B4 %s" % (tag, exp), got == exp, got)
    elif name in RT:
        kT = 4 if m in (6, 12) else 3
        ct = RT[name]["catalogo"]["por_f"]
        got = [[cat["por_f"][str(f)]["validas"], cat["por_f"][str(f)]["asignaciones"]] for f in range(1, kT + 1)]
        exp = [[ct[str(f)]["validas"], ct[str(f)]["asignaciones"]] for f in range(1, kT + 1)]
        check("C-T20 %s: validas y asignaciones por f (f <= %d) = T2-0" % (tag, kT), got == exp,
              {"obtenido": got, "T2-0": exp})
        got = [cat["polimeros_por_k"][str(k)] for k in range(1, kT + 1)]
        exp = [RT[name]["estructura"][str(k)]["polimeros"] for k in range(1, kT + 1)]
        check("C-T20 %s: polimeros por k (k <= %d) = T2-0 %s" % (tag, kT, exp), got == exp, got)
        if m == 6:
            cat["negativo"] = find_negative(cat)
            ng = RT[name]["regla_de_suma_negativo"]
            mine = [cat["negativo"][x] for x in ("soporte_hex", "clave", "recuento_antes")]
            exp = [ng["soporte_hex"], ng["clave"], ng["recuento_antes"]]
            check("C-T20 %s: decoracion quitada por la regla de T2-0 (soporte, clave, recuento) = %s" % (tag, exp),
                  mine == exp, mine)
            o["negativo"] = {x: cat["negativo"][x] for x in ("soporte_hex", "clave", "recuento_antes")}
    I[tag] = o
    return cat


def find_negative(cat):
    """remove_one de T2-0: soporte minimo (etiquetas originales) con una decoracion de 4 sitios, clave (na, nb)
    minima, una decoracion menos (n_all de la menor firma con esa clave; con odds iguales da el mismo U)."""
    Y, YR = cat["neg"]
    g, data = load_section(cat, anchor_of(YR))
    i = data["YR"].index(YR)
    keys = {}
    for j in range(data["off"][i], data["off"][i + 1]):
        s = data["sig"][j]
        if s < 1 << 18:
            continue
        syms = sig_syms(s)
        na = sum(1 for sy in syms if cat["cx"]["lts"][sy // 8][0] == "d")
        keys.setdefault((na, len(syms) - na), []).append((s, data["nA"][j]))
    key = min(keys)
    return {"soporte_hex": hex(Y), "clave": list(key), "recuento_antes": sum(n for _, n in keys[key]),
            "quitar": (g, YR, min(s for s, n in keys[key] if n > 0))}


# ------------------------------------------------------------------ trabajos de la DP
def make_specs(C, grid12, prueba):
    deadline = time.time() + (SEG_MAX - (time.perf_counter() - T0))
    S = []

    def add(cat, lab, k, mode, cells=(), exact=False, limit=None, quitar=None, prio=0.0):
        S.append({"etiqueta": lab, "dir": cat["dir"], "manifiesto": cat["manifiesto"], "NS": cat["NS"],
                  "gsize": cat["gsize"], "ngroups": cat["ngroups"], "k": k, "modo": mode, "celdas": list(cells),
                  "exacto": exact, "limite_grupos": limit, "quitar": quitar, "firmas": cat["firmas"],
                  "plazo": deadline, "prioridad": prio})

    for (kind, m), cat in C.items():
        if kind != "fen":
            continue
        cx = cat["cx"]
        cells = [{"clave": "%s|%s" % (p, q), "odds": odds_fen(cx, p, q)[0]} for p, q in grid12]
        add(cat, "fen_%d_k3" % m, 3, "celdas", cells, prio=3e4 * m * 13)
        ns = FEN_SPLIT[m]
        per = -(-len(cells) // ns)
        for j in range(ns):
            part = cells[j * per:(j + 1) * per]
            add(cat, "fen_%d_k4_%d" % (m, j), 4, "celdas", part, prio=4e6 * m * (1 + 2 * len(part)))
    for (kind, m), cat in C.items():
        if kind != "circ":
            continue
        cx = cat["cx"]
        kT = 4 if m in (6, 12) else 3
        lab = "circ_%d" % m
        if m <= 24:
            for k in range(1, kT):
                add(cat, "%s_est_k%d" % (lab, k), k, "estructura", prio=1e3 * m * k)
            add(cat, "%s_suma_k%d" % (lab, kT), kT, "suma", [{"clave": "suma_D", "odds": odds_sum(cx)}],
                prio=(2.3e7 if kT == 4 else 6e5) * m * 3)
        if m == 6:
            add(cat, "%s_negativo_k4" % lab, 4, "suma", [{"clave": "suma_D", "odds": odds_sum(cx)}],
                quitar=cat["negativo"]["quitar"], prio=2.3e7 * 6 * 3)
            for j, (rc, rm) in enumerate(CELLS_EXACT):
                add(cat, "%s_exacto_%d" % (lab, j), 3, "celdas", [{"clave": "%s|%s" % (rc, rm),
                                                                   "odds": odds_circ(cx, rc, rm)[0]}],
                    exact=True, prio=1e9)
        if m in (6, 12):
            add(cat, "%s_dom_k4" % lab, 4, "celdas", [{"clave": "%s|%s" % CELL_DOM,
                                                       "odds": odds_circ(cx, *CELL_DOM)[0]}],
                prio=2.3e7 * m * 4)
        if m in (12, 48) or prueba:
            add(cat, "%s_cal" % lab, 4, "celdas", [{"clave": "%s|%s" % CELL_CAL, "odds": odds_circ(cx, *CELL_CAL)[0]}],
                limit=CAL_GROUPS, prio=5e9)
        if m in (24, 48):
            add(cat, "%s_est_k4" % lab, 4, "estructura", prio=4e7 * m * 8)
        if m == 48:
            add(cat, "%s_est_k3" % lab, 3, "estructura", prio=1e6 * m)
    return S


def run_jobs(ctx, specs):
    R = {}
    specs = sorted(specs, key=lambda s: -s["prioridad"])
    print("== DP: %d trabajos en %d procesos" % (len(specs), PROCS), flush=True)
    with ctx.Pool(PROCS, maxtasksperchild=1) as pool:
        mon = Monitor(pool, "dp")
        try:
            for r in pool.imap_unordered(dp_job, specs, chunksize=1):
                lab = r["etiqueta"]
                R[lab] = r
                if "error" in r:
                    LIMITS.append("%s NO_MEDIDA: %s" % (lab, r["error"]))
                    print("  %s NO_MEDIDA: %s" % (lab, r["error"]), flush=True)
                else:
                    b4.note_peak("dp_" + lab, r)
                    print("  %s: k=%d %s estados=%s terminos=%d polimeros=%d (%.1fs; carga %.1fs, pesos %.1fs)" % (
                        lab, r["k"], r["modo"], r.get("estados", r.get("estados_procesados")), r["terminos"],
                        r["polimeros"], r["segundos"], r["segundos_carga"], r["segundos_pesos"]), flush=True)
                guard("trabajos de la DP")
        finally:
            mon.stop()
    return R


# ------------------------------------------------------------------ controles con los resultados de la DP
def job_summary(r):
    keys = ("k", "modo", "completo", "grupos", "estados", "estados_procesados", "terminos", "polimeros",
            "polimeros_por_ancla_max", "pendientes_max", "anclas_malas", "entradas", "quitadas", "bytes_cargados",
            "bytes_por_estado_pendiente", "segundos", "segundos_carga", "segundos_pesos", "segundos_dp", "error")
    return {x: r[x] for x in keys if x in r}


class NoMedida(Exception):
    pass


def need(R, lab):
    r = R.get(lab)
    if r is None or "error" in r:
        raise NoMedida(lab + (": " + r["error"] if r else ": sin resultado"))
    return r


def eval_fen(C, R, RB4, I, m, ctl):
    name = "TIRA_%dx6" % m
    tag = "fen_" + name
    cat = C[("fen", m)]
    cx = cat["cx"]
    o = I[tag]
    o["celdas"] = {}
    for k in (3, 4):
        labs = ["fen_%d_k3" % m] if k == 3 else ["fen_%d_k4_%d" % (m, j) for j in range(FEN_SPLIT[m])]
        rows = [need(R, lb) for lb in labs]
        gk = RB4[name]["gas_por_k"][str(k)]
        exp = [gk["polimeros"], gk["estados"], gk["terminos"]]
        got = [[r["polimeros"], r["estados"], r["terminos"]] for r in rows]
        check("C-B4 %s k=%d: polimeros, estados y terminos = B4 %s (%d trabajos)" % (tag, k, exp, len(rows)),
              all(g == exp for g in got), got)
        fails, maxw, ncell = [], F(0), 0
        for r in rows:
            for key, vals in zip(r["claves"], r["valores"]):
                p, q = (F(x) for x in key.split("|"))
                B = odds_fen(cx, p, q)[1]
                iv = cell_intervals(B, vals)
                sealed = RB4[name]["celdas"][key]["por_k"][str(k)]
                row = o["celdas"].setdefault(key, {})
                for nm, sk in (("eps_comp", "eps_comp"), ("L_comp", "L_comp"), ("P_0", "P_0_cota_inferior")):
                    ex = b4.frac(sealed[sk]["exacto"])
                    lo, hi = iv[nm]
                    if not lo <= ex <= hi:
                        fails.append([key, k, nm])
                    w = relw(lo, hi)
                    maxw = max(maxw, w)
                    row["%s_k%d" % (nm, k)] = [dec(lo, False), dec(hi, True)]
                ncell += 1
        check("C-B4 %s k=%d: los intervalos de eps_comp, L_comp y P_0 de las %d celdas contienen la fraccion "
              "sellada en B4" % (tag, k, ncell), not fails and ncell == 12, fails)
        check("C-anch %s k=%d: anchura relativa de los %d intervalos <= 1e-30 (maxima %.3e)" % (
            tag, k, 3 * ncell, float(maxw)), maxw <= REL_W, float(maxw))
        ctl.setdefault("C-B4", {})["%s_k%d" % (tag, k)] = {"estructura": got[0], "celdas": ncell,
                                                           "anchura_relativa_max": float(maxw)}


def sum_interval(cx, k, vals):
    S = 1 << PREC
    Ul, Uh = vals
    NS = len(cx["nclass"])
    DN = t2.D ** NS
    X0_lo = (DN * Ul) >> PREC
    X0_hi = -((-(DN * Uh)) >> PREC)
    width = F(Uh - Ul, S)
    return t2.criterion(cx["nclass"], k, X0_hi), t2.criterion(cx["nclass"], k, X0_lo), width


def eval_sum(tag, cx, k, r, ref, want, ctl, key):
    c_hi, c_lo, width = sum_interval(cx, k, r["valores"][0])
    verdicts = [c_hi["veredicto"], c_lo["veredicto"]]
    check("C-T20 %s k=%d %s: veredicto de la regla de suma en los dos extremos = %s" % (tag, k, key, want),
          verdicts == [want, want], verdicts)
    dm = max(abs(c["margen_log10"] - ref["margen_log10"]) for c in (c_hi, c_lo))
    dp = max(abs(c["peso_orden_siguiente_aprox"] / ref["peso_orden_siguiente_aprox"] - 1) for c in (c_hi, c_lo))
    check("C-T20 %s k=%d %s: margen (dif. abs. %.2e <= 1e-6) y peso del orden siguiente (dif. rel. %.2e <= 1e-9) "
          "= T2-0" % (tag, k, key, dm, dp), dm <= TOL_MARGEN and dp <= TOL_PESO,
          {"margen": [c_hi["margen_log10"], c_lo["margen_log10"]], "T2-0": ref["margen_log10"]})
    lim = SUM_W * F(1, t2.D ** k)
    check("C-anch %s k=%d %s: anchura de U <= 1e-12 x^k (%.3e)" % (tag, k, key, float(width / lim) * 1e-12),
          width <= lim, float(width))
    ctl.setdefault("C-T20", {})["%s_%s_k%d" % (tag, key, k)] = {
        "veredicto": want, "margen_log10": [c_hi["margen_log10"], c_lo["margen_log10"]],
        "peso_orden_siguiente_aprox": [c_hi["peso_orden_siguiente_aprox"], c_lo["peso_orden_siguiente_aprox"]],
        "dif_margen": dm, "dif_rel_peso": dp, "anchura_U": float(width)}


def eval_circ(C, R, RT, I, m, ctl):
    name = "TIRA_%dx6" % m
    tag = "circ_" + name
    cat = C[("circ", m)]
    cx = cat["cx"]
    kT = 4 if m in (6, 12) else 3
    rt = RT[name]
    for k in range(1, kT + 1):
        lab = "circ_%d_est_k%d" % (m, k) if k < kT else "circ_%d_suma_k%d" % (m, k)
        r = need(R, lab)
        e = rt["estructura"][str(k)]
        exp = [e["polimeros"], e["estados"], e["terminos"]]
        got = [r["polimeros"], r["estados"], r["terminos"]]
        check("C-T20 %s k=%d: polimeros, estados y terminos = T2-0 %s" % (tag, k, exp), got == exp, got)
    eval_sum(tag, cx, kT, need(R, "circ_%d_suma_k%d" % (m, kT)), rt["regla_de_suma"], "COMPLETO", ctl, "suma")
    if m == 6:
        r = need(R, "circ_6_negativo_k4")
        check("C-T20 %s: el negativo quita exactamente una decoracion" % tag, r["quitadas"] == 1, r["quitadas"])
        eval_sum(tag, cx, 4, r, rt["regla_de_suma_negativo"]["criterio"], "INCOMPLETO", ctl, "negativo")
        S = 1 << PREC
        for j, (rc, rm) in enumerate(CELLS_EXACT):
            r = need(R, "circ_6_exacto_%d" % j)
            vals, ex = r["valores"][0], r["exactos"][0]
            ok, maxw = True, F(0)
            for (lo, hi), x in zip(((vals[0], vals[1]), (vals[2], vals[3]), (vals[4], vals[5])), ex):
                ok = ok and F(lo, S) <= x <= F(hi, S)
                maxw = max(maxw, relw(F(lo, S), F(hi, S)))
            check("C-exacto %s k=3 celda %s|%s: U, E y O exactos dentro de los intervalos de punto fijo" % (
                tag, rc, rm), ok)
            check("C-anch %s exacto %s|%s: anchura relativa <= 1e-30 (%.3e)" % (tag, rc, rm, float(maxw)),
                  maxw <= REL_W, float(maxw))
            ctl.setdefault("C-exacto", {})["%s|%s" % (rc, rm)] = {"estados": r["estados"], "terminos": r["terminos"],
                                                                 "anchura_relativa_max": float(maxw)}
    if m in (6, 12):
        r = need(R, "circ_%d_dom_k4" % m)
        e = rt["estructura"]["4"]
        check("C-dom %s: estructura del paso de celda = T2-0" % tag,
              [r["polimeros"], r["estados"], r["terminos"]] == [e["polimeros"], e["estados"], e["terminos"]])
        B = odds_circ(cx, *CELL_DOM)[1]
        L = cell_intervals(B, r["valores"][0])["L_comp"]
        dom = F(rt["P_sil_dominante"]["%s|%s" % CELL_DOM]["P_sil_dominante"])
        rlo, rhi = L[0] / dom, L[1] / dom
        check("C-dom %s k=4 celda %s|%s: L_comp / termino dominante de T2-0 en [1 - 1e-12, 1 + 1e-2] (%.12f, %.12f)"
              % (tag, CELL_DOM[0], CELL_DOM[1], float(rlo), float(rhi)), DOM_LO <= rlo and rhi <= DOM_HI,
              [float(rlo), float(rhi)])
        w = relw(*L)
        check("C-anch %s dom: anchura relativa de L_comp <= 1e-30 (%.3e)" % (tag, float(w)), w <= REL_W, float(w))
        ctl.setdefault("C-dom", {})[tag] = {"cociente": [dec(rlo, False), dec(rhi, True)], "anchura_relativa":
                                            float(w)}


def measures(C, R, I, out):
    med = {}
    for m in (24, 48):
        cat = C.get(("circ", m))
        if cat is None:
            continue
        row = {"catalogo": I["circ_TIRA_%dx6" % m]["catalogo"]}
        for k in ((4,) if m == 24 else (3, 4)):
            r = R.get("circ_%d_est_k%d" % (m, k))
            row["estructura_k%d" % k] = job_summary(r) if r else None
        med["TIRA_%dx6" % m] = row
    cal = {}
    for m in (12, 48):
        r = R.get("circ_%d_cal" % m)
        if r is None or "error" in r:
            cal[str(m)] = None if r is None else {"error": r["error"]}
            continue
        cal[str(m)] = {"grupos": r["grupos"], "terminos": r["terminos"], "entradas": r["entradas"],
                       "bytes_cargados": r["bytes_cargados"], "segundos_dp": r["segundos_dp"],
                       "segundos_pesos": r["segundos_pesos"], "segundos_carga": r["segundos_carga"],
                       "t_celda_s_por_termino": r["segundos_dp"] / max(1, r["terminos"]),
                       "t_peso_s_por_entrada": r["segundos_pesos"] / max(1, r["entradas"]),
                       "t_carga_s_por_byte": r["segundos_carga"] / max(1, r["bytes_cargados"]),
                       "bytes_por_estado_pendiente": r["bytes_por_estado_pendiente"], "pico_bytes": r["peak"],
                       "nota": "valores de la celda descartados; solo tiempos"}
    med["calibracion_tiempo"] = cal
    out["medidas"] = med


def decide_m2(C, R, out):
    try:
        terms = {}
        for m in (6, 12):
            terms[m] = need(R, "circ_%d_suma_k4" % m)["terminos"]
        for m in (24, 48):
            terms[m] = need(R, "circ_%d_est_k4" % m)["terminos"]
        cal = {12: out["medidas"]["calibracion_tiempo"]["12"], 48: out["medidas"]["calibracion_tiempo"]["48"]}
        if cal[12] is None or cal[48] is None or "error" in cal[12] or "error" in cal[48]:
            raise NoMedida("calibracion de tiempo")
        est48 = need(R, "circ_48_est_k4")
    except NoMedida as ex:
        out["decision_D_m2"] = {"decision": "NO_MEDIDA", "motivo": str(ex)}
        return
    main_peak = b4.mem_now(peak=True) or 0
    mem_proc = cal[48]["pico_bytes"] + est48["pendientes_max"] * cal[48]["bytes_por_estado_pendiente"]
    procs = max(0, min(PROCS, int((MEM_TOTAL - main_peak) // max(1, mem_proc))))
    rows, cpu, cat_s, disk = {}, 0.0, 0.0, 0
    for m in (6, 12, 24, 48):
        c = cal[12] if m <= 12 else cal[48]
        cat = C[("circ", m)]
        ent = cat["entradas_por_k"]["4"]
        s = terms[m] * c["t_celda_s_por_termino"] + ent * c["t_peso_s_por_entrada"] + \
            cat["bytes"] * c["t_carga_s_por_byte"]
        rows[str(m)] = {"terminos_4": terms[m], "entradas_4": ent, "bytes": cat["bytes"], "segundos_por_celda": s}
        cpu += s * T21_CELLS
        cat_s += cat["segundos_catalogo"] + cat["segundos_fusion"]
        disk += cat["bytes"]
    T = (cpu / procs if procs else float("inf")) + cat_s
    if mem_proc > MEM_DP or procs < 1 or disk > DISK_MAX:
        d = "NO_CABE"
    elif T <= HOURS_OK * 3600:
        d = "CABE"
    else:
        d = "CABE_CON_MAS_TIEMPO"
    out["decision_D_m2"] = {
        "decision": d, "horas_previstas": T / 3600, "procesos": procs, "memoria_prevista_por_proceso_bytes": mem_proc,
        "memoria_principal_bytes": main_peak, "disco_bytes": disk, "segundos_catalogos": cat_s, "por_m": rows,
        "informativo_un_trabajo_de_m48_horas": rows["48"]["segundos_por_celda"] / 3600,
        "regla": "suma_m (terminos_4 t_celda + entradas_4 t_peso + bytes t_carga) x 8 / procesos + catalogos; "
                 "CABE si <= 6 h y memoria y disco dentro de los techos"}


# ------------------------------------------------------------------ corrida
def run(out, args):
    ce, rect, b1 = b4.load_all(with_b1=True)
    ctx = mp.get_context("spawn")
    RB4 = json.load(open(B4D + "RESULTADO_B4.json", encoding="utf-8"))["instancias"]
    RT = json.load(open(T20 + "RESULTADO_T2_0.json", encoding="utf-8"))["instancias"]
    grid8 = [(rc, rm) for rc in ce.GRID_P for rm in (rc, F(0))]
    grid12 = [(p, q) for p in ce.GRID_P for q in (p, p / 2, F(0))]
    I = out["instancias"]
    wroot = args.trabajo
    os.makedirs(wroot, exist_ok=False)
    fen_ms = (6, 12) if args.prueba else (6, 12, 24, 48)
    circ_ms = (6,) if args.prueba else (6, 12, 24, 48)
    C = {}
    t = time.perf_counter()
    for m in fen_ms:
        C[("fen", m)] = catalogue_instance(ctx, ce, rect, m, True, wroot, RB4, RT, grid8, I)
    for m in circ_ms:
        C[("circ", m)] = catalogue_instance(ctx, ce, rect, m, False, wroot, RB4, RT, grid8, I)
    TIMES["catalogos"] = round(time.perf_counter() - t, 3)
    man = {cat["nombre"]: cat["manifiesto"] for cat in C.values()}
    mpath = os.path.join(wroot, "MANIFIESTO.json")
    with open(mpath, "w", encoding="ascii", newline="\n") as fh:
        fh.write(json.dumps(man, ensure_ascii=True, indent=1, sort_keys=True) + "\n")
    out["trabajo"] = {"manifiesto": mpath.replace("\\", "/"), "manifiesto_sha256": sha_file(mpath),
                      "bytes": sum(cat["bytes"] for cat in C.values()), "pico_disco_bytes": DISK[1]}
    t = time.perf_counter()
    R = run_jobs(ctx, make_specs(C, grid12, args.prueba))
    TIMES["dp"] = round(time.perf_counter() - t, 3)
    out["trabajos"] = {lb: job_summary(r) for lb, r in sorted(R.items())}
    ctl = out["controles"] = {}
    nomed = []
    bad = sum(r.get("anclas_malas", 0) for r in R.values() if "error" not in r)
    check("C-anc: en los %d trabajos, ningun hijo con ancla <= la de su padre" % len(R), bad == 0, bad)
    for m in fen_ms:
        try:
            eval_fen(C, R, RB4, I, m, ctl)
        except NoMedida as ex:
            nomed.append("C-B4 m=%d: %s" % (m, ex))
    for m in [x for x in circ_ms if x <= 24]:
        try:
            eval_circ(C, R, RT, I, m, ctl)
        except NoMedida as ex:
            nomed.append("C-T20/C-exacto/C-dom m=%d: %s" % (m, ex))
    ctl["controles_no_medidos"] = nomed
    for x in nomed:
        LIMITS.append("control NO_MEDIDO: " + x)
    out["decision_D_m1"] = "VALIDADO" if not nomed else "NO_VALIDADO (controles no medidos)"
    if args.prueba:
        r = R.get("circ_6_cal")
        out["medidas"] = {"prueba_calibracion_m6": job_summary(r) if r else None}
        out["decision_D_m2"] = {"decision": "NO_MEDIDA", "motivo": "prueba de desarrollo"}
    else:
        measures(C, R, I, out)
        decide_m2(C, R, out)
    out["lectura_global"] = "D-m1 %s; D-m2 %s" % (out["decision_D_m1"], out["decision_D_m2"]["decision"])
    print("LECTURA:", out["lectura_global"], flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--recursos", required=True)
    ap.add_argument("--trabajo", required=True)
    ap.add_argument("--prueba", action="store_true")
    args = ap.parse_args()
    out = {"unit": "T2-1M", "plan": PLAN, "plan_sha256": SHA[PLAN], "referencias_sha256": SHA,
           "programa_sha256": sha_file(os.path.abspath(__file__)), "prueba": bool(args.prueba), "P_bits": PREC,
           "grupos": NGROUPS, "D_regla_de_suma": str(t2.D), "instancias": {}}
    status = "FAIL"
    try:
        run(out, args)
        status = "PASS" if not LIMITS else "PASS_CON_NO_MEDIDAS"
    except b4.CheckFail as ex:
        out["detenido_en"] = "control fallido: " + str(ex)
        out["decision_D_m1"] = "NO_VALIDADO (control fallido)"
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
                             "tiempos": dict(TIMES, total=round(time.perf_counter() - T0, 3)),
                             "pico_memoria_principal_bytes": b4.mem_now(peak=True),
                             "pico_memoria_total_bytes_por_etapa": MEMPEAK,
                             "pico_memoria_procesos_bytes": b4.WPEAK, "pico_disco_bytes": DISK[1],
                             "llamadas_totales": b4.CALLS[0]}, ensure_ascii=True, indent=1, sort_keys=True) + "\n")
    print("%s controles=%d llamadas=%d t=%.1fs" % (status, len(b4.CHECKS), b4.CALLS[0], time.perf_counter() - T0),
          flush=True)


if __name__ == "__main__":
    main()
