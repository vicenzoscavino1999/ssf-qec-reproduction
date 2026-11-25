#!/usr/bin/env python3
"""B4 (PLAN_B4.txt, sellado en DEC-048, con ENMIENDA_1.txt, sellada en DEC-049): escalado en N a orden
silencioso fijo con tiras HGP(camino_m, camino_6), m = 6, 12, 24, 48, k = 1..4, 12 celdas.

Importa, tras comprobar su SHA-256 y sin cambios: control_exacto.py de D1-CE (decoder, logicos, truncacion
ingenua, show/err), b1_exacto_k4.py de B1 (DP del gas: prep_gas, anchor_masks, explore, eval_fly,
finish_cell) y calibracion/rect.py de esta unidad (construccion rectangular). Piezas nuevas: catalogo por la
regla de crecimiento de B2 con las mismas pruebas que cat_task de B1, recuento P3 del rival ingenuo,
controles C+, C- y V, compuerta, metricas y decisiones del plan y de la enmienda.

Uso (desde D:\\paper ia cuantica\\qec_revision_matematica\\05_COMPROBACIONES):
  python -I -B B4_TIRAS_ESCALADO_N_20261113/b4_tiras.py --output <nuevo.json> --recursos <nuevo.json> --tmp <carpeta>
  (--prueba: prueba de desarrollo con C+, C- y TIRA_8x6 de calibracion; no toca V ni E)
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

sys.set_int_max_str_digits(0)
T0 = time.perf_counter()
UDIR = "B4_TIRAS_ESCALADO_N_20261113/"
CE_PATH = "D1_CONTROL_EXACTO_Y_METODOLOGIA_20261113/control_exacto.py"
B1_PATH = "B1_K4_EXACTO_S1_20261113/b1_exacto_k4.py"
RECT_PATH = UDIR + "calibracion/rect.py"
RB1 = "B1_K4_EXACTO_S1_20261113/RESULTADO_B1.json"
RCE = "D1_CONTROL_EXACTO_Y_METODOLOGIA_20261113/RESULTADO_CONTROL.json"
FASEA = "D1_VIABILIDAD_ALCANCE_CERTIFICADO_20261112/FASE_A.json"
PLAN = UDIR + "PLAN_B4.txt"
ENM = UDIR + "ENMIENDA_1.txt"
CALOUT = UDIR + "calibracion/SALIDA_calib.txt"
SHA = {
    CE_PATH: "c003050183a14364750a33784c7c815c0c07b5814f552a8b52a14d142db3306b",
    B1_PATH: "3f0c3796181d2b4c7be752f2d44fd9fe489a9f3f71667073a3d21936c19dbc92",
    RECT_PATH: "358442f672f616c63861080ad32051d23cf470d8c46ef981aa50d2666b1a6254",
    RB1: "da8cd31148c91cfec5200af3e1ec4166c6bd010f7e5111657adb8905795653a8",
    RCE: "334d6485968c914c299bbc7e5362cf49f4b975b740d47d78bccfaec47255024d",
    FASEA: "e68435581f661be720cc03be80aae62261b8810fdc0aa09ea5b3af934b87de2c",
    PLAN: "08b3bc64e59c772674aff1f3b97a5944d0e15562787d833125fa69f3d7991238",
    ENM: "d36d53316e9262f87df454ff7462cae8103b7838960807160faeb7f553b5de25",
    CALOUT: "9e31253ca2c04b6dbdf60e6993022a8f0dd90ed26c0d8cceae5a1b1341374c6e",
}
EXH8 = {1: 123, 2: 1736, 3: 30471, 4: 572227}     # exhaustivo de TIRA_8x6 (calibracion/SALIDA_calib.txt)
NSHORT = 6
M_GRID = (6, 12, 24, 48)
KMAX = 4
GIB = 1 << 30
SEG_MAX = 5 * 3600.0
LLAMADAS_MAX = 80_000_000
MEM_TOTAL = 18 * GIB
MEM_MAIN = 9 * GIB
MEM_W = GIB // 2
PROCS_CAT = 14
W_GAS_MAX = 3
GATE_STATES = 5_000_000
GATE_TERMS = 2_500_000_000
CAL_TERMS = 20_000_000
SMALL_DP_S = 300.0
TARGETS = (F(1, 20), F(1, 10), F(1, 4))
CALLS = [0]
CHECKS = []
LIMITS = []
TIMES = {}
WPEAK = {}
MODS = {}


class Stop(Exception):
    """Techo de recursos."""


class CheckFail(Exception):
    """Control fallido."""


def sha_file(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def load_mod(path, name):
    if sha_file(path) != SHA[path]:
        raise SystemExit("SHA-256 distinto del sellado: " + path)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_all(with_b1=True):
    if "ce" not in MODS:
        MODS["ce"] = load_mod(CE_PATH, "control_exacto")
        MODS["rect"] = load_mod(RECT_PATH, "rect_b4")
    if with_b1 and "b1" not in MODS:
        MODS["b1"] = load_mod(B1_PATH, "b1_exacto_k4")
        MODS["b1"].CE[0] = MODS["ce"]
    return MODS["ce"], MODS["rect"], MODS.get("b1")


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


def guard(where):
    if time.perf_counter() - T0 > SEG_MAX:
        raise Stop("techo de tiempo alcanzado en " + where)
    m = mem_now()
    if m is not None and m > MEM_MAIN:
        raise Stop("techo de memoria del proceso principal alcanzado en " + where)


def check(desc, ok, detail=None, verbose=True):
    CHECKS.append({"control": desc, "pasa": bool(ok), "detalle": detail})
    if verbose or not ok:
        print("  [%s] %s" % ("PASS" if ok else "FALLA", desc), flush=True)
    if not ok:
        raise CheckFail(desc)


def frac(s):
    a, b = s.split("/")
    return F(int(a), int(b))


def bits(v):
    while v:
        b = v & -v
        yield b.bit_length() - 1
        v ^= b


def group_of(reps):
    g = []
    for msk in range(1, 1 << len(reps)):
        u = 0
        for i in range(len(reps)):
            if msk >> i & 1:
                u ^= reps[i]
        g.append(u)
    return g


def make_inst(ce, rect, kind, m):
    if kind == "tira":
        r1, a1 = ce.path_seed(m)
        r2, a2 = ce.path_seed(NSHORT)
        return rect.build_rect(ce, "TIRA_%dx%d" % (m, NSHORT), r1, a1, r2, a2)
    rows, a = ce.regular_seed_34()
    return rect.build_rect(ce, "S2_34_6x8", rows, a, rows, a)


def mats_sha(inst):
    return hashlib.sha256(json.dumps({"hx": inst["hx"], "hz": inst["hz"]}, sort_keys=True).encode()).hexdigest()


# ------------------------------------------------------------------ procesos de trabajo del catalogo
_W = {}


def w_init(kind, m):
    sys.set_int_max_str_digits(0)
    ce, rect, _ = load_all(with_b1=False)
    inst = make_inst(ce, rect, kind, m)
    dec, RES = ce.make_decoder(inst)
    _W.update(ce=ce, inst=inst, dec=dec, RES=RES, group=group_of(ce.logical_reps(inst)))


def classify(sm):
    """Mismas pruebas que cat_task de B1: pasa, silenciosa, caracteres logicos y validez de Y."""
    ce, inst, dec, group = _W["ce"], _W["inst"], _W["dec"], _W["group"]
    psi, cols, bz, adj, N = inst["psi"], inst["cols"], inst["bz"], inst["adjK"], inst["N"]
    s = 0
    for j in bits(sm):
        s ^= psi[j]
    cor, touched = dec(s)
    E = sm & ((1 << N) - 1)
    na = E.bit_count()
    key = (na, sm.bit_count() - na)
    r = E ^ cor
    passes = ce.lin(cols, r) == 0
    silent, chis, bad = False, None, 0
    if passes:
        silent = not ce.inside(r, bz)
        chis = tuple(-1 if (u & r).bit_count() & 1 else 1 for u in group)
        bad = int(silent != any(c < 0 for c in chis))
    Y = sm | touched
    seen = 0
    todo = Y & -Y
    while todo:
        bb = todo & -todo
        todo ^= bb
        seen |= bb
        todo |= adj[bb.bit_length() - 1] & Y & ~seen
    return key, passes, silent, chis, (Y if seen == Y else None), bad


def w_grow(chunk):
    out, bad = [], 0
    for sm in chunk:
        if len(_W["RES"]) > 300000:
            _W["RES"].clear()
        key, passes, silent, chis, Y, b = classify(sm)
        bad += b
        if Y is not None:
            out.append((sm, Y, key, passes, silent, chis))
    _W["RES"].clear()
    m = mem_now()
    if m is not None and m > MEM_W:
        raise RuntimeError("techo de memoria de un proceso del catalogo (%d bytes)" % m)
    return {"out": out, "bad": bad, "calls": len(chunk), "peak": mem_now(peak=True), "pid": os.getpid()}


def w_exh(task):
    f, firsts = task
    NS = _W["inst"]["NS"]
    valid, sil, hist, bad = [], {}, 0, 0
    for first in firsts:
        for rest in combinations(range(first + 1, NS), f - 1):
            sm = 1 << first
            for i in rest:
                sm |= 1 << i
            hist += 1
            if len(_W["RES"]) > 300000:
                _W["RES"].clear()
            key, passes, silent, chis, Y, b = classify(sm)
            bad += b
            if silent:
                sil[key] = sil.get(key, 0) + 1
            if Y is not None:
                valid.append(sm)
        m = mem_now()
        if m is not None and m > MEM_W:
            raise RuntimeError("techo de memoria de un proceso del exhaustivo (%d bytes)" % m)
    _W["RES"].clear()
    return {"valid": valid, "sil": sil, "hist": hist, "bad": bad, "peak": mem_now(peak=True), "pid": os.getpid()}


def note_peak(stage, res):
    pk = WPEAK.setdefault(stage, {})
    pk[str(res["pid"])] = max(pk.get(str(res["pid"]), 0), res["peak"] or 0)


def split_firsts(NS, f, nbins):
    items = sorted(((math.comb(NS - 1 - i, f - 1), i) for i in range(NS - f + 1)), reverse=True)
    bins = [[0, []] for _ in range(nbins)]
    for w, i in items:
        b = min(bins, key=lambda x: x[0])
        b[0] += w
        b[1].append(i)
    return [(f, sorted(b[1])) for b in bins if b[1]]


# ------------------------------------------------------------------ catalogo por crecimiento (regla de B2)
def grow(pool, inst, name, kmax, keep_last_set, nchar):
    NS, adj = inst["NS"], inst["adjK"]
    rec, POLY, silc, info = {}, {}, {}, {}
    last_set = set() if keep_last_set else None
    bad = 0
    for f in range(1, kmax + 1):
        t = time.perf_counter()
        if f == 1:
            cand = [1 << i for i in range(NS)]
        else:
            cs = set()
            for A, val in rec[f - 1].items():
                nb = val[0]
                for v in bits(val[0]):
                    nb |= adj[v]
                nb &= ~A
                for x in bits(nb):
                    cs.add(A | (1 << x))
            cand = sorted(cs)
            del cs
        if CALLS[0] + len(cand) > LLAMADAS_MAX:
            raise Stop("techo de llamadas: catalogo de %s con f=%d no cabe" % (name, f))
        n = max(1, len(cand) // (8 * PROCS_CAT))
        chunks = [cand[i:i + n] for i in range(0, len(cand), n)]
        got = {} if f < kmax else None
        sf = silc.setdefault(f, {})
        nval, calls = 0, 0
        try:
            for res in pool.imap_unordered(w_grow, chunks):
                calls += res["calls"]
                bad += res["bad"]
                note_peak("catalogo_" + name, res)
                for sm, Y, key, passes, silent, chis in res["out"]:
                    nval += 1
                    slot = POLY.get(Y)
                    if slot is None:
                        slot = {}
                        POLY[Y] = slot
                    ent = slot.get(key)
                    if ent is None:
                        ent = [0] * nchar
                        slot[key] = ent
                    ent[0] += 1
                    if passes:
                        ent[1] += 1
                        for i, c in enumerate(chis):
                            ent[2 + i] += c
                    if silent:
                        sf[key] = sf.get(key, 0) + 1
                    if got is not None:
                        got[sm] = (Y, key, passes, silent)
                    elif last_set is not None:
                        last_set.add(sm)
                res["out"] = None
                guard("catalogo de %s f=%d" % (name, f))
        except (Stop, CheckFail):
            raise
        except Exception as ex:
            raise Stop("catalogo de %s f=%d: %s" % (name, f, ex))
        CALLS[0] += calls
        check("%s f=%d: cada candidato distinto se decodifica una vez (llamadas = |C_f| = %d)"
              % (name, f, len(cand)), calls == len(cand))
        if f == 1:
            check("%s f=1: se prueban los NS = %d sitios" % (name, NS), len(cand) == NS)
        if got is not None:
            rec[f] = got
        info[str(f)] = {"candidatos": len(cand), "validas": nval, "C(NS,f)": math.comb(NS, f),
                        "segundos": round(time.perf_counter() - t, 3)}
        print("  %s f=%d: %d candidatos, %d validas (%.1fs)" % (name, f, len(cand), nval,
                                                                 time.perf_counter() - t), flush=True)
        del cand, chunks
    check("I4 %s: silenciosa <=> algun caracter logico vale -1 (discrepancias = 0)" % name, bad == 0, bad)
    return rec, POLY, silc, info, last_set


def exhaustive(pool, inst, name, kmax):
    NS = inst["NS"]
    valid, sil, bad = {}, {}, 0
    for f in range(1, kmax + 1):
        t = time.perf_counter()
        c = math.comb(NS, f)
        if CALLS[0] + c > LLAMADAS_MAX:
            raise Stop("techo de llamadas: exhaustivo de %s con f=%d no cabe" % (name, f))
        vs, hist = set(), 0
        try:
            for res in pool.imap_unordered(w_exh, split_firsts(NS, f, 4 * PROCS_CAT)):
                vs.update(res["valid"])
                hist += res["hist"]
                bad += res["bad"]
                for kk, v in res["sil"].items():
                    sil[kk] = sil.get(kk, 0) + v
                note_peak("exhaustivo_" + name, res)
                guard("exhaustivo de %s f=%d" % (name, f))
        except (Stop, CheckFail):
            raise
        except Exception as ex:
            raise Stop("exhaustivo de %s f=%d: %s" % (name, f, ex))
        CALLS[0] += hist
        check("%s exhaustivo f=%d: historias = C(NS,f) = %d" % (name, f, c), hist == c)
        valid[f] = vs
        print("  %s exhaustivo f=%d: %d historias, %d validas (%.1fs)" % (name, f, hist, len(vs),
                                                                        time.perf_counter() - t), flush=True)
    check("I4 %s exhaustivo: silenciosa <=> algun caracter logico vale -1" % name, bad == 0, bad)
    return valid, sil


# ------------------------------------------------------------------ recuento P3 del rival ingenuo
def p3_count(inst, rec, silc, k):
    adj = inst["adjK"]
    nonempty = [f for f in silc if f <= k and silc[f]]
    fs = min(nonempty) if nonempty else None
    info = {"f_s": fs, "k": k}
    check("%s P3 aplicable con k=%d: ninguna decoracion valida silenciosa con < k-1 fallos (f_s = %s)"
          % (inst["name"], k, fs), fs is None or fs >= k - 1, fs)
    naive = {}
    for f in range(1, k + 1):
        for key, c in silc.get(f, {}).items():
            naive[key] = naive.get(key, 0) + c
    pairs = 0
    ben = [(Y, key) for sm, (Y, key, passes, silent) in rec[1].items() if passes and not silent]
    info["singletes_que_pasan_con_clase_trivial"] = len(ben)
    info["singletes_que_no_pasan"] = sum(1 for v in rec[1].values() if not v[2])
    if fs is not None and fs + 1 <= k:
        for sm, (Y, key, passes, silent) in rec[fs].items():
            if not silent:
                continue
            nb = Y
            for v in bits(Y):
                nb |= adj[v]
            for Yx, kx in ben:
                if Yx & nb == 0:
                    kk = (key[0] + kx[0], key[1] + kx[1])
                    naive[kk] = naive.get(kk, 0) + 1
                    pairs += 1
    info["pares_compatibles"] = pairs
    return naive, info


def keyed(d):
    return {"%d,%d" % kk: v for kk, v in sorted(d.items())}


# ------------------------------------------------------------------ DP del gas (B1) y compuerta
def gas_init(path):
    sys.set_int_max_str_digits(0)
    load_all(with_b1=True)
    with open(path, "rb") as fh:
        G, sts = pickle.load(fh)
    _W["G"], _W["sts"] = G, sts
    _W["index"] = {J: i for i, J in enumerate(sts)}
    _W["A"] = MODS["b1"].anchor_masks(G)


def gas_task(task):
    p, q, NV, mem_limit = task
    t = time.perf_counter()
    X0, XT, T = MODS["b1"].eval_fly(_W["G"], _W["A"], _W["sts"], _W["index"], p, q, NV, mem_limit=mem_limit)
    return {"p": p, "q": q, "X0": X0, "XT": XT, "T": T, "segundos": round(time.perf_counter() - t, 3),
            "peak": mem_now(peak=True), "pid": os.getpid()}


def gate(b1, G, sts, terms, NV, name):
    """Formula de la compuerta de B1 (gate_check) con los techos de B4."""
    ns = len(sts)
    info = {"estados": ns, "terminos": terms, "limite_estados": GATE_STATES, "limite_terminos": GATE_TERMS}
    if ns > GATE_STATES or terms > GATE_TERMS:
        info.update(pasa=False, motivo="estados o terminos por encima de la compuerta")
        return 0, None, info
    A = b1.anchor_masks(G)
    index = {J: i for i, J in enumerate(sts)}
    p = q = F(1, 10 ** 6)
    cal = b1.eval_fly(G, A, sts, index, p, q, NV, cal_terms=min(CAL_TERMS, terms))
    t_term = cal["segundos"] / max(cal["terminos"], 1)
    del A, index
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
                                 "segundos": round(cal["segundos"], 3), "segundos_por_termino": t_term},
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
    print("  compuerta %s k=4: pasa; %d procesos; previsto %.0f s" % (name, nproc, projected), flush=True)
    return nproc, mem_limit, info


def eval_main(b1, G, sts, cells, NV, klog, name, k):
    out, secs = {}, {}
    A = b1.anchor_masks(G)
    index = {J: i for i, J in enumerate(sts)}
    for (p, q) in cells:
        t = time.perf_counter()
        X0, XT, T = b1.eval_fly(G, A, sts, index, p, q, NV)
        kk = "%s|%s" % (p, q)
        out[kk] = b1.finish_cell(p, q, X0, XT, T, klog)
        secs[kk] = round(time.perf_counter() - t, 3)
        guard("DP de %s k=%d" % (name, k))
    return out, secs


def dump_gas(G, sts, tmpdir, name, k):
    path = os.path.join(tmpdir, "b4_gas_%s_k%d.pkl" % (name, k))
    with open(path, "wb") as fh:
        pickle.dump((G, sts), fh, protocol=pickle.HIGHEST_PROTOCOL)
    return path


def eval_pool(ctx, b1, path, cells, NV, klog, nproc, mem_limit, name, k):
    out, secs = {}, {}
    try:
        with ctx.Pool(nproc, initializer=gas_init, initargs=(path,)) as pool:
            try:
                for res in pool.imap_unordered(gas_task, [(p, q, NV, mem_limit) for (p, q) in cells]):
                    kk = "%s|%s" % (res["p"], res["q"])
                    out[kk] = b1.finish_cell(res["p"], res["q"], res["X0"], res["XT"], res["T"], klog)
                    secs[kk] = res["segundos"]
                    note_peak("dp_%s_k%d" % (name, k), res)
                    print("    %s k=%d celda %s lista (%.0fs)" % (name, k, kk, res["segundos"]), flush=True)
                    guard("DP de %s k=%d" % (name, k))
            except (Stop, CheckFail):
                raise
            except Exception as ex:
                raise Stop("DP de %s k=%d: %s" % (name, k, ex))
    finally:
        os.remove(path)
    return out, secs


# ------------------------------------------------------------------ una instancia completa
def run_instance(ctx, ce, b1, rect, kind, m, kmax, cells, order_kind, tmpdir, want_exh=False,
                 gate_k4=False, ref=None, ref_kind=None):
    t_inst = time.perf_counter()
    inst = make_inst(ce, rect, kind, m)
    name = inst["name"]
    NS, N, M = inst["NS"], inst["N"], inst["M"]
    reps = ce.logical_reps(inst)
    check("%s: %d representantes logicos = dimension logica %d" % (name, len(reps), inst["klog"]),
          len(reps) == inst["klog"])
    klog = len(reps)
    NV = 1 + (1 << klog)
    if kind == "tira":
        check("%s: N = 11m-5, M = 5m, NS = 16m-5 y un qubit logico" % name,
              (N, M, NS, klog) == (11 * m - 5, 5 * m, 16 * m - 5, 1), [N, M, NS, klog])
    deg = max(x.bit_count() for x in inst["adjK"])
    oi = {"N": N, "M": M, "sitios": NS, "qubits_logicos": klog, "matrices_sha256": mats_sha(inst),
          "grado_max_K": deg, "grado_sobre_NS": round(deg / NS, 4), "orden_dp": order_kind}
    print("%s: %d sitios, grado maximo de K %d (%.3f de NS)" % (name, NS, deg, deg / NS), flush=True)
    if ref_kind == "B1":
        check("%s matrices iguales a S1_L6 de B1" % name, oi["matrices_sha256"] == ref["matrices_sha256"])
    elif ref_kind == "CE":
        check("%s matrices iguales a FASE_A de D1" % name, oi["matrices_sha256"] == ref["fase_a_sha"])

    # ---- catalogo por crecimiento
    t = time.perf_counter()
    c0 = CALLS[0]
    with ctx.Pool(PROCS_CAT, initializer=w_init, initargs=(kind, m)) as pool:
        rec, POLY, silc, cinfo, last_set = grow(pool, inst, name, kmax, want_exh, NV)
        TIMES[name + "_catalogo"] = round(time.perf_counter() - t, 3)
        calls_growth = {k: sum(cinfo[str(f)]["candidatos"] for f in range(1, k + 1)) for k in range(1, kmax + 1)}
        n3 = sum(1 for slot in POLY.values() if any(kk[0] + kk[1] <= 3 for kk in slot))
        naive_sil, pinfo = p3_count(inst, rec, silc, kmax)
        oi["catalogo"] = {"por_f": cinfo, "polimeros_distintos": len(POLY), "polimeros_distintos_f_le_3": n3,
                          "llamadas": CALLS[0] - c0,
                          "silenciosas_conexas_por_f": {str(f): keyed(v) for f, v in sorted(silc.items())},
                          "P3": dict(pinfo, silenciosas_todas=keyed(naive_sil))}
        print("  %s catalogo listo (%.1fs); polimeros %d; P3 %s" % (name, TIMES[name + "_catalogo"], len(POLY),
                                                                    keyed(naive_sil)), flush=True)
        # ---- validacion exhaustiva (V)
        if want_exh:
            t = time.perf_counter()
            valid_x, sil_x = exhaustive(pool, inst, name, kmax)
            TIMES[name + "_exhaustivo"] = round(time.perf_counter() - t, 3)
            got = {f: (set(rec[f]) if f < kmax else last_set) for f in range(1, kmax + 1)}
            miss = {str(f): len(valid_x[f] - got[f]) for f in got}
            extra = {str(f): len(got[f] - valid_x[f]) for f in got}
            oi["validacion_exhaustiva"] = {"validas_por_f": {str(f): len(valid_x[f]) for f in valid_x},
                                           "faltan": miss, "sobran": extra,
                                           "silenciosas_todas_exhaustivo": keyed(sil_x)}
            check("V %s: catalogo por crecimiento = exhaustivo, conjunto a conjunto, f = 1..%d" % (name, kmax),
                  all(valid_x[f] == got[f] for f in got), {"faltan": miss, "sobran": extra})
            check("V %s: P3 = recuento exhaustivo de historias silenciosas por (na, nb)" % name, sil_x == naive_sil,
                  {"exhaustivo": keyed(sil_x), "P3": keyed(naive_sil)})
            oi["_validas_exh"] = {f: len(valid_x[f]) for f in valid_x}
            del valid_x, got
    del last_set

    if ref_kind == "B1":
        rc = ref["catalogo"]
        check("C+ %s: decoraciones validas por f iguales a B1 (S1_L6)" % name,
              all(cinfo[f]["validas"] == rc["por_f"][f]["decoraciones_validas"] for f in cinfo))
        check("C+ %s: polimeros distintos %d y con f <= 3 %d (B1)" % (name, rc["polimeros_distintos"],
                                                                     rc["polimeros_distintos_f_le_3"]),
              len(POLY) == rc["polimeros_distintos"] and n3 == rc["polimeros_distintos_f_le_3"])
        check("C+ %s: silenciosas conexas por f y (na, nb) iguales a B1" % name,
              all(keyed(silc.get(int(f), {})) == rc["por_f"][f]["silenciosas_conexas_por_na_nb"] for f in rc["por_f"]))
        check("C+ %s: P3 = silenciosas todas de B1 %s" % (name, rc["silenciosas_todas"]),
              keyed(naive_sil) == rc["silenciosas_todas"])
    elif ref_kind == "CE":
        rc = ref["catalogo"]
        check("C- %s: decoraciones validas por f iguales a D1-CE" % name,
              all(cinfo[f]["validas"] == rc["por_f"][f]["decoraciones_validas"] for f in cinfo))
        check("C- %s: polimeros distintos %d (D1-CE)" % (name, rc["polimeros_distintos"]),
              len(POLY) == rc["polimeros_distintos"])
        check("C- %s: P3 = silenciosas todas de D1-CE %s" % (name, rc["silenciosas_todas"]),
              keyed(naive_sil) == rc["silenciosas_todas"])
    del rec

    # ---- orden del barrido
    coord = inst["coord"]
    order_prod = sorted(range(NS), key=lambda v: (coord[v][0], coord[v][1], v))
    order_A = ce.phase_a_order(inst) if (order_kind == "fase_A" or ref_kind == "B1") else None
    order = order_prod if order_kind == "producto" else order_A
    oi["frontera"] = ce.frontier(order, inst["adjK"])
    if ref_kind is not None:
        t = time.perf_counter()
        for k in range(1, kmax + 1):
            G = b1.prep_gas(inst, POLY, order_A, k)
            sts, terms = b1.explore(G, b1.anchor_masks(G))
            gk = ref["gas_por_k"][str(k)]
            exp = (gk["polimeros"], gk["estados"], gk["terminos"])
            check("%s %s k=%d orden de la fase A: polimeros, estados y terminos = %s (sellados)"
                  % ("C+" if ref_kind == "B1" else "C-", name, k, list(exp)),
                  (len(G["polys"]), len(sts), terms) == exp, [len(G["polys"]), len(sts), terms])
            del G, sts
        TIMES[name + "_estructura_fase_A"] = round(time.perf_counter() - t, 3)

    # ---- DP por k
    polyholder = [POLY]
    del POLY
    table, gas_k, evaluated = {}, {}, {}
    for k in range(1, kmax + 1):
        t = time.perf_counter()
        G = b1.prep_gas(inst, polyholder[0], order, k)
        if k == kmax:
            polyholder[0] = None
        sts, terms = b1.explore(G, b1.anchor_masks(G))
        TIMES["%s_estructura_k%d" % (name, k)] = round(time.perf_counter() - t, 3)
        row = {"polimeros": len(G["polys"]), "estados": len(sts), "terminos": terms}
        print("  %s k=%d: %d polimeros, %d estados, %d terminos (%.1fs)" % (name, k, len(G["polys"]), len(sts), terms,
                                                                          time.perf_counter() - t), flush=True)
        guard("estructura de %s k=%d" % (name, k))
        comp = None
        t2 = time.perf_counter()
        try:
            if gate_k4 and k == kmax:
                nproc, mem_limit, ginfo = gate(b1, G, sts, terms, NV, name)
                row["compuerta"] = ginfo
                if nproc == 0:
                    row["no_evaluada"] = "compuerta: " + ginfo.get("motivo", "")
                    LIMITS.append("%s k=%d NO_EVALUADA por la compuerta: %s" % (name, k, ginfo.get("motivo")))
            else:
                nproc, mem_limit = 1, None
            if nproc == 1:
                comp, secs = eval_main(b1, G, sts, cells, NV, klog, name, k)
            elif nproc > 1:
                path = dump_gas(G, sts, tmpdir, name, k)
                G = sts = None
                comp, secs = eval_pool(ctx, b1, path, cells, NV, klog, nproc, mem_limit, name, k)
        except Stop as ex:
            comp = None
            row["no_evaluada"] = "limite de recursos: " + str(ex)
            LIMITS.append("%s k=%d NO_EVALUADA: %s" % (name, k, ex))
            print("  %s k=%d NO_EVALUADA: %s" % (name, k, ex), flush=True)
        if comp is not None:
            TIMES["%s_gas_k%d" % (name, k)] = round(time.perf_counter() - t2, 3)
            TIMES["%s_gas_k%d_por_celda" % (name, k)] = secs
        evaluated[k] = comp is not None
        cellrow = {}
        for (p, q) in cells:
            kk = "%s|%s" % (p, q)
            Li, eps_i = ce.naive(inst, p, q, naive_sil, k)
            c = {"eps_ing": eps_i, "L_ing": Li, "eps_comp": None, "L_comp": None, "P0": None}
            if comp is not None:
                eps_c, Lc, P0 = comp[kk]
                c.update(eps_comp=eps_c, L_comp=Lc, P0=P0)
                check("I2 %s k=%d %s: 0 <= eps_comp <= eps_ing" % (name, k, kk), 0 <= eps_c <= eps_i,
                      verbose=False)
                check("I5 %s k=%d %s: L_ing <= L_comp" % (name, k, kk), Li <= Lc, verbose=False)
                ec, ei = ce.err(Lc, eps_c), ce.err(Li, eps_i)
                if ec is not None and ei is not None:
                    check("I6 %s k=%d %s: e_comp <= e_ing" % (name, k, kk), ec <= ei, verbose=False)
            cellrow[kk] = c
        if ref_kind == "B1" and comp is not None:
            rcl = ref["celdas"]
            same = all(cellrow[kk][a] == frac(rcl[kk]["por_k"][str(k)][b]["exacto"]) for kk in cellrow
                       for a, b in (("eps_comp", "eps_comp"), ("L_comp", "L_comp"), ("eps_ing", "eps_ing"),
                                    ("L_ing", "L_ing"), ("P0", "P_0_cota_inferior")))
            check("C+ %s k=%d orden de producto: eps_comp, L_comp, P_0, eps_ing y L_ing identicos a B1 (12 celdas)"
                  % (name, k), same)
        if ref_kind == "CE" and comp is not None:
            rcl = ref["celdas"]
            same = all(cellrow[kk][a] == frac(rcl[kk]["por_k"][str(k)][b]["exacto"]) for kk in cellrow
                       for a, b in (("eps_comp", "eps_comp"), ("L_comp", "L_comp"), ("eps_ing", "eps_ing"),
                                    ("L_ing", "L_ing"), ("P0", "P_0_cota_inferior")))
            check("C- %s k=%d: eps_comp, L_comp, P_0, eps_ing y L_ing identicos a D1-CE (12 celdas)" % (name, k),
                  same)
        table[k] = cellrow
        gas_k[str(k)] = row
        del G, sts
        guard("fin de k=%d en %s" % (k, name))
    for (p, q) in cells:
        kk = "%s|%s" % (p, q)
        ks = [k for k in range(1, kmax + 1) if evaluated[k]]
        Ls = [table[k][kk]["L_comp"] for k in ks]
        Es = [table[k][kk]["eps_comp"] for k in ks]
        Li = [table[k][kk]["L_ing"] for k in range(1, kmax + 1)]
        Ei = [table[k][kk]["eps_ing"] for k in range(1, kmax + 1)]
        check("I3 %s %s: L_comp y L_ing no decrecen; eps_comp y eps_ing no crecen con k" % (name, kk),
              all(Ls[j] <= Ls[j + 1] for j in range(len(Ls) - 1))
              and all(Es[j] >= Es[j + 1] for j in range(len(Es) - 1))
              and all(Li[j] <= Li[j + 1] for j in range(kmax - 1))
              and all(Ei[j] >= Ei[j + 1] for j in range(kmax - 1)), verbose=False)
    nck = sum(1 for c in CHECKS if c["control"].startswith(("I2 " + name, "I5 " + name, "I6 " + name,
                                                             "I3 " + name)))
    print("  %s: %d controles I2, I3, I5 e I6 por celda, todos PASS" % (name, nck), flush=True)
    oi["gas_por_k"] = gas_k
    oi["evaluada_por_k"] = {str(k): evaluated[k] for k in evaluated}
    oi["llamadas_por_k"] = {str(k): {"catalogo_por_crecimiento": calls_growth[k],
                                     "ingenua_convencional_sum_C(NS,f)": sum(math.comb(NS, f) for f in range(1, k + 1))}
                            for k in range(1, kmax + 1)}
    oi["celdas"] = summarize(ce, inst, table, cells, kmax, evaluated, calls_growth)
    oi["muestreo_por_subconjuntos"] = subset_sampling(inst, naive_sil, calls_growth, kmax)
    oi["segundos"] = round(time.perf_counter() - t_inst, 3)
    oi["_tabla"] = table
    oi["_eval"] = evaluated
    print("%s ok %.1fs llamadas=%d" % (name, oi["segundos"], CALLS[0]), flush=True)
    return oi


# ------------------------------------------------------------------ metricas
def summarize(ce, inst, table, cells, kmax, evaluated, calls_growth):
    NS = inst["NS"]
    res = {}
    for (p, q) in cells:
        kk = "%s|%s" % (p, q)
        rows = {}
        for k in range(1, kmax + 1):
            c = table[k][kk]
            ok = evaluated[k]
            ec = ce.err(c["L_comp"], c["eps_comp"]) if ok else None
            ei = ce.err(c["L_ing"], c["eps_ing"])
            R = c["eps_ing"] / c["eps_comp"] if ok and c["eps_comp"] > 0 else None
            rows[str(k)] = {"eps_comp": ce.show(c["eps_comp"]) if ok else None, "eps_ing": ce.show(c["eps_ing"]),
                            "L_comp": ce.show(c["L_comp"]) if ok else None, "L_ing": ce.show(c["L_ing"]),
                            "P_0_cota_inferior": ce.show(c["P0"]) if ok else None,
                            "e_comp": ce.show(ec) if ec is not None else None,
                            "e_ing": ce.show(ei) if ei is not None else None,
                            "R": ce.show(R) if R is not None else None, "componentes_evaluadas": ok,
                            "llamadas_catalogo_por_crecimiento": calls_growth[k],
                            "llamadas_ingenua_convencional": sum(math.comb(NS, f) for f in range(1, k + 1)),
                            "_e_comp": ec, "_e_ing": ei, "_R": R}
        _, epsn = ce.naive(inst, p, q, {}, kmax + 1)
        Lk = table[kmax][kk]["L_ing"]
        boundn = epsn / (2 * Lk + epsn) if Lk > 0 else None
        kst = {}
        for tg in TARGETS:
            kc = next((k for k in range(1, kmax + 1) if rows[str(k)]["_e_comp"] is not None
                       and rows[str(k)]["_e_comp"] <= tg), None)
            ki = next((k for k in range(1, kmax + 1) if rows[str(k)]["_e_ing"] is not None
                       and rows[str(k)]["_e_ing"] <= tg), None)
            if ki is None:
                ki = (kmax + 1) if (boundn is not None and boundn <= tg) else ">%d" % (kmax + 1)
            kst[str(tg)] = {"componentes": kc if kc is not None else ">%d" % kmax, "ingenua": ki}
        res[kk] = {"p": str(p), "q": str(q), "por_k": rows, "eps_ing_%d" % (kmax + 1): ce.show(epsn),
                   "cota_e_ing_%d_informativa" % (kmax + 1): ce.show(boundn) if boundn is not None else None,
                   "k_estrella": kst}
    return res


def subset_sampling(inst, naive_sil, calls_growth, kmax):
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
    return {"w_estrella": w, "n_w_estrella": n, "C_NS_w": C, "pi": float(pi), "muestras_10_por_ciento_95": nss,
            "llamadas_catalogo_k_igual_w": sum(math.comb(NS, f) for f in range(1, w + 1)),
            "llamadas_crecimiento_hasta_w": calls_growth[w], "llamadas_crecimiento_hasta_k4": calls_growth[kmax],
            "nota": "solo q = p; el muestreo no da certificado ni acota la cola"}


def strip_private(oi):
    for cell in oi.get("celdas", {}).values():
        for row in cell["por_k"].values():
            for kk in [x for x in row if x.startswith("_")]:
                del row[kk]
    for kk in [x for x in oi if x.startswith("_")]:
        del oi[kk]


def m8_category(vals):
    n_rel = sum(1 for v in vals if v >= 10)
    n_nul = sum(1 for v in vals if v < 2)
    cat = "RELEVANTE" if n_rel == len(vals) else ("NULA" if n_nul > 0 else "ZONA_GRIS")
    return {"categoria": cat, "celdas_mayor_o_igual_10": n_rel, "celdas_menor_que_2": n_nul,
            "celdas_entre_2_y_10": len(vals) - n_rel - n_nul}


def decide(ce, S, keys, s2):
    """S: {m: oi} de las tiras evaluadas; s2: oi de S2."""
    out = {}
    R4 = {m: ({kk: S[m]["celdas"][kk]["por_k"]["4"]["_R"] for kk in keys} if m in S and S[m]["_eval"][4] else None)
          for m in M_GRID}
    t1 = {}
    for m in (12, 24, 48):
        if R4[m] is None:
            t1[str(m)] = {"evaluada": False}
            continue
        vals = list(R4[m].values())
        row = m8_category(vals)
        row.update(evaluada=True, R4_min=ce.dec(min(vals)), R4_max=ce.dec(max(vals)))
        t1[str(m)] = row
    steps = {}
    all_grow = True
    for a, b in ((6, 12), (12, 24), (24, 48)):
        if R4[a] is None or R4[b] is None:
            steps["%d->%d" % (a, b)] = "NO_DECIDIDO"
            all_grow = None if all_grow is not False else False
            continue
        n = sum(1 for kk in keys if R4[b][kk] > R4[a][kk])
        steps["%d->%d" % (a, b)] = {"celdas_que_crecen": n, "celdas": len(keys)}
        if n < len(keys):
            all_grow = False
    t1["pasos"] = steps
    t1["lectura"] = ("CRECE_CON_N" if all_grow is True else ("NO_CRECE" if all_grow is False else "NO_DECIDIDO"))
    slopes = {}
    for a, b in ((6, 12), (12, 24), (24, 48)):
        if R4[a] is not None and R4[b] is not None:
            vals = [math.log2(ce.dec(R4[b][kk] / R4[a][kk])) for kk in keys]
            slopes["%d->%d" % (a, b)] = {"min": round(min(vals), 4), "max": round(max(vals), 4),
                                         "por_celda": {kk: round(v, 4) for kk, v in zip(keys, vals)}}
    t1["pendientes_log2_informativas"] = slopes
    out["D_t1"] = t1
    t2 = {}
    if 12 in S and S[12]["_eval"][4]:
        base = {kk: S[12]["celdas"][kk]["por_k"]["4"] for kk in keys}
        stable = True
        for m in (24, 48):
            if m not in S or not S[m]["_eval"][4]:
                t2[str(m)] = "NO_DECIDIDO"
                stable = None if stable is not False else False
                continue
            def ratio(a, b):
                return None if a is None or b is None or b == 0 else a / b

            rc = [ratio(S[m]["celdas"][kk]["por_k"]["4"]["_e_comp"], base[kk]["_e_comp"]) for kk in keys]
            ri = [ratio(S[m]["celdas"][kk]["por_k"]["4"]["_e_ing"], base[kk]["_e_ing"]) for kk in keys]
            if any(v is None for v in rc):
                t2[str(m)] = {"e_comp_indefinido_en_alguna_celda": True}
                stable = None if stable is not False else False
                continue
            n = sum(1 for v in rc if v >= 2)
            rid = [v for v in ri if v is not None]
            t2[str(m)] = {"celdas_con_cociente_e_comp_mayor_o_igual_2": n, "cociente_e_comp_min": ce.dec(min(rc)),
                          "cociente_e_comp_max": ce.dec(max(rc)),
                          "cociente_e_ing_min": ce.dec(min(rid)) if rid else None,
                          "cociente_e_ing_max": ce.dec(max(rid)) if rid else None,
                          "cociente_e_comp_por_celda": {kk: ce.dec(v) for kk, v in zip(keys, rc)},
                          "cociente_e_ing_por_celda": {kk: (ce.dec(v) if v is not None else None)
                                                       for kk, v in zip(keys, ri)}}
            if n > 0:
                stable = False
        t2["lectura"] = "ESTABLE" if stable is True else ("CRECE" if stable is False else "NO_DECIDIDO")
    else:
        t2["lectura"] = "NO_DECIDIDO"
    out["D_t2"] = t2
    t3 = {}
    for m in M_GRID:
        if m not in S or not S[m]["_eval"][4]:
            t3[str(m)] = {"evaluada": False}
            continue
        row = {"evaluada": True}
        for tg in TARGETS:
            adv = [kk for kk in keys if S[m]["celdas"][kk]["por_k"]["4"]["_e_comp"] is not None
                   and S[m]["celdas"][kk]["por_k"]["4"]["_e_comp"] <= tg
                   and (S[m]["celdas"][kk]["por_k"]["4"]["_e_ing"] is None
                        or S[m]["celdas"][kk]["por_k"]["4"]["_e_ing"] > tg)]
            row["celdas_con_ventaja_a_igual_k_%s" % tg] = len(adv)
            row["lista_%s" % tg] = adv
        t3[str(m)] = row
    out["D_t3_ventaja_a_igual_k"] = t3
    r3 = [s2["celdas"][kk]["por_k"]["3"]["_R"] for kk in keys]
    out["D_t4"] = {"R3_min": ce.dec(min(r3)), "R3_max": ce.dec(max(r3)),
                   "lectura": "PASA" if all(v < 2 for v in r3) else "FALLA"}
    e5 = {}
    for m in M_GRID:
        if m in S and S[m]["_eval"][4] and S[m]["_eval"][3]:
            e5[str(m)] = sum(1 for kk in keys if S[m]["celdas"][kk]["por_k"]["4"]["_R"]
                             > S[m]["celdas"][kk]["por_k"]["3"]["_R"])
    out["E5_celdas_R4_mayor_que_R3"] = e5
    out["E3_regimen_host"] = {str(m): {"grado_max_K": S[m]["grado_max_K"], "grado_sobre_NS": S[m]["grado_sobre_NS"]}
                              for m in M_GRID if m in S}
    return out


def affine_check(S):
    """Completitud por recuento en m = 24 y 48 bajo la hipotesis afin (informativa; no decide)."""
    if 12 not in S or "_validas_exh" not in S[12]:
        return {"evaluada": False}
    c12 = S[12]["_validas_exh"]
    out = {"evaluada": True, "hipotesis": "recuentos afines en m desde m = 8 (medido en m = 8, 10, 14; no demostrado)",
           "exhaustivo_m8": {str(f): EXH8[f] for f in EXH8}, "exhaustivo_m12": {str(f): c12[f] for f in c12}}
    for m in (24, 48):
        if m not in S:
            out[str(m)] = None
            continue
        row = {}
        for f in range(1, KMAX + 1):
            pred = F(c12[f]) + F(m - 12) * F(c12[f] - EXH8[f], 4)
            got = S[m]["catalogo"]["por_f"][str(f)]["validas"]
            row[str(f)] = {"prevista": str(pred), "crecimiento": got, "igual": pred == got}
        row["todas_iguales"] = all(row[str(f)]["igual"] for f in range(1, KMAX + 1))
        out[str(m)] = row
    return out


# ------------------------------------------------------------------ programa principal
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--recursos", required=True)
    ap.add_argument("--tmp", required=True)
    ap.add_argument("--prueba", action="store_true",
                    help="prueba de desarrollo: C+, C- y TIRA_8x6 (calibracion); no toca V ni E")
    args = ap.parse_args()
    for pth in (args.output, args.recursos):
        if os.path.exists(pth):
            raise SystemExit("La salida ya existe; use un nombre nuevo: " + pth)
    if not os.path.isdir(args.tmp):
        raise SystemExit("No existe la carpeta temporal: " + args.tmp)
    for pth in SHA:
        if sha_file(pth) != SHA[pth]:
            raise SystemExit("SHA-256 distinto del sellado: " + pth)
    out = {"unit": "B4", "plan": PLAN, "plan_sha256": SHA[PLAN], "enmienda": ENM, "enmienda_sha256": SHA[ENM],
           "programa_sha256": sha_file(os.path.abspath(__file__)), "referencias_sha256": SHA,
           "prueba": args.prueba, "instancias": {}}
    status = "PASS"
    try:
        run(out, args)
    except CheckFail as ex:
        status = "FAIL"
        out["detenido_en"] = "control fallido: " + str(ex)
    except Stop as ex:
        status = "LIMITE_DE_RECURSOS"
        out["detenido_en"] = str(ex)
    if status == "PASS" and LIMITS:
        status = "PASS_CON_NO_EVALUADAS"
    out["status"] = status
    out["limites_alcanzados"] = LIMITS
    out["checks"] = CHECKS
    out["llamadas_totales"] = CALLS[0]
    for oi in out["instancias"].values():
        strip_private(oi)
    with open(args.output, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(json.dumps(out, ensure_ascii=False, indent=1, sort_keys=True, default=str) + "\n")
    TIMES["total"] = round(time.perf_counter() - T0, 3)
    with open(args.recursos, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(json.dumps({"nota": "tiempos de reloj en segundos y memoria; no deterministas; no son resultado",
                             "tiempos": TIMES, "pico_memoria_principal_bytes": mem_now(peak=True),
                             "pico_memoria_procesos_bytes": WPEAK, "llamadas_totales": CALLS[0]},
                            ensure_ascii=False, indent=1, default=str) + "\n")
    print("%s controles=%d llamadas=%d t=%.1fs" % (status, len(CHECKS), CALLS[0], TIMES["total"]), flush=True)
    for x in LIMITS:
        print("LIMITE:", x, flush=True)
    if status not in ("PASS", "PASS_CON_NO_EVALUADAS"):
        print("DETENIDO EN:", out.get("detenido_en"), flush=True)
        sys.exit(1)


def run(out, args):
    ce, rect, b1 = load_all(with_b1=True)
    ctx = mp.get_context("spawn")
    for name, L in (("S1_L6", 6), ("S2_34_6x8", None)):
        if L is not None:
            rows, a = ce.path_seed(L)
        else:
            rows, a = ce.regular_seed_34()
        ref = ce.build(name, rows, a)
        new = rect.build_rect(ce, name, rows, a, rows, a)
        check("build_rect con h1 = h2 igual a build() en %s (%d campos)" % (name, len(ref)),
              all(new[x] == ref[x] for x in ref))
    cal = open(CALOUT, encoding="ascii").read()
    check("exhaustivo de TIRA_8x6 en la calibracion sellada: %s" % EXH8,
          "exhaustivo: catalogo igual True {1: 123, 2: 1736, 3: 30471, 4: 572227}" in cal)
    RB = json.load(open(RB1, encoding="utf-8"))["instancias"]["S1_L6"]
    RC = json.load(open(RCE, encoding="utf-8"))["instancias"]["S2_34_6x8"]
    RC["fase_a_sha"] = json.load(open(FASEA, encoding="utf-8"))["instancias"]["S2_34_6x8"]["matrices_sha256"]
    grid = [(p, q) for p in ce.GRID_P for q in (p, p / 2, F(0))]
    keys = ["%s|%s" % (p, q) for (p, q) in grid]
    S = {}

    print("== C+ (control positivo): TIRA_6x6 = S1 L=6 de B1", flush=True)
    S[6] = run_instance(ctx, ce, b1, rect, "tira", 6, KMAX, grid, "producto", args.tmp, ref=RB, ref_kind="B1")
    out["instancias"]["TIRA_6x6"] = S[6]

    print("== C- (control negativo): S2 de D1-CE, k <= 3", flush=True)
    s2 = run_instance(ctx, ce, b1, rect, "s2", None, 3, grid, "fase_A", args.tmp, ref=RC, ref_kind="CE")
    out["instancias"]["S2_34_6x8"] = s2
    r3 = [s2["celdas"][kk]["por_k"]["3"]["_R"] for kk in keys]
    d4 = "PASA" if all(v < 2 for v in r3) else "FALLA"
    print("  D-t4 (control negativo): %s; R_3 entre %.4f y %.4f" % (d4, ce.dec(min(r3)), ce.dec(max(r3))), flush=True)
    if d4 != "PASA":
        out["lectura_detenida"] = "D-t4 FALLA: R_3 >= 2 en S2; se detiene la lectura (M7)"
        return

    if args.prueba:
        print("== PRUEBA: TIRA_8x6 (calibracion) con exhaustivo, compuerta y DP en paralelo", flush=True)
        o8 = run_instance(ctx, ce, b1, rect, "tira", 8, KMAX, grid, "producto", args.tmp, want_exh=True, gate_k4=True)
        out["instancias"]["TIRA_8x6_prueba"] = o8
        check("prueba: exhaustivo de TIRA_8x6 igual al de la calibracion sellada",
              o8["_validas_exh"] == EXH8, o8["_validas_exh"])
        # solo para ejercitar el codigo de decisiones y completitud: TIRA_8x6 ocupa los lugares de m = 12 y 24
        out["prueba_codigo_decisiones"] = decide(ce, {6: S[6], 12: o8}, keys, s2)
        out["prueba_codigo_completitud"] = affine_check({12: o8, 24: o8})
        return

    print("== V y E: TIRA_12x6 (validacion exhaustiva)", flush=True)
    S[12] = run_instance(ctx, ce, b1, rect, "tira", 12, KMAX, grid, "producto", args.tmp, want_exh=True,
                         gate_k4=True)
    out["instancias"]["TIRA_12x6"] = S[12]
    for m in (24, 48):
        print("== E: TIRA_%dx6" % m, flush=True)
        try:
            S[m] = run_instance(ctx, ce, b1, rect, "tira", m, KMAX, grid, "producto", args.tmp, gate_k4=True)
            out["instancias"]["TIRA_%dx6" % m] = S[m]
        except Stop as ex:
            LIMITS.append("TIRA_%dx6 NO_EVALUADA: %s" % (m, ex))
            out["instancias"]["TIRA_%dx6" % m] = {"no_evaluada": "limite de recursos: " + str(ex)}
            print("  TIRA_%dx6 NO_EVALUADA: %s" % (m, ex), flush=True)
            break
    out["completitud_afin"] = affine_check(S)
    out["decisiones"] = decide(ce, S, keys, s2)
    d = out["decisiones"]
    t1, t2 = d["D_t1"], d["D_t2"]
    cats = ", ".join("m=%s %s" % (m, t1[m]["categoria"] if t1[m].get("evaluada") else "NO_EVALUADA")
                     for m in ("12", "24", "48"))
    adv = ", ".join("m=%s %s" % (m, ("%d/12" % d["D_t3_ventaja_a_igual_k"][m]["celdas_con_ventaja_a_igual_k_1/10"])
                                 if d["D_t3_ventaja_a_igual_k"][m].get("evaluada") else "NO_EVALUADA")
                    for m in ("6", "12", "24", "48"))
    out["lectura_global"] = ("D-t1 %s (%s); D-t2 %s; D-t3 ventaja a igual k al 10%%: %s; D-t4 %s; completitud afin "
                             "m=24 %s, m=48 %s" % (t1["lectura"], cats, t2["lectura"], adv, d["D_t4"]["lectura"],
                                                   (out["completitud_afin"].get("24") or {}).get("todas_iguales"),
                                                   (out["completitud_afin"].get("48") or {}).get("todas_iguales")))
    print("LECTURA:", out["lectura_global"], flush=True)


if __name__ == "__main__":
    main()
