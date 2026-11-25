#!/usr/bin/env python3
"""H6-U2: circuito de una ronda sobre PATH25 con ocho ubicaciones ruidosas.
Contrasta la ejecucion global (propagacion literal + decoder global) con la
recomposicion primitiva, el oraculo de primera salida y el gas de polimeros
del apendice C del manuscrito. Solo biblioteca estandar; aritmetica exacta.

Uso:  python -I -B verificar_h6_u2.py --output <nuevo.json> [--tiempos <nuevo.txt>]
      python -I -B verificar_h6_u2.py --output <nuevo.json> --parcial 0,3,6
--parcial restringe la lista de ubicaciones ruidosas (prueba de depuracion).
El plan esta en PLAN_UNIDAD.txt; este archivo no cambia instancia ni tasa."""
import argparse
import hashlib
import json
import os
import sys
import time
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import combinations, product

T0 = time.perf_counter()
WALL_CAP = 1150.0
MAX_ORACLE_SUPPORTS = 20000
HERE = os.path.dirname(os.path.abspath(__file__))
PILOT = os.path.join(os.path.dirname(HERE), "PILOTO_CONECTADO_PATH25_20261028")
FROZEN = {
    "preflight": os.path.join(PILOT, "inputs", "PREFLIGHT.json"),
    "table": os.path.join(PILOT, "attempt_01", "REFERENCE_TABLE.json"),
}

A_COLS, SEED = 4, (3, 6, 12)
B_ROWS = len(SEED)
N = A_COLS * A_COLS + B_ROWS * B_ROWS
M = A_COLS * B_ROWS
U0 = 4369
RATES = (F(1, 100), F(1, 1000))
# ubicaciones ruidosas fijadas en el plan (orden fijo)
ALL_LOCS = (
    ("dwait", 1, 0), ("dwait", 2, 0), ("dwait", 3, 0), ("dwait", 16, 0),
    ("dwait", 9, 0), ("dwait", 10, 0), ("cnot", 10, 8), ("read", 1),
)
LOC_NAMES = ("W1", "W2", "W3", "W16", "W9", "W10", "C10_8", "R1")

CHECKS = []
COUNT = {}
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


def need(cond, name, **detail):
    """Comprobacion dentro de un bucle: solo deja rastro si falla."""
    if not cond:
        CHECKS.append({"name": name, "pass": False, "detail": detail})
        raise Stop(name)


def guard():
    if time.perf_counter() - T0 > WALL_CAP:
        raise Stop("techo de tiempo de reloj alcanzado")


def tick(label, t):
    TIMES[label] = round(time.perf_counter() - t, 3)


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


def sha_file(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


# ---------------------------------------------------------------- codigo
def build_masks():
    a, H, b = A_COLS, SEED, B_ROWS
    hx = [sum(1 << (i * a + k) for k in bits(H[j])) + sum(1 << (a * a + l * b + j) for l in range(b) if H[l] >> i & 1)
          for i in range(a) for j in range(b)]
    hz = [sum(1 << (k * a + j) for k in bits(H[i])) + sum(1 << (a * a + i * b + l) for l in range(b) if H[l] >> j & 1)
          for i in range(b) for j in range(a)]
    cols = [sum(1 << c for c in range(M) if hx[c] >> v & 1) for v in range(N)]
    cand = set()
    for g in hz:
        f = g
        while f:
            if lin(cols, f):
                cand.add(f)
            f = (f - 1) & g
    actions = sorted(cand, key=lambda v: (bin(v).count("1"), v))
    return hx, hz, cols, actions, [lin(cols, f) for f in actions]


def build_sets():
    """Construccion independiente: productos de Kronecker y conjuntos."""
    a, b = A_COLS, B_ROWS
    h = [[(row >> i) & 1 for i in range(a)] for row in SEED]
    ht = [list(c) for c in zip(*h)]

    def kron(A, B):
        return [[x * y for x in ra for y in rb] for ra in A for rb in B]

    def eye(n):
        return [[int(i == j) for j in range(n)] for i in range(n)]

    mx = [x + y for x, y in zip(kron(eye(a), h), kron(ht, eye(b)))]
    mz = [x + y for x, y in zip(kron(h, eye(a)), kron(eye(b), ht))]
    hxs = [frozenset(i for i, v in enumerate(row) if v) for row in mx]
    hzs = [frozenset(i for i, v in enumerate(row) if v) for row in mz]
    colsets = [frozenset(j for j, row in enumerate(mx) if row[i]) for i in range(N)]

    def foot(fs):
        out = set()
        for q in fs:
            out.symmetric_difference_update(colsets[q])
        return frozenset(out)

    cand = set()
    for g in hzs:
        for k in range(1, len(g) + 1):
            for tup in combinations(sorted(g), k):
                fs = frozenset(tup)
                if foot(fs):
                    cand.add(fs)
    enc = lambda fs: sum(1 << i for i in fs)
    actions = sorted(cand, key=lambda fs: (len(fs), enc(fs)))
    return hxs, hzs, colsets, actions, [foot(fs) for fs in actions], enc


def rank_basis(vectors):
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


def in_span(v, piv):
    while v:
        k = v.bit_length() - 1
        if k not in piv:
            return False
        v ^= piv[k]
    return True


# ---------------------------------------------------------------- decoders
def transitions_bits(foot):
    """Motor 1: mascaras de bits y comparacion por productos cruzados."""
    widths = [bin(p).count("1") for p in foot]
    table = []
    nc = len(foot)
    for s in range(1 << M):
        best, num, den = -1, 0, 1
        for j in range(nc):
            d = widths[j]
            g = 2 * bin(s & foot[j]).count("1") - d
            if best < 0 or g * den > num * d:
                best, num, den = j, g, d
        if best < 0 or 2 * num < den:
            table.append((-1, s))
        else:
            table.append((best, s ^ foot[best]))
    COUNT["motor1_scores"] = (1 << M) * nc
    return table


def transitions_sets(footsets, enc):
    """Motor 2: conjuntos y fracciones; escrito aparte."""
    table = []
    half = F(1, 2)
    for s in range(1 << M):
        cur = frozenset(bits(s))
        win, sc = None, None
        for j, ft in enumerate(footsets):
            val = F(2 * len(cur & ft) - len(ft), len(ft))
            if sc is None or val > sc:
                win, sc = j, val
        if win is None or sc < half:
            table.append((-1, s))
        else:
            table.append((win, enc(cur.symmetric_difference(footsets[win]))))
    COUNT["motor2_scores"] = (1 << M) * len(footsets)
    return table


def resolve(table, actions):
    size = len(table)
    cor, end, traces, touched = [0] * size, list(range(size)), [()] * size, [0] * size
    for s in sorted(range(size), key=lambda x: bin(x).count("1")):
        j, nxt = table[s]
        if j >= 0:
            need(bin(nxt).count("1") < bin(s).count("1"), "descenso estricto del sindrome", s=s)
            cor[s] = actions[j] ^ cor[nxt]
            end[s] = end[nxt]
            traces[s] = (j,) + traces[nxt]
            touched[s] = actions[j] | touched[nxt]
    return cor, end, traces, touched


# ---------------------------------------------------------------- host K
def host_K(actionsets, footsets, colsets):
    nv = N + M
    psi = [colsets[v] for v in range(N)] + [frozenset((c,)) for c in range(M)]
    adj = [0] * nv
    for fs, phi in zip(actionsets, footsets):
        readers = [x for x in range(nv) if psi[x] & phi]
        rm = sum(1 << x for x in readers)
        am = sum(1 << y for y in fs)
        for x in readers:
            adj[x] |= rm | am          # (K1) y (K2)
        for y in fs:
            adj[y] |= rm               # (K2) simetrico
    for x in range(nv):
        adj[x] &= ~(1 << x)
    return adj


def connected_parts(Y, adj):
    out = []
    rem = Y
    while rem:
        seen = 0
        todo = rem & -rem
        while todo:
            b = todo & -todo
            todo ^= b
            seen |= b
            todo |= adj[b.bit_length() - 1] & rem & ~seen
        rem &= ~seen
        out.append(seen)
    return out


# ---------------------------------------------------------------- circuito
def greedy_schedule(hx):
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
    return layer


def x_of(p):
    return 1 if p in (1, 2) else 0


def z_of(p):
    return 1 if p in (2, 3) else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--tiempos")
    ap.add_argument("--parcial")
    args = ap.parse_args()
    if os.path.exists(args.output):
        raise SystemExit("La salida ya existe; use un nombre nuevo.")
    active = tuple(range(len(ALL_LOCS)))
    if args.parcial:
        active = tuple(sorted(set(int(x) for x in args.parcial.split(","))))
    full_run = active == tuple(range(len(ALL_LOCS)))
    result = {"unit": "H6-U2", "mode": "COMPLETA" if full_run else "PARCIAL",
              "active_locations": [LOC_NAMES[i] for i in active]}
    status = "PASS"
    try:
        run(active, full_run, result)
    except Stop as ex:
        status = "FAIL"
        result["stopped_at"] = str(ex)
    result["status"] = status
    result["checks_total"] = len(CHECKS)
    result["checks_failed"] = sum(1 for c in CHECKS if not c["pass"])
    result["checks"] = CHECKS
    result["counters"] = COUNT
    with open(args.output, "w", encoding="utf-8", newline="\n") as out:
        out.write(json.dumps(result, ensure_ascii=False, indent=1, sort_keys=True) + "\n")
    TIMES["total"] = round(time.perf_counter() - T0, 3)
    if args.tiempos:
        with open(args.tiempos, "w", encoding="utf-8", newline="\n") as out:
            out.write("Tiempos de reloj en segundos (no deterministas; no son resultado)\n")
            for k, v in TIMES.items():
                out.write("%-40s %10.3f\n" % (k, v))
    print("%s  controles=%d  fallidos=%d  t=%.1fs" % (status, result["checks_total"], result["checks_failed"],
                                                      TIMES["total"]))
    if status != "PASS":
        print("DETENIDO EN:", result.get("stopped_at"))
        sys.exit(1)


def run(active, full_run, result):
    # ------------------------------------------------ 0. codigo y decoders
    t = time.perf_counter()
    hx, hz, cols, actions, foot = build_masks()
    hxs, hzs, colsets, actionsets, footsets, enc = build_sets()
    check("A0 dos construcciones del codigo coinciden (H_X, H_Z, candidatos, huellas)",
          [enc(x) for x in hxs] == hx and [enc(x) for x in hzs] == hz
          and [enc(x) for x in actionsets] == actions and [enc(x) for x in footsets] == foot)
    check("A1 conmutacion H_X H_Z^T = 0", all(bin(x & z).count("1") % 2 == 0 for x in hx for z in hz))
    bx, bz = rank_basis(hx), rank_basis(hz)
    check("A2 n=25, 12 checks, rangos 12 y 12, k=1, 103 candidatos, 43 huellas",
          N == 25 and M == 12 and len(bx) == 12 and len(bz) == 12 and len(actions) == 103
          and len(set(foot)) == 43)
    check("A3 u_0=4369 es un representante logico protegido (conmuta con H_Z, fuera del espacio de filas de H_X)",
          all(bin(U0 & g).count("1") % 2 == 0 for g in hz) and not in_span(U0, bx))
    frozen_sha = {k: sha_file(p) for k, p in FROZEN.items()}
    pf = json.loads(open(FROZEN["preflight"], encoding="utf-8").read())
    tab = json.loads(open(FROZEN["table"], encoding="utf-8").read())
    check("A4 H_X y H_Z iguales a los del preflight congelado del piloto PATH25", pf["HX"] == hx and pf["HZ"] == hz)
    check("A5 candidatos y huellas iguales a la tabla congelada del piloto",
          tab["actions"] == actions and tab["footprints"] == foot)
    tick("codigo", t)

    t = time.perf_counter()
    T1 = transitions_bits(foot)
    T2 = transitions_sets(footsets, enc)
    check("B1 los dos motores dan la misma transicion en los 4096 sindromes", T1 == T2)
    check("B2 transiciones iguales a la tabla congelada del piloto en los 4096 sindromes",
          [list(x) for x in T1] == tab["transitions"])
    COR1, END1, TR1, TOUCH1 = resolve(T1, actions)
    COR2, END2, TR2, TOUCH2 = resolve(T2, actions)
    check("B3 correcciones, sindromes terminales y trazas iguales a la tabla congelada",
          COR1 == tab["correction"] and END1 == tab["terminal"] and [list(x) for x in TR1] == tab["trace"]
          and TOUCH1 == tab["touched"] and COR2 == COR1 and TR2 == TR1)
    maxlen = max(len(x) for x in TR1)
    multi = sum(1 for x in TR1 if len(x) >= 2)
    check("B4 traza maxima de 8 acciones y 4005 sindromes con dos o mas acciones (cifras del manuscrito)",
          maxlen == 8 and multi == 4005, {"max": maxlen, "dos_o_mas": multi})
    tick("decoders", t)

    # ------------------------------------------------ 1. host K
    t = time.perf_counter()
    adjK = host_K(actionsets, footsets, colsets)
    data_mask = (1 << N) - 1
    missing = sum(1 for x in range(N) for y in range(x + 1, N) if not adjK[x] >> y & 1)
    check("C1 host K igual al del preflight congelado", adjK == pf["host_adjacency"])
    check("C2 host K: 37 vertices, grado maximo 36, conexo; host de datos conexo con 111 pares no adyacentes",
          max(bin(x).count("1") for x in adjK) == 36 and len(connected_parts((1 << (N + M)) - 1, adjK)) == 1
          and len(connected_parts(data_mask, adjK)) == 1 and missing == 111, {"faltan": missing})
    tick("host", t)

    # ------------------------------------------------ 2. circuito
    t = time.perf_counter()
    layer = greedy_schedule(hx)
    r_w = max(bin(x).count("1") for x in hx)
    d_w = max(bin(x).count("1") for x in cols)
    L = r_w + d_w - 1
    h_nnz = sum(bin(x).count("1") for x in hx)
    data_locs = 0
    for v in range(N):
        data_locs += 2 + L          # espera inicial, L capas (puerta o espera), espera final
    anc_locs = 0
    for c in range(M):
        own = set(layer[(c, v)] for v in bits(hx[c]))
        anc_locs += 2 + (L - len(own))   # preparacion, lectura y esperas; las puertas ya se contaron
    n_loc = data_locs + anc_locs
    check("S0a r=4, d=2, L=5, 42 CNOT, colores dentro de L, N_loc=217 por recuento y por formula, 242 vertices",
          r_w == 4 and d_w == 2 and L == 5 and h_nnz == 42 and max(layer.values()) <= L
          and n_loc == 217 == (L + 2) * (N + M) - h_nnz and n_loc + N == 242,
          {"L": L, "N_loc": n_loc, "capas_usadas": max(layer.values())})
    check("S0b t(check 8, dato 10)=1 y t(check 7, dato 10)=2 (prediccion del plan)",
          layer[(8, 10)] == 1 and layer[(7, 10)] == 2, {"t_8_10": layer[(8, 10)], "t_7_10": layer[(7, 10)]})

    def Tmask(v, tt):
        return sum(1 << c for c in bits(cols[v]) if layer[(c, v)] > tt)

    locs = [ALL_LOCS[i] for i in active]
    names = [LOC_NAMES[i] for i in active]
    NL = len(locs)
    LOCMASK = (1 << NL) - 1
    alph = [4 if l[0] == "dwait" else 16 if l[0] == "cnot" else 8 for l in locs]
    total = 1
    for x in alph:
        total *= x
    need(total <= 524288, "techo de asignaciones", total=total)

    # estados estructurales y respuestas por la formula del apendice C
    state_of, n_states, ERESP, SRESP, base, loc_time = [], [], [], [], [], []
    for l in locs:
        if l[0] == "dwait":
            v, tt = l[1], l[2]
            state_of.append([0, 2, 2, 1])
            n_states.append(3)
            ERESP.append([0, 0, 1 << v])
            SRESP.append([0, 0, Tmask(v, tt)])
            base.append(v)
            loc_time.append(tt)
        elif l[0] == "cnot":
            v, c = l[1], l[2]
            tt = layer[(c, v)]
            so = []
            for q in range(16):
                a_, b_ = x_of(q >> 2), x_of(q & 3)
                so.append(0 if q == 0 else 1 + a_ + 2 * b_)
            state_of.append(so)
            n_states.append(5)
            ERESP.append([0] + [a_ << v for b_ in (0, 1) for a_ in (0, 1)])
            SRESP.append([0] + [(Tmask(v, tt) if a_ else 0) ^ ((1 << c) if b_ else 0) for b_ in (0, 1) for a_ in (0, 1)])
            base.append(v)
            loc_time.append(tt)
        else:
            c = l[1]
            state_of.append([0 if r == 0 else 1 + (x_of(r & 3) ^ (r >> 2)) for r in range(8)])
            n_states.append(3)
            ERESP.append([0, 0, 0])
            SRESP.append([0, 0, 1 << c])
            base.append(N + c)
            loc_time.append(L + 2)

    # programa del propagador literal: recorre el calendario puerta a puerta
    prog = []
    for k, l in enumerate(locs):
        if l[0] == "dwait" and l[2] == 0:
            prog.append((1, k, 1 << l[1]))
    for tt in range(1, L + 1):
        gates = tuple((1 << v, 1 << c) for (c, v), lt in sorted(layer.items()) if lt == tt)
        prog.append((0, gates))
        for k, l in enumerate(locs):
            if l[0] == "cnot" and layer[(l[2], l[1])] == tt:
                prog.append((2, k, 1 << l[1], 1 << l[2]))
            if l[0] == "dwait" and l[2] == tt:
                prog.append((1, k, 1 << l[1]))
    for k, l in enumerate(locs):
        if l[0] == "dwait" and l[2] == L + 1:
            prog.append((1, k, 1 << l[1]))
    reads = [(k, 1 << l[1]) for k, l in enumerate(locs) if l[0] == "read"]
    check("S0c el programa literal contiene las 42 puertas",
          sum(len(st[1]) for st in prog if st[0] == 0) == 42)

    def literal(cats):
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
                    s ^= cb          # Pauli antes de medir
                if r >> 2:
                    s ^= cb          # flip clasico despues
        return xd, s

    # host primitivo restringido a las ubicaciones activas y los 25 vertices de accion
    NV = NL + N
    vbase = base + list(range(N))
    ADJH = [0] * NV
    for x in range(NV):
        for y in range(NV):
            if x != y and (vbase[x] == vbase[y] or adjK[vbase[x]] >> vbase[y] & 1):
                ADJH[x] |= 1 << y
    psi_loc = [cols[bv] if bv < N else 1 << (bv - N) for bv in base]
    tick("circuito", t)

    # ------------------------------------------------ 3. leyes de probabilidad
    def cat_prob(kind, cat, e):
        if kind == "dwait":
            return 1 - e if cat == 0 else e / 3
        if kind == "cnot":
            return 1 - e if cat == 0 else e / 15
        return (1 - e if cat & 3 == 0 else e / 3) * (e if cat >> 2 else 1 - e)

    class_count, class_w = [], {e: [] for e in RATES}
    for k, l in enumerate(locs):
        cc = [0] * n_states[k]
        for cat in range(alph[k]):
            cc[state_of[k][cat]] += 1
        class_count.append(cc)
        for e in RATES:
            ww = [F(0)] * n_states[k]
            for cat in range(alph[k]):
                ww[state_of[k][cat]] += cat_prob(l[0], cat, e)
            class_w[e].append([x / ww[0] for x in ww])   # cociente con la categoria base
    Bc = {e: F(1) for e in RATES}
    for l in locs:
        for e in RATES:
            Bc[e] *= cat_prob(l[0], 0, e)

    # tipo de asignacion: (esperas no base, CNOT no base, Pauli de lectura no trivial, flip)
    n_wait = sum(1 for l in locs if l[0] == "dwait")
    sig_of = []
    for k, l in enumerate(locs):
        if l[0] == "dwait":
            sig_of.append([0, 1, 1, 1])
        elif l[0] == "cnot":
            sig_of.append([0] + [n_wait + 1] * 15)
        else:
            sig_of.append([((n_wait + 1) * 2) * ((1 if r & 3 else 0) + 2 * (r >> 2)) for r in range(8)])
    NSIG = (n_wait + 1) * 2 * 4

    def sig_prob(sig, e):
        nw = sig % (n_wait + 1)
        cn = (sig // (n_wait + 1)) % 2
        rp = (sig // ((n_wait + 1) * 2)) % 2
        rb = (sig // ((n_wait + 1) * 4)) % 2
        pr = (e / 3) ** nw * (1 - e) ** (n_wait - nw)
        if any(l[0] == "cnot" for l in locs):
            pr *= e / 15 if cn else 1 - e
        if any(l[0] == "read" for l in locs):
            pr *= (e / 3 if rp else 1 - e) * (e if rb else 1 - e)
        return pr

    # ------------------------------------------------ 4. R2: clases estructurales
    t = time.perf_counter()
    half = F(1, 2)

    def score2(s, j):
        ft = footsets[j]
        return F(2 * len(frozenset(bits(s)) & ft) - len(ft), len(ft))

    radix = []
    mul = 1
    for k in range(NL):
        radix.append(mul)
        mul *= n_states[k]
    NKEYS = mul
    RECS = [None] * NKEYS
    REF = {}          # (S) -> Y -> [N, {e: [W0, WT, WTu]}]
    stat = {"clases": 0, "asignaciones": 0, "por_componentes": {}, "dos_o_mas_componentes_con_acciones": 0,
            "intercalacion_estricta": 0, "cancelacion_tras_acciones": 0, "componente_con_dos_o_mas_acciones": 0,
            "fuente_con_respuesta_X_nula_en_soporte": 0, "max_acciones": 0, "max_componentes": 0}
    for states in product(*[range(n) for n in n_states]):
        key = sum(st * radix[k] for k, st in enumerate(states))
        src = E = s = 0
        mult = 1
        wts = {e: F(1) for e in RATES}
        for k, st in enumerate(states):
            mult *= class_count[k][st]
            if st:
                src |= 1 << k
                E ^= ERESP[k][st]
                s ^= SRESP[k][st]
                for e in RATES:
                    wts[e] *= class_w[e][k][st]
        tr = TR2[s]
        Y = src | (TOUCH2[s] << NL)
        R = E ^ COR2[s]
        parts = connected_parts(Y, ADJH) if Y else []
        hall = xr = 0
        assigned = 0
        local = []
        for P in parts:
            sp, apm = P & LOCMASK, P >> NL
            need(sp != 0, "R2 componente sin fuentes", states=states, P=P)
            Ep = sP = 0
            for k in bits(sp):
                Ep ^= ERESP[k][states[k]]
                sP ^= SRESP[k][states[k]]
            exp = [j for j in tr if actions[j] & apm]
            need(all(actions[j] & ~apm == 0 for j in exp), "R2 accion que cruza componentes", states=states, P=P)
            need((sp | (TOUCH2[sP] << NL)) == P, "R2 restriccion: el soporte aislado no es la componente",
                 states=states, P=P, aislado=sp | (TOUCH2[sP] << NL))
            need(list(TR2[sP]) == exp, "R2 restriccion: la traza local no es la subsucesion global",
                 states=states, P=P)
            rp = Ep ^ COR2[sP]
            hp = lin(cols, rp)
            need(hp & hall == 0, "R2 huellas fisicas de dos componentes se cortan", states=states, P=P)
            hall |= hp
            xr ^= rp
            assigned += len(exp)
            local.append((P, sP, rp, hp, len(exp)))
        need(assigned == len(tr), "R2 una accion no pertenece a exactamente una componente", states=states)
        need(xr == R, "R2 el residuo no es el XOR de los residuos locales", states=states)
        need(hall == lin(cols, R), "R2 el sindrome fisico no es la union de los locales", states=states)
        ss = [x[1] for x in local]
        merged, labels = [], []
        while True:
            best = None
            for hidx, sl in enumerate(ss):
                j, nxt = T2[sl]
                if j >= 0:
                    kk = (score2(sl, j), -j)
                    if best is None or kk > best[0]:
                        best = (kk, hidx, j, nxt)
            if best is None:
                break
            merged.append(best[2])
            labels.append(best[1])
            ss[best[1]] = best[3]
            need(len(merged) <= 12, "R2 la recomposicion no termina", states=states)
        need(merged == list(tr), "R2 recomposicion canonica distinta de la traza global", states=states,
             merged=merged, traza=list(tr))
        zero = hall == 0
        mark = -1 if bin(U0 & R).count("1") & 1 else 1
        silent = zero and mark < 0
        closed, cur, strict = set(), None, False
        for lb in labels:
            if lb != cur:
                if lb in closed:
                    strict = True
                if cur is not None:
                    closed.add(cur)
                cur = lb
        n_active = sum(1 for x in local if x[4])
        bad_local = sum(1 for x in local if x[3] == 0 and bin(U0 & x[2]).count("1") & 1 and x[4])
        cancel = bad_local >= 2 and zero and mark > 0
        RECS[key] = (E, s, R, zero, silent, mark)
        stat["clases"] += 1
        stat["asignaciones"] += mult
        nc = str(len(parts))
        stat["por_componentes"][nc] = stat["por_componentes"].get(nc, 0) + mult
        if n_active >= 2:
            stat["dos_o_mas_componentes_con_acciones"] += mult
        if strict:
            stat["intercalacion_estricta"] += mult
        if cancel:
            stat["cancelacion_tras_acciones"] += mult
        if any(x[4] >= 2 for x in local):
            stat["componente_con_dos_o_mas_acciones"] += mult
        if any(st and ERESP[k][st] == 0 and SRESP[k][st] == 0 for k, st in enumerate(states)):
            stat["fuente_con_respuesta_X_nula_en_soporte"] += mult
        stat["max_acciones"] = max(stat["max_acciones"], len(tr))
        stat["max_componentes"] = max(stat["max_componentes"], len(parts))
        if src:
            slot = REF.setdefault(src, {}).setdefault(Y, [0, {e: [F(0), F(0), F(0)] for e in RATES}, len(parts)])
            slot[0] += mult
            for e in RATES:
                slot[1][e][0] += wts[e]
                if zero:
                    slot[1][e][1] += wts[e]
                    slot[1][e][2] += mark * wts[e]
    check("D1 R2: restriccion, trazas locales, XOR de residuos, filtro fisico local y recomposicion canonica "
          "en todas las clases estructurales", True, {"clases": stat["clases"]})
    check("D2 las clases estructurales suman todas las asignaciones", stat["asignaciones"] == total,
          {"asignaciones": stat["asignaciones"]})
    tick("clases_estructurales", t)

    # ------------------------------------------------ 5. R1: todas las asignaciones
    t = time.perf_counter()
    keymul = [[state_of[k][c] * radix[k] for c in range(alph[k])] for k in range(NL)]
    c_tot, c_zero, c_sil, c_chT, c_ch = ([0] * NSIG for _ in range(5))
    memoR = {}
    done = 0
    rng = range(NL)
    for cats in product(*[range(x) for x in alph]):
        E1, s1 = literal(cats)
        key = sg = 0
        for k in rng:
            c = cats[k]
            key += keymul[k][c]
            sg += sig_of[k][c]
        rec = RECS[key]
        if E1 != rec[0] or s1 != rec[1]:
            need(False, "R1/R2 respuesta conjunta distinta de la formula", cats=cats, literal=[E1, s1],
                 formula=[rec[0], rec[1]])
        R = E1 ^ COR1[s1]
        fl = memoR.get(R)
        if fl is None:
            z0 = lin(cols, R) == 0
            fl = (z0, z0 and not in_span(R, bz), -1 if bin(U0 & R).count("1") & 1 else 1)
            memoR[R] = fl
        if R != rec[2] or fl[0] != rec[3] or fl[1] != rec[4]:
            need(False, "R1/R2 residuo, filtro fisico o fallo silencioso distintos", cats=cats,
                 global_=[R, fl[0], fl[1]], componentes=[rec[2], rec[3], rec[4]])
        c_tot[sg] += 1
        c_ch[sg] += fl[2]
        if fl[0]:
            c_zero[sg] += 1
            c_chT[sg] += fl[2]
            if fl[1]:
                c_sil[sg] += 1
        done += 1
        if done & 0xFFFF == 0:
            guard()
    COUNT["asignaciones_propagadas_literalmente"] = done
    COUNT["puertas_CNOT_aplicadas"] = done * 42
    COUNT["residuos_distintos"] = len(memoR)
    check("E1 R1 recorre todas las asignaciones", done == total, {"asignaciones": done})
    check("E2 respuesta conjunta literal = formula del apendice C en todas las asignaciones", True)
    check("E3 residuo, filtro fisico y fallo silencioso globales = los de la recomposicion en todas las "
          "asignaciones (clase por rango de H_Z frente a marca u_0)", True)
    tick("referencia_global", t)

    probs = {}
    for e in RATES:
        sp = [sig_prob(g, e) for g in range(NSIG)]
        tot = sum(c_tot[g] * sp[g] for g in range(NSIG))
        P0 = sum(c_zero[g] * sp[g] for g in range(NSIG))
        Ps = sum(c_sil[g] * sp[g] for g in range(NSIG))
        chT = sum(c_chT[g] * sp[g] for g in range(NSIG))
        ch = sum(c_ch[g] * sp[g] for g in range(NSIG))
        probs[e] = (tot, P0, Ps, chT, ch)
    check("E4 normalizacion: las probabilidades de todas las asignaciones suman 1 (dos tasas)",
          all(probs[e][0] == 1 for e in RATES))
    check("E5 P_0 - E[(-1)^{u.R}; H_X R=0] = 2 P_sil (dos tasas)",
          all(probs[e][1] - probs[e][3] == 2 * probs[e][2] for e in RATES))

    # ------------------------------------------------ 6. R3: oraculo de primera salida
    t = time.perf_counter()
    inverse = [[] for _ in range(M)]
    for j, ft in enumerate(foot):
        for c in bits(ft):
            inverse[c].append(j)
    widths = [bin(p).count("1") for p in foot]
    polymers = sorted(Y for S in REF for Y, slot in REF[S].items() if slot[2] == 1)
    variants = set()
    for Y in polymers:
        A = Y >> NL
        for v in bits(A):
            variants.add(Y & ~(1 << (NL + v)))
        nb = 0
        for x in bits(Y):
            nb |= ADJH[x]
        for x in bits(nb & ~Y & ~LOCMASK):
            variants.add(Y | (1 << x))
    variants = sorted(variants - set(polymers))
    need(len(polymers) <= MAX_ORACLE_SUPPORTS, "techo de soportes consultados", n=len(polymers))
    room = MAX_ORACLE_SUPPORTS - len(polymers)
    COUNT["oraculo_variantes_generadas"] = len(variants)
    COUNT["oraculo_variantes_omitidas_por_techo"] = max(0, len(variants) - room)
    queries = sorted(set(polymers) | set(variants[:room]))
    rid = {}
    memo_run = {}
    outcomes = {"accept": 0, "underfilled": 0, "exit": 0}
    oracle_states = 0
    ORW = {}

    def restricted_run(s, A, readers, rk):
        kk = (s, A, rk)
        got = memo_run.get(kk)
        if got is not None:
            return got
        tr, c, seen = [], 0, 0
        while True:
            best, num, den = -1, 0, 1
            for j in readers:
                d = widths[j]
                g = 2 * bin(s & foot[j]).count("1") - d
                if best < 0 or g * den > num * d:
                    best, num, den = j, g, d
            COUNT["oraculo_scores"] = COUNT.get("oraculo_scores", 0) + len(readers)
            if best < 0 or 2 * num < den:
                got = ("accept" if seen == A else "underfilled", tuple(tr), c)
                break
            a_ = actions[best]
            tr.append(best)
            if a_ & ~A:
                got = ("exit", tuple(tr), c)
                break
            c ^= a_
            seen |= a_
            s ^= foot[best]
        memo_run[kk] = got
        return got

    def expected_prefix(s, A):
        tr, c, seen = [], 0, 0
        for j in TR1[s]:
            tr.append(j)
            a_ = actions[j]
            if a_ & ~A:
                return ("exit", tuple(tr), c)
            c ^= a_
            seen |= a_
        return ("accept" if seen == A else "underfilled", tuple(tr), c)

    for Y in queries:
        S, A = Y & LOCMASK, Y >> NL
        PY = 0
        for k in bits(S):
            PY |= psi_loc[k]
        for v in bits(A):
            PY |= cols[v]
        readers = tuple(sorted({j for c in bits(PY) for j in inverse[c]}))
        rk = rid.setdefault(readers, len(rid))
        ks = list(bits(S))
        agg = [0, {e: [F(0), F(0), F(0)] for e in RATES}]
        for st in product(*[range(1, n_states[k]) for k in ks]):
            E = s = 0
            mult = 1
            for k, sk in zip(ks, st):
                E ^= ERESP[k][sk]
                s ^= SRESP[k][sk]
                mult *= class_count[k][sk]
            got = restricted_run(s, A, readers, rk)
            if got != expected_prefix(s, A):
                need(False, "R3 el oraculo no reproduce el prefijo de la ejecucion completa", Y=Y, estados=st,
                     oraculo=got, completa=expected_prefix(s, A))
            outcomes[got[0]] += 1
            oracle_states += 1
            if got[0] == "accept":
                r_ = E ^ got[2]
                agg[0] += mult
                zero = lin(cols, r_) == 0
                mk = -1 if bin(U0 & r_).count("1") & 1 else 1
                for e in RATES:
                    w = F(1)
                    for k, sk in zip(ks, st):
                        w *= class_w[e][k][sk]
                    agg[1][e][0] += w
                    if zero:
                        agg[1][e][1] += w
                        agg[1][e][2] += mk * w
        ref = REF.get(S, {}).get(Y)
        refN = ref[0] if ref else 0
        refW = ref[1] if ref else {e: [F(0), F(0), F(0)] for e in RATES}
        if agg[0] != refN or any(agg[1][e] != refW[e] for e in RATES):
            need(False, "R3 pesos del oraculo distintos de los de las ejecuciones completas", Y=Y,
                 oraculo_N=agg[0], completa_N=refN)
        ORW[Y] = agg
        guard()
    COUNT["oraculo_soportes"] = len(queries)
    COUNT["oraculo_consultas_fuente_soporte"] = oracle_states
    COUNT["oraculo_ejecuciones_restringidas_distintas"] = len(memo_run)
    COUNT["oraculo_desenlaces"] = outcomes
    check("F1 oraculo de primera salida: desenlace y traza iguales al prefijo de la ejecucion completa en todas "
          "las consultas", True, {"soportes": len(queries), "consultas": oracle_states, "desenlaces": dict(outcomes)})
    check("F2 pesos W_0, W_T y W_{T,u} del oraculo iguales a los de las ejecuciones completas en todos los "
          "soportes consultados (dos tasas)", True)
    check("F3 hay consultas con los tres desenlaces", all(outcomes[k] > 0 for k in outcomes) or not full_run,
          dict(outcomes))
    tick("oraculo", t)

    # ------------------------------------------------ 7. R4: gas de polimeros
    t = time.perf_counter()
    nbm = []
    for Y in polymers:
        m_ = Y
        for x in bits(Y):
            m_ |= ADJH[x]
        nbm.append(m_)
    pw = [(ORW[Y][0], [ORW[Y][1][e] for e in RATES]) for Y in polymers]
    npoly = len(polymers)
    acc = {"N": 0, "fam": 0, "maxsize": 0, "X": [[F(0), F(0), F(0)] for _ in RATES]}
    one = [[F(1), F(1), F(1)] for _ in RATES]

    def rec_fam(start, blocked, cnt, ws, size):
        acc["N"] += cnt
        acc["fam"] += 1
        if size > acc["maxsize"]:
            acc["maxsize"] = size
        for ri in range(len(RATES)):
            for q in range(3):
                acc["X"][ri][q] += ws[ri][q]
        for i in range(start, npoly):
            if polymers[i] & blocked == 0:
                n_i, w_i = pw[i]
                rec_fam(i + 1, blocked | nbm[i], cnt * n_i,
                        [[ws[ri][q] * w_i[ri][q] for q in range(3)] for ri in range(len(RATES))], size + 1)
        if acc["fam"] & 0x3FFF == 0:
            guard()

    rec_fam(0, 0, 1, one, 0)
    COUNT["polimeros"] = npoly
    COUNT["familias_compatibles"] = acc["fam"]
    COUNT["tamano_maximo_de_familia"] = acc["maxsize"]
    check("G1 identidad de recuento: familias compatibles con multiplicidad = numero de asignaciones",
          acc["N"] == total, {"familias": acc["fam"], "suma": acc["N"], "polimeros": npoly})
    check("G2 B_c * Xi_0 = 1 (dos tasas)", all(Bc[e] * acc["X"][ri][0] == 1 for ri, e in enumerate(RATES)))
    check("G3 P_0 = B_c * Xi_T frente a la suma global de R1 (dos tasas)",
          all(Bc[e] * acc["X"][ri][1] == probs[e][1] for ri, e in enumerate(RATES)))
    check("G4 E[(-1)^{u.R}; H_X R=0] = B_c * Xi_{T,u} frente a R1 (dos tasas)",
          all(Bc[e] * acc["X"][ri][2] == probs[e][3] for ri, e in enumerate(RATES)))
    check("G5 P_sil = B_c (Xi_T - Xi_{T,u}) / 2 frente a R1 (dos tasas)",
          all(Bc[e] * (acc["X"][ri][1] - acc["X"][ri][2]) / 2 == probs[e][2] for ri, e in enumerate(RATES)))
    tick("gas", t)

    # ------------------------------------------------ 8. existencia y predicciones
    if full_run:
        check("H1 existen asignaciones con dos o mas componentes con acciones",
              stat["dos_o_mas_componentes_con_acciones"] > 0, stat["dos_o_mas_componentes_con_acciones"])
        check("H2 existen asignaciones con intercalacion estricta", stat["intercalacion_estricta"] > 0,
              stat["intercalacion_estricta"])
        check("H3 existen asignaciones con cancelacion entre componentes que han ejecutado acciones",
              stat["cancelacion_tras_acciones"] > 0, stat["cancelacion_tras_acciones"])

        idx = {nm: k for k, nm in enumerate(names)}

        def case(**kw):
            cats = [0] * NL
            for nm, c in kw.items():
                cats[idx[nm]] = c
            cats = tuple(cats)
            E1, s1 = literal(cats)
            tr = [actions[j] for j in TR1[s1]]
            R = E1 ^ COR1[s1]
            src = sum(1 << k for k in range(NL) if cats[k])
            Y = src | (TOUCH1[s1] << NL)
            parts = connected_parts(Y, ADJH) if Y else []
            desc = sorted(sorted([names[k] for k in bits(P & LOCMASK)] + ["a%d" % v for v in bits(P >> NL)])
                          for P in parts)
            local_res = []
            for P in parts:
                Ep = sP = 0
                for k in bits(P & LOCMASK):
                    stt = state_of[k][cats[k]]
                    Ep ^= ERESP[k][stt]
                    sP ^= SRESP[k][stt]
                local_res.append(Ep ^ COR1[sP])
            return {"E": E1, "s": s1, "traza": tr, "R": R, "fisico": lin(cols, R), "en_S": in_span(R, bz),
                    "componentes": desc, "residuos_locales": sorted(local_res)}

        p1 = case(W16=1, W3=1)
        check("P1 X en esperas 16 y 3: sindrome {0,2,3}, traza 0,3,4, dos componentes, intercalacion, residuo "
              "65553 estabilizador",
              p1["s"] == 13 and p1["traza"] == [1, 8, 16] and p1["R"] == 65553 and p1["en_S"]
              and p1["componentes"] == [["W16", "a0", "a4"], ["W3", "a3"]], p1)
        p2 = case(W1=1, W2=1)
        check("P2 X en esperas 1 y 2: traza 0,3; una componente; residuo 15; fallo silencioso",
              p2["traza"] == [1, 8] and p2["R"] == 15 and p2["fisico"] == 0 and not p2["en_S"]
              and len(p2["componentes"]) == 1, p2)
        p3 = case(W1=1, W2=1, W9=1, W10=1)
        check("P3 X en esperas 1,2,9,10: traza 0,3,8,11; dos componentes con residuos 15 y 3840; suma 3855 en S_X",
              p3["traza"] == [1, 8, 256, 2048] and p3["R"] == 3855 and p3["fisico"] == 0 and p3["en_S"]
              and p3["residuos_locales"] == [15, 3840] and len(p3["componentes"]) == 2
              and lin(cols, 15) == 0 and lin(cols, 3840) == 0 and not in_span(15, bz) and not in_span(3840, bz), p3)
        p4 = case(W1=3)
        check("P4 Z en la espera 1: respuesta X nula, traza vacia, componente {W1}, exito",
              p4["E"] == 0 and p4["s"] == 0 and p4["traza"] == [] and p4["componentes"] == [["W1"]] and p4["en_S"], p4)
        p5a, p5b = case(R1=4), case(R1=1)
        check("P5 lectura del check 1 con (I,flip) y con (X,sin flip): sindrome {1}, sin accion, residuo 0, "
              "componente {R1}",
              all(x["s"] == 2 and x["traza"] == [] and x["R"] == 0 and x["componentes"] == [["R1"]]
                  for x in (p5a, p5b)), [p5a, p5b])
        p6a, p6b, p6c = case(C10_8=4), case(C10_8=5), case(C10_8=1)
        check("P6 CNOT ruidoso: (1,0) sindrome {7}, sin accion, fisico {7,8}; (1,1) volteo de 10, residuo 0; "
              "(0,1) volteo de 11, residuo 2048, fisico {8}",
              p6a["E"] == 1024 and p6a["s"] == 128 and p6a["traza"] == [] and p6a["fisico"] == 384
              and p6b["s"] == 384 and p6b["traza"] == [1024] and p6b["R"] == 0
              and p6c["E"] == 0 and p6c["s"] == 256 and p6c["traza"] == [2048] and p6c["R"] == 2048
              and p6c["fisico"] == 256, [p6a, p6b, p6c])

    # ------------------------------------------------ 9. salida
    def show(x):
        with localcontext() as ctx:
            ctx.prec = 30
            return {"exacto": str(x), "decimal": format(Decimal(x.numerator) / Decimal(x.denominator), ".12E")}

    result["instancia"] = {
        "codigo": "PATH25, semilla (3,6,12), 4 columnas", "n": N, "checks": M, "candidatos": len(actions),
        "beta": "1/2", "orden": "(|F|, mascara) ascendente", "L": L, "N_loc": n_loc, "vertices_host_primitivo": n_loc + N,
        "ubicaciones_ruidosas": [{"nombre": nm, "ubicacion": list(l), "categorias": al}
                                 for nm, l, al in zip(names, locs, alph)],
        "asignaciones": total, "clases_estructurales": NKEYS,
    }
    result["archivos_congelados_sha256"] = frozen_sha
    result["estadistica_con_multiplicidad"] = stat
    result["recuentos_por_tipo"] = {
        "tipo": "indice = esperas_no_base + (n_esperas+1)*(cnot_no_base + 2*(pauli_lectura + 2*flip))",
        "n_esperas": n_wait, "total": c_tot, "filtro_fisico_nulo": c_zero, "fallo_silencioso": c_sil,
        "caracter_filtrado": c_chT, "caracter_sin_filtrar": c_ch,
    }
    result["probabilidades"] = {}
    for ri, e in enumerate(RATES):
        tot, P0, Ps, chT, ch = probs[e]
        result["probabilidades"][str(e)] = {
            "B_c": show(Bc[e]), "P_0": show(P0), "1-P_0": show(1 - P0), "P_sil": show(Ps),
            "Q_sil": show(Ps / P0), "lambda_T": show(chT / P0), "lambda_sin_filtro": show(ch),
            "Xi_0": show(acc["X"][ri][0]), "Xi_T": show(acc["X"][ri][1]), "Xi_Tu": show(acc["X"][ri][2]),
        }
    result["alcance"] = ("una instancia, ocho ubicaciones ruidosas y el resto ideal; contrato generico no regular; "
                         "no es un estimador ni una medida de ventaja")


if __name__ == "__main__":
    main()
