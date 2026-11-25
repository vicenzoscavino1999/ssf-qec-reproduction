"""Construccion HGP rectangular: primer factor h1 (b1 x a1), segundo h2 (b2 x a2).
Con h1 = h2 reproduce build() de control_exacto.py; la cola es copia literal."""
import hashlib
import importlib.util

BASE = "D:/paper ia cuantica/qec_revision_matematica/05_COMPROBACIONES/"
CE_PATH = BASE + "D1_CONTROL_EXACTO_Y_METODOLOGIA_20261113/control_exacto.py"
SHA_CE = "c003050183a14364750a33784c7c815c0c07b5814f552a8b52a14d142db3306b"


def load_ce():
    with open(CE_PATH, "rb") as fh:
        assert hashlib.sha256(fh.read()).hexdigest() == SHA_CE
    spec = importlib.util.spec_from_file_location("control_exacto", CE_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_rect(ce, name, rows1, a1, rows2, a2):
    bits, lin, span = ce.bits, ce.lin, ce.span
    b1, b2 = len(rows1), len(rows2)
    h1 = [[(rows1[r] >> c) & 1 for c in range(a1)] for r in range(b1)]
    h2 = [[(rows2[r] >> c) & 1 for c in range(a2)] for r in range(b2)]
    N, M = a1 * a2 + b1 * b2, a1 * b2
    NS = N + M
    hx, hz = [0] * (a1 * b2), [0] * (b1 * a2)
    for i in range(a1):
        for j in range(b2):
            row = 0
            for k in range(a2):
                if h2[j][k]:
                    row |= 1 << (i * a2 + k)
            for l in range(b1):
                if h1[l][i]:
                    row |= 1 << (a1 * a2 + l * b2 + j)
            hx[i * b2 + j] = row
    for i in range(b1):
        for j in range(a2):
            row = 0
            for k in range(a1):
                if h1[i][k]:
                    row |= 1 << (k * a2 + j)
            for l in range(b2):
                if h2[l][j]:
                    row |= 1 << (a1 * a2 + i * b2 + l)
            hz[i * a2 + j] = row
    # ---- cola: copia literal de build() de control_exacto.py
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
    # ---- coordenadas: las de build() con (a, b) de cada factor
    coord = [None] * NS
    for i in range(a1):
        for k in range(a2):
            coord[i * a2 + k] = (2 * i, 2 * k)
    for l in range(b1):
        for j in range(b2):
            coord[a1 * a2 + l * b2 + j] = (2 * l + 1, 2 * j + 1)
    for i in range(a1):
        for j in range(b2):
            coord[N + i * b2 + j] = (2 * i, 2 * j + 1)
    return {"name": name, "rows": rows1, "a": a1, "b": b1, "rows2": rows2, "a2": a2, "b2": b2,
            "N": N, "M": M, "NS": NS, "hx": hx, "hz": hz,
            "cols": cols, "bx": bx, "bz": bz, "klog": klog, "ACT": ACT, "FOOT": FOOT, "WID": WID, "INV": INV,
            "psi": psi, "adjK": adjK, "coord": coord}


def strip(ce, m, n):
    r1, a1 = ce.path_seed(m)
    r2, a2 = ce.path_seed(n)
    return build_rect(ce, "TIRA_%dx%d" % (m, n), r1, a1, r2, a2)
