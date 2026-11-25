#!/usr/bin/env python3
"""H3 - Busqueda implicita exacta integrada en el oraculo de primera salida.

Uso:  python -I -B verificar_h3.py --output <nuevo.json>
Solo biblioteca estandar y aritmetica entera exacta. La salida estandar y el
JSON son deterministas; los tiempos de reloj van a la salida de errores.

Partes (fijadas en PLAN_UNIDAD.txt antes de calcular):
 A. Rectangulos pequenos: ganador implicito contra fuerza bruta.
 B. Planos afines AG(2,q), q=2,3,5: lista explicita contra busqueda implicita
    (recuentos, trazas del decoder, oraculo de primera salida, host).
 C. F17 real: recuentos antes/despues, rectangulos densos y seis consultas
    fuente-soporte de respuesta conocida, sin materializar candidatos.
No se ejecuta ningun censo, muestreo sobre F17, SAT ni CNF.
"""
import argparse
import json
import os
import sys
import time
from math import comb

M64 = (1 << 64) - 1
CHECKS = []
EVALS = [0]          # evaluaciones de representante de la busqueda implicita
BN, BD = 1, 100      # beta = 1/100


def check(name, cond):
    CHECKS.append(name)
    if not cond:
        raise SystemExit("FALLA: " + name)


def note(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


class Rng:
    """xorshift64* determinista; no depende de la version de Python."""

    def __init__(self, seed):
        self.s = ((seed * 0x9E3779B97F4A7C15) + 1) & M64 or 1

    def next(self):
        x = self.s
        x ^= x >> 12
        x ^= (x << 25) & M64
        x ^= x >> 27
        self.s = x
        return (x * 0x2545F4914F6CDD1D) & M64

    def below(self, n):
        return (self.next() >> 11) % n

    def sample(self, pop, k):
        pop = list(pop)
        out = []
        for _ in range(k):
            out.append(pop.pop(self.below(len(pop))))
        return out


def bits(x):
    while x:
        low = x & -x
        yield low.bit_length() - 1
        x ^= low


# =====================================================================
# Busqueda implicita en un rectangulo (resultado heredado de 2R9, mas la
# compresion por lineas sin sindrome)
# =====================================================================
_TABLES = {}


def group_table(nS, nE, Af, kmode):
    """tab[e][s] = (huella, kappa) del grupo con e lineas enumeradas y s
    ordenadas, o None si el grupo no es candidato (vacio, huella vacia o
    filtro 2|phi| >= Af|F|). kmode 0: kappa=-|phi|; kmode 1: kappa=|F|."""
    key = (nS, nE, Af, kmode)
    t = _TABLES.get(key)
    if t is None:
        t = []
        for e in range(nE + 1):
            row = []
            for s in range(nS + 1):
                v = s * nE + e * nS - 2 * s * e
                if e + s == 0 or v <= 0 or 2 * v < Af * (e + s):
                    row.append(None)
                else:
                    row.append((v, -v if kmode == 0 else e + s))
            t.append(row)
        _TABLES[key] = t
    return t


def _members(sid, eid, inact, sub, y0, rows):
    out = [sid[i] for i in rows]
    j = 0
    while sub:
        if sub & 1:
            out.append(eid[j])
        sub >>= 1
        j += 1
    out.extend(inact[:y0])
    out.sort(reverse=True)
    return tuple(out)


def core(nS, nE, sid, eid, MS, Af, kmode, compress):
    """Primer maximizador del score en un rectangulo.
    Lineas ordenadas: nS, indices fisicos sid. Lineas enumeradas: nE, indices
    eid. MS[i] = bits de la linea ordenada i sobre las enumeradas.
    Devuelve (h, |phi|, kappa, miembros en orden descendente)."""
    tab = group_table(nS, nE, Af, kmode)
    r = [m.bit_count() for m in MS]
    if compress:
        act = 0
        for m in MS:
            act |= m
    else:
        act = (1 << nE) - 1
    inact = sorted(eid[j] for j in range(nE) if not (act >> j) & 1)
    ni = len(inact)
    bh = -1
    bv = bk = 0
    bdesc = None
    bmem = None
    rng = range(nS)
    srange = range(nS + 1)
    ev = 0
    sub = act
    while True:
        ya = sub.bit_count()
        d = []
        h0 = 0
        for i in rng:
            t = (MS[i] & sub).bit_count()
            h0 += t
            d.append((2 * t - r[i], sid[i], i))
        d.sort()
        hs = [h0]
        h = h0
        for item in d:
            h -= item[0]
            hs.append(h)
        for y0 in range(ni + 1):
            row = tab[ya + y0]
            for s in srange:
                g = row[s]
                if g is None:
                    continue
                ev += 1
                v, k = g
                h = hs[s]
                mem = None
                if bh >= 0:
                    c = h * bv - bh * v
                    if c < 0:
                        continue
                    if c == 0:
                        if k > bk:
                            continue
                        if k == bk:
                            mem = _members(sid, eid, inact, sub, y0, [x[2] for x in d[:s]])
                            if bmem is None:
                                bmem = _members(sid, eid, inact, *bdesc)
                            if not mem < bmem:
                                continue
                bh, bv, bk = h, v, k
                bdesc = (sub, y0, [x[2] for x in d[:s]])
                bmem = mem
        if sub == 0:
            break
        sub = (sub - 1) & act
    EVALS[0] += ev
    if bh < 0:
        return None
    if bmem is None:
        bmem = _members(sid, eid, inact, *bdesc)
    return (bh, bv, bk, bmem)


def rect_best(rid, cid, M, Af, kmode, orient, compress):
    """rid/cid: indices fisicos de filas (sitios VV) y columnas (sitios CC);
    M[i]: bits de la fila i sobre las columnas. orient: 'col' enumera
    columnas (2R9), 'row' enumera filas, 'auto' el lado mas barato."""
    R, C = len(rid), len(cid)
    if orient == 'auto':
        actc = 0
        kr = 0
        for m in M:
            actc |= m
            if m:
                kr += 1
        kc = actc.bit_count()
        costc = (1 << kc) * (C - kc + 1) * (R + 1)
        costr = (1 << kr) * (R - kr + 1) * (C + 1)
        orient = 'col' if costc <= costr else 'row'
    if orient == 'col':
        return core(R, C, rid, cid, M, Af, kmode, compress)
    MT = [0] * C
    for i, m in enumerate(M):
        for j in bits(m):
            MT[j] |= 1 << i
    return core(C, R, cid, rid, MT, Af, kmode, compress)


def better(a, b):
    """a precede a b: score mayor, luego kappa menor, luego mascara menor."""
    c = a[0] * b[1] - b[0] * a[1]
    if c:
        return c > 0
    if a[2] != b[2]:
        return a[2] < b[2]
    return a[3] < b[3]


def rect_brute(rid, cid, M, Af, kmode):
    """Fuerza bruta sobre todos los subconjuntos del generador."""
    R, C = len(rid), len(cid)
    fullc = (1 << C) - 1
    best = None
    for X in range(1 << R):
        x = X.bit_count()
        for Y in range(1 << C):
            y = Y.bit_count()
            if x + y == 0:
                continue
            v = x * C + y * R - 2 * x * y
            if v <= 0 or 2 * v < Af * (x + y):
                continue
            h = 0
            nY = fullc & ~Y
            for i in range(R):
                h += (M[i] & (nY if (X >> i) & 1 else Y)).bit_count()
            k = -v if kmode == 0 else x + y
            if best is not None:
                c = h * best[1] - best[0] * v
                if c < 0 or (c == 0 and k > best[2]):
                    continue
            mem = tuple(sorted([rid[i] for i in range(R) if (X >> i) & 1]
                               + [cid[j] for j in range(C) if (Y >> j) & 1], reverse=True))
            cand = (h, v, k, mem)
            if best is None or better(cand, best):
                best = cand
    return best


# =====================================================================
# Parte A
# =====================================================================
VARIANTS = (('col', False), ('row', False), ('col', True), ('row', True), ('auto', True))


def part_a():
    t0 = time.time()
    rng = Rng(20261111)
    out = {"exhaustive": [], "sampled": []}
    total = 0

    def layouts(R, C):
        nat = (list(range(R)), list(range(R, R + C)))
        perm = rng.sample(range(3 * (R + C)), R + C)
        return (nat, (perm[:R], perm[R:]))

    def run(R, C, mats):
        nonlocal total
        cnt = 0
        for lay in layouts(R, C):
            rid, cid = lay
            for kmode in (0, 1):
                for Af in (C, 0):
                    for M in mats:
                        ref = rect_brute(rid, cid, M, Af, kmode)
                        for orient, comp in VARIANTS:
                            got = rect_best(rid, cid, M, Af, kmode, orient, comp)
                            if got != ref:
                                raise SystemExit("FALLA A: %r %r %r %r %r %r -> %r != %r" % (
                                    R, C, M, Af, kmode, (orient, comp), got, ref))
                            cnt += 1
        total += cnt
        return cnt

    for R, C in ((1, 1), (1, 2), (2, 1), (2, 2), (1, 3), (3, 1), (2, 3), (3, 2), (3, 3),
                 (2, 4), (4, 2), (3, 4), (4, 3)):
        mats = []
        for code in range(1 << (R * C)):
            mats.append([(code >> (C * i)) & ((1 << C) - 1) for i in range(R)])
        cnt = run(R, C, mats)
        out["exhaustive"].append({"rows": R, "cols": C, "matrices": len(mats), "comparisons": cnt})
    for R, C in ((4, 4), (5, 4), (4, 5), (5, 5), (6, 5), (5, 6)):
        mats = []
        full = (1 << C) - 1
        for t in range(100):
            mode = t % 5
            if mode == 0:        # muy dispersa
                M = [0] * R
                for _ in range(1 + rng.below(3)):
                    M[rng.below(R)] |= 1 << rng.below(C)
            elif mode == 4:      # estructurada: filas o columnas completas
                M = [0] * R
                for _ in range(1 + rng.below(2)):
                    M[rng.below(R)] = full
                if rng.below(2):
                    j = rng.below(C)
                    M = [m ^ (1 << j) for m in M]
            else:
                num = (3, 10, 17)[mode - 1] if mode < 4 else 10
                M = [sum(1 << j for j in range(C) if rng.below(20) < num) for _ in range(R)]
            mats.append(M)
        mats.append([full] * R)
        mats.append([0] * R)
        cnt = run(R, C, mats)
        out["sampled"].append({"rows": R, "cols": C, "matrices": len(mats), "comparisons": cnt})
    check("A: todas las variantes implicitas coinciden con la fuerza bruta", total > 0)
    out["total_comparisons"] = total
    out["variants"] = ["columnas sin compresion (2R9)", "filas sin compresion", "columnas con compresion",
                       "filas con compresion", "lado automatico con compresion"]
    out["orders"] = ["(-|phi|, mascara ascendente)", "(|F|, mascara ascendente)"]
    out["filters"] = ["FGL 2|phi|>=A|F|", "sin filtro (solo huella no vacia)"]
    note("parte A: %.1f s" % (time.time() - t0))
    return out


# =====================================================================
# Planos afines y codigo producto
# =====================================================================
class Code:
    """HGP de la incidencia punto-recta de AG(2,q), q primo, con la
    vectorizacion del manuscrito: punto (x,y)->q*x+y; recta y=mx+c->q*m+c;
    recta x=c->q*q+c; VV(i,l)=a*i+l; CC(c,d)=a*a+b*c+d; check(i,d)=b*i+d."""

    def __init__(self, q):
        self.q = q
        self.b = q * q
        self.a = q * q + q
        self.A = q
        self.B = q + 1
        lp = []
        for m in range(q):
            for c in range(q):
                lp.append(sorted(q * x + (m * x + c) % q for x in range(q)))
        for c in range(q):
            lp.append([q * c + y for y in range(q)])
        self.line_pts = lp
        pl = [[] for _ in range(self.b)]
        for l, pts in enumerate(lp):
            for d in pts:
                pl[d].append(l)
        self.pt_lines = pl
        self.n = self.a * self.a + self.b * self.b
        self.m = self.a * self.b
        self.G = self.a * self.b
        self.classes = [list(range(q * m, q * m + q)) for m in range(q)] + [list(range(q * q, q * q + q))]

    def fp_site(self, v):
        a, b = self.a, self.b
        if v < a * a:
            i, l = divmod(v, a)
            base = b * i
            return [base + d for d in self.line_pts[l]]
        c, d = divmod(v - a * a, b)
        return [b * i + d for i in self.pt_lines[c]]

    def name(self, v):
        a, b = self.a, self.b
        if v < a * a:
            return "VV(%d,%d)" % divmod(v, a)
        return "CC(%d,%d)" % divmod(v - a * a, b)


def eligible(h, v):
    return (2 * h - v) * BD >= BN * v


def implicit_step(code, s, cache, kmode, stats):
    """Paso del decoder sobre el conjunto de checks s: mejor representante
    entre los padres cuyo rectangulo toca s. Sin lista de candidatos."""
    a, b = code.a, code.b
    lp, pl = code.line_pts, code.pt_lines
    parents = set()
    for ch in s:
        i, d = divmod(ch, b)
        for r in lp[i]:
            base = a * r
            for j in pl[d]:
                parents.add(base + j)
    best = None
    for g in parents:
        r, j = divmod(g, a)
        rows = pl[r]
        cols = lp[j]
        M = []
        for i in rows:
            base = b * i
            mrow = 0
            bit = 1
            for d in cols:
                if base + d in s:
                    mrow |= bit
                bit <<= 1
            M.append(mrow)
        pat = tuple(M)
        ent = cache.get(g)
        if ent is not None and ent[0] == pat:
            w = ent[1]
        else:
            rid = [a * i + j for i in rows]
            cid = [a * a + b * r + d for d in cols]
            w = rect_best(rid, cid, M, code.A, kmode, 'auto', True)
            cache[g] = (pat, w)
            stats["rectangles"] += 1
        if best is None or better(w, best):
            best = w
    stats["parents"] += len(parents)
    stats["steps"] += 1
    return best


def run_implicit(code, E, D, Yv=None, kmode=0, stats=None):
    """Decoder completo (Yv None) o consulta de primera salida sobre Yv."""
    if stats is None:
        stats = {"rectangles": 0, "parents": 0, "steps": 0}
    s = set()
    for v in E:
        s.symmetric_difference_update(code.fp_site(v))
    s.symmetric_difference_update(D)
    cache = {}
    trace = []
    touched = set(E)
    last = None
    while s:
        w = implicit_step(code, s, cache, kmode, stats)
        last = w
        if not eligible(w[0], w[1]):
            break
        F = w[3]
        if Yv is not None and not Yv.issuperset(F):
            return {"status": "exit", "trace": trace, "exit": (F, w[0], w[1]), "s": s, "touched": touched}
        for v in F:
            s.symmetric_difference_update(code.fp_site(v))
        touched.update(F)
        trace.append((F, w[0], w[1]))
    return {"status": "stop", "trace": trace, "s": s, "touched": touched, "last": last}


# ---------------------------------------------------------------- lista explicita
def fp_masks(code):
    return [sum(1 << c for c in code.fp_site(v)) for v in range(code.n)]


def parent_sites(code, r, j):
    a, b = code.a, code.b
    return [a * i + j for i in code.pt_lines[r]] + [a * a + b * r + d for d in code.line_pts[j]]


def materialize(code, fpm, kmode):
    """Familia FGL deduplicada en el orden del decoder. Solo q pequeno."""
    A = code.A
    seen = {}
    inc = 0
    s1 = 0
    for r in range(code.b):
        for j in range(code.a):
            sites = parent_sites(code, r, j)
            w = len(sites)
            fpa = [0] * (1 << w)
            for sub in range(1, 1 << w):
                low = sub & -sub
                fp = fpa[sub ^ low] ^ fpm[sites[low.bit_length() - 1]]
                fpa[sub] = fp
                if not fp:
                    continue
                v = fp.bit_count()
                if 2 * v < A * sub.bit_count():
                    continue
                inc += 1
                s1 += v
                mem = tuple(sorted((sites[t] for t in range(w) if (sub >> t) & 1), reverse=True))
                seen[mem] = fp
    if kmode == 0:
        items = sorted(seen.items(), key=lambda kv: (-kv[1].bit_count(), kv[0]))
    else:
        items = sorted(seen.items(), key=lambda kv: (len(kv[0]), kv[0]))
    cands = [(mem, fp, fp.bit_count()) for mem, fp in items]
    return cands, inc, s1


def explicit_step(cands, sint, restrict):
    """Primer maximizador del score en la lista ordenada; restrict limita a
    los candidatos que leen P_Y (Proposicion del oraculo local)."""
    best = -1
    bh = bv = 0
    idx = 0
    for mem, fp, v in cands:
        if restrict is None or fp & restrict:
            h = (fp & sint).bit_count()
            if best < 0 or h * bv > bh * v:
                best = idx
                bh = h
                bv = v
        idx += 1
    return best, bh, bv


def run_explicit(code, fpm, cands, E, D, Yv=None, restrict=None):
    sint = 0
    for v in E:
        sint ^= fpm[v]
    for c in D:
        sint ^= 1 << c
    trace = []
    touched = set(E)
    while True:
        idx, h, v = explicit_step(cands, sint, restrict)
        if idx < 0 or not eligible(h, v):
            break
        mem = cands[idx][0]
        if Yv is not None and not Yv.issuperset(mem):
            return {"status": "exit", "trace": trace, "exit": (mem, h, v), "sint": sint, "touched": touched}
        sint ^= cands[idx][1]
        touched.update(mem)
        trace.append((mem, h, v))
    return {"status": "stop", "trace": trace, "sint": sint, "touched": touched}


def explicit_scan_run(code, fpm, E, D, kmode):
    """Decoder de referencia sin almacenamiento: en cada paso recorre todos
    los subconjuntos validos de todos los generadores (q=5)."""
    A = code.A
    sint = 0
    for v in E:
        sint ^= fpm[v]
    for c in D:
        sint ^= 1 << c
    trace = []
    plist = [parent_sites(code, r, j) for r in range(code.b) for j in range(code.a)]
    while sint:
        best = None
        for sites in plist:
            w = len(sites)
            fpa = [0] * (1 << w)
            for sub in range(1, 1 << w):
                low = sub & -sub
                fp = fpa[sub ^ low] ^ fpm[sites[low.bit_length() - 1]]
                fpa[sub] = fp
                if not fp:
                    continue
                v = fp.bit_count()
                nf = sub.bit_count()
                if 2 * v < A * nf:
                    continue
                h = (fp & sint).bit_count()
                k = -v if kmode == 0 else nf
                if best is not None:
                    c = h * best[1] - best[0] * v
                    if c < 0 or (c == 0 and k > best[2]):
                        continue
                mem = tuple(sorted((sites[t] for t in range(w) if (sub >> t) & 1), reverse=True))
                cand = (h, v, k, mem, fp)
                if best is None or better(cand, best):
                    best = cand
        if best is None or not eligible(best[0], best[1]):
            break
        sint ^= best[4]
        trace.append((best[3], best[0], best[1]))
    return trace, sint


# ---------------------------------------------------------------- host
def host_closed(code):
    """K = union de cliques de T_g (forma cerrada del manuscrito)."""
    a, b, n = code.a, code.b, code.n
    lp, pl = code.line_pts, code.pt_lines
    V2 = []
    for j in range(a):
        s = {j}
        for d in lp[j]:
            s.update(pl[d])
        V2.append(sorted(s))
    C2 = []
    for r in range(b):
        s = {r}
        for i in pl[r]:
            s.update(lp[i])
        C2.append(sorted(s))
    adj = [0] * (n + code.m)
    for r in range(b):
        for j in range(a):
            T = 0
            for i in pl[r]:
                for l in V2[j]:
                    T |= 1 << (a * i + l)
                for d in lp[j]:
                    T |= 1 << (n + b * i + d)
            for c in C2[r]:
                for d in lp[j]:
                    T |= 1 << (a * a + b * c + d)
            for x in bits(T):
                adj[x] |= T
    for x in range(len(adj)):
        adj[x] &= ~(1 << x)
    return adj


def host_from_candidates(code, cands):
    """K por las reglas (K1) y (K2) a partir de la lista materializada."""
    n, m = code.n, code.m
    rd = [1 << (n + c) for c in range(m)]
    for v in range(n):
        for c in code.fp_site(v):
            rd[c] |= 1 << v
    adj = [0] * (n + m)
    for mem, fp, v in cands:
        Rm = 0
        for c in bits(fp):
            Rm |= rd[c]
        Fm = 0
        for y in mem:
            Fm |= 1 << y
        tot = Rm | Fm
        for x in bits(Rm):
            adj[x] |= tot
        for y in mem:
            adj[y] |= Rm
    for x in range(len(adj)):
        adj[x] &= ~(1 << x)
    return adj


# ---------------------------------------------------------------- formulas de recuento
def valid_pairs(A, B):
    out = []
    for x in range(B + 1):
        for y in range(A + 1):
            v = A * x + B * y - 2 * x * y
            if x + y > 0 and v > 0 and 2 * v >= A * (x + y):
                out.append((x, y, v))
    return out


def counts(A, B, a, b):
    vp = valid_pairs(A, B)
    K = sum(comb(B, x) * comb(A, y) for x, y, v in vp)
    S1 = sum(comb(B, x) * comb(A, y) * v for x, y, v in vp)
    G = a * b
    n = a * a + b * b
    ncand = G * K - a * a * (A - 1) - b * b * (B - 1)
    sinv = G * S1 - a * a * A * (A - 1) - b * b * B * (B - 1)
    # candidatos que NO leen una fila completa del rectangulo (fila i fija)
    nr = 0
    for x, y, v in vp:
        if y == A and x >= 1:
            nr += comb(B - 1, x - 1)
        if y == 0:
            nr += comb(B - 1, x)
    read_row = A * a * (K - nr) - a * (A - 1) - A * b * (B - 1)
    return {"valid_pairs_xy": len(vp), "K_valid_per_parent": K, "S1_per_parent": S1, "parents": G, "n": n,
            "dedup_candidates": ncand, "inverse_list_entries": sinv, "not_reading_full_row_per_parent": nr,
            "distinct_candidates_reading_full_row": read_row, "touching_parents_full_row": A * a}


# =====================================================================
# Parte B
# =====================================================================
def protected_mask(code, l0):
    """u = (1_{P_0}+1_{P_inf}) e_{l0}^T en el bloque VV; se comprueba u en ker H_Z."""
    a = code.a
    u = set(a * i + l0 for i in code.classes[0] + code.classes[-1])
    ok = True
    for r in range(code.b):
        for j in range(a):
            ov = sum(1 for v in parent_sites(code, r, j) if v in u)
            if ov % 2:
                ok = False
            if j != l0 and ov:
                ok = False
    check("mascara protegida q=%d: paridad par con todos los generadores (u en ker H_Z)" % code.q, ok)
    return u


def gen_inputs(code, rng, n_small, n_heavy):
    n, m = code.n, code.m
    inputs = []
    for _ in range(n_small):
        ne = 2 + rng.below(2)
        nd = rng.below(2)
        inputs.append((tuple(sorted(rng.sample(range(n), ne))), tuple(sorted(rng.sample(range(m), nd)))))
    for _ in range(n_heavy):
        ne = 4 + rng.below(7)
        nd = rng.below(4)
        inputs.append((tuple(sorted(rng.sample(range(n), ne))), tuple(sorted(rng.sample(range(m), nd)))))
    return inputs


def class_inputs(code, rows):
    a = code.a
    return [(tuple(a * i + l for l in cl), ()) for i in rows for cl in code.classes]


def compare_traces(code, fpm, cands, inputs, kmode, tag):
    nsteps = 0
    nstop_nonzero = 0
    for E, D in inputs:
        ri = run_implicit(code, E, D, None, kmode)
        re = run_explicit(code, fpm, cands, E, D)
        if ri["trace"] != re["trace"]:
            raise SystemExit("FALLA %s: trazas distintas para E=%r D=%r" % (tag, E, D))
        sint = 0
        for c in ri["s"]:
            sint |= 1 << c
        if sint != re["sint"]:
            raise SystemExit("FALLA %s: sindrome final distinto para E=%r D=%r" % (tag, E, D))
        nsteps += len(ri["trace"])
        if sint:
            nstop_nonzero += 1
    return {"inputs": len(inputs), "accepted_actions": nsteps, "stops_with_nonzero_syndrome": nstop_nonzero}


def oracle_block(code, fpm, cands, adj, u, rng, n_pairs, n_rand, struct_rows):
    """Consultas fuente-soporte: implicito, explicito restringido a P_Y
    (Proposicion del oraculo local) y decoder global explicito."""
    n, m, a = code.n, code.m, code.a
    N = n + m
    supports = []
    for x in range(N):
        supports.append((x,))
    for _ in range(n_pairs):
        x = rng.below(N)
        nb = list(bits(adj[x]))
        supports.append(tuple(sorted((x, nb[rng.below(len(nb))]))))
    for _ in range(n_rand):
        size = 3 + rng.below(4)
        Y = {rng.below(N)}
        guard = 0
        while len(Y) < size and guard < 200:
            guard += 1
            x = sorted(Y)[rng.below(len(Y))]
            nb = list(bits(adj[x]))
            y = nb[rng.below(len(nb))]
            if y < n and sum(1 for z in Y if z < n) >= 6:
                continue
            Y.add(y)
        supports.append(tuple(sorted(Y)))
    for i in struct_rows:
        for cl in code.classes:
            supports.append(tuple(a * i + l for l in cl))
        for cl in code.classes[1:]:
            supports.append(tuple(sorted(a * i + l for l in cl + code.classes[0])))
    stat = {"accept": 0, "smaller": 0, "exit": 0}
    nq = 0
    wsum = {"W0": {}, "WT": {}, "Wu": {}}
    for Y in supports:
        Yv = set(x for x in Y if x < n)
        Yc = tuple(x - n for x in Y if x >= n)
        PY = 0
        for v in Yv:
            PY |= fpm[v]
        for c in Yc:
            PY |= 1 << c
        Dmask = 0
        for c in Yc:
            Dmask |= 1 << c
        yv = sorted(Yv)
        Wi = {"W0": {}, "WT": {}, "Wu": {}}
        We = {"W0": {}, "WT": {}, "Wu": {}}
        for sub in range(1 << len(yv)):
            E = tuple(yv[t] for t in range(len(yv)) if (sub >> t) & 1)
            ri = run_implicit(code, E, Yc, Yv)
            re = run_explicit(code, fpm, cands, E, Yc, Yv, PY)
            rg = run_explicit(code, fpm, cands, E, Yc)
            # desenlace derivado del decoder global
            gt = set(E)
            gstat = None
            gtrace = []
            for item in rg["trace"]:
                if not Yv.issuperset(item[0]):
                    gstat = ("exit", item)
                    break
                gt.update(item[0])
                gtrace.append(item)
            if gstat is None:
                gstat = ("accept" if gt == Yv else "smaller", None)
            outs = []
            for res in (ri, re):
                if res["status"] == "exit":
                    outs.append(("exit", res["exit"], res["trace"]))
                else:
                    outs.append(("accept" if res["touched"] == Yv else "smaller", None, res["trace"]))
            if outs[0] != outs[1] or outs[0] != (gstat[0], gstat[1], gtrace):
                raise SystemExit("FALLA B3: desenlaces distintos en Y=%r E=%r" % (Y, E))
            nq += 1
            stat[outs[0][0]] += 1
            if outs[0][0] == "accept":
                for res, W in ((ri, Wi), (re, We)):
                    R = set(E)
                    for item in res["trace"]:
                        R.symmetric_difference_update(item[0])
                    if "s" in res:
                        fin = 0
                        for c in res["s"]:
                            fin |= 1 << c
                    else:
                        fin = res["sint"]
                    e = len(E)
                    W["W0"][e] = W["W0"].get(e, 0) + 1
                    if fin == Dmask:                       # H_X R = 0
                        W["WT"][e] = W["WT"].get(e, 0) + 1
                        sign = -1 if len(R & u) % 2 else 1
                        W["Wu"][e] = W["Wu"].get(e, 0) + sign
        if Wi != We:
            raise SystemExit("FALLA B3: pesos distintos en Y=%r" % (Y,))
        for key in Wi:
            for e, c in Wi[key].items():
                wsum[key][e] = wsum[key].get(e, 0) + c
    return {"supports": len(supports), "queries": nq, "outcomes": stat,
            "sum_of_weight_coefficients_by_sources": {k: {str(e): c for e, c in sorted(v.items())}
                                                      for k, v in wsum.items()}}


def part_b():
    out = {}
    rng = Rng(3111)
    for q in (2, 3):
        t0 = time.time()
        code = Code(q)
        fpm = fp_masks(code)
        cands, inc, s1 = materialize(code, fpm, 0)
        cf = counts(code.A, code.B, code.a, code.b)
        check("B1 q=%d: incidencias validas = G*K_valid" % q, inc == cf["parents"] * cf["K_valid_per_parent"])
        check("B1 q=%d: suma de huellas = G*S1" % q, s1 == cf["parents"] * cf["S1_per_parent"])
        check("B1 q=%d: candidatos deduplicados = formula" % q, len(cands) == cf["dedup_candidates"])
        check("B1 q=%d: entradas de listas inversas = formula" % q,
              sum(c[2] for c in cands) == cf["inverse_list_entries"])
        i0 = 1
        row = sum(1 << (code.b * i0 + d) for d in range(code.b))
        nread = sum(1 for c in cands if c[1] & row)
        check("B1 q=%d: candidatos que leen una fila completa = formula" % q,
              nread == cf["distinct_candidates_reading_full_row"])
        res = {"n": code.n, "checks": code.m, "parents": code.G, "A": code.A, "B": code.B,
               "counts": cf, "explicit_list_built": len(cands)}
        # B4 host
        adj = host_closed(code)
        adj2 = host_from_candidates(code, cands)
        check("B4 q=%d: host por (K1)-(K2) = union de cliques de T_g" % q, adj == adj2)
        res["host"] = {"vertices": len(adj), "edges": sum(x.bit_count() for x in adj) // 2,
                       "max_degree": max(x.bit_count() for x in adj)}
        # B2 trazas
        n, m = code.n, code.m
        inputs = [((v,), ()) for v in range(n)] + [((), (c,)) for c in range(m)]
        if q == 2:
            inputs += [((v, w), ()) for v in range(n) for w in range(v + 1, n)]
            inputs += [((v,), (c,)) for v in range(n) for c in range(m)]
            inputs += [((), (c, e)) for c in range(m) for e in range(c + 1, m)]
            inputs += gen_inputs(code, rng, 300, 300)
            rows = (0, 3, 5)
        else:
            inputs += gen_inputs(code, rng, 1200, 300)
            rows = (0, 5, 11)
        inputs += class_inputs(code, rows)
        res["traces_order_footprint_mask"] = compare_traces(code, fpm, cands, inputs, 0, "B2 q=%d" % q)
        check("B2 q=%d: trazas completas iguales, orden (-|phi|,mascara)" % q, True)
        cands1, _, _ = materialize(code, fpm, 1)
        sub_inputs = inputs[::7][:400] + class_inputs(code, rows)
        res["traces_order_size_mask"] = compare_traces(code, fpm, cands1, sub_inputs, 1, "B2' q=%d" % q)
        check("B2 q=%d: trazas completas iguales, orden (|F|,mascara)" % q, True)
        del cands1
        # historia de clase: prediccion analoga a Q4/Q6
        a = code.a
        for i in rows:
            for mcls, cl in enumerate(code.classes):
                E = tuple(a * i + l for l in cl)
                ri = run_implicit(code, E, (), None, 0)
                exp = [((a * i + l,), code.A, code.A) for l in code.classes[0]]
                check("B5 q=%d: clase %d de la fila %d -> el decoder retira P_0 en orden" % (q, mcls, i),
                      ri["trace"] == exp and not ri["s"])
        # B3 oraculo
        u = protected_mask(code, q)
        if q == 2:
            ob = oracle_block(code, fpm, cands, adj, u, rng, 300, 200, rows)
        else:
            ob = oracle_block(code, fpm, cands, adj, u, rng, 200, 80, rows)
        check("B3 q=%d: desenlaces, trazas y pesos W_0, W_T, W_u iguales" % q, True)
        res["oracle"] = ob
        out["q%d" % q] = res
        note("parte B q=%d: %.1f s" % (q, time.time() - t0))
    # q=5: sin lista almacenada
    t0 = time.time()
    code = Code(5)
    fpm = fp_masks(code)
    cf = counts(code.A, code.B, code.a, code.b)
    inc = s1 = 0
    plist = [parent_sites(code, r, j) for r in range(code.b) for j in range(code.a)]
    for sites in plist:
        w = len(sites)
        fpa = [0] * (1 << w)
        for sub in range(1, 1 << w):
            low = sub & -sub
            fp = fpa[sub ^ low] ^ fpm[sites[low.bit_length() - 1]]
            fpa[sub] = fp
            if fp and 2 * fp.bit_count() >= code.A * sub.bit_count():
                inc += 1
                s1 += fp.bit_count()
    check("B1 q=5: incidencias validas = G*K_valid", inc == cf["parents"] * cf["K_valid_per_parent"])
    check("B1 q=5: suma de huellas = G*S1", s1 == cf["parents"] * cf["S1_per_parent"])
    psets = [frozenset(s) for s in plist]
    where = {}
    for g, s in enumerate(plist):
        for v in s:
            where.setdefault(v, []).append(g)
    maxshare = 0
    for g, s in enumerate(psets):
        seen = {}
        for v in s:
            for h in where[v]:
                if h != g:
                    seen[h] = seen.get(h, 0) + 1
        if seen:
            maxshare = max(maxshare, max(seen.values()))
    check("B1 q=5: dos generadores distintos comparten a lo sumo un qubit", maxshare == 1)
    rng5 = Rng(555)
    n, m, a = code.n, code.m, code.a
    inputs = [((v,), ()) for v in (0, 437, 899, 900, 1300, 1524)]
    inputs += gen_inputs(code, rng5, 4, 3)
    inputs += [(tuple(a * 7 + l for l in code.classes[2]), ()), (tuple(a * 7 + l for l in code.classes[0]), ()),
               ((), (123,))]
    nst = 0
    for E, D in inputs:
        ri = run_implicit(code, E, D, None, 0)
        tr, sint = explicit_scan_run(code, fpm, E, D, 0)
        fin = 0
        for c in ri["s"]:
            fin |= 1 << c
        if ri["trace"] != tr or fin != sint:
            raise SystemExit("FALLA B2 q=5: E=%r D=%r" % (E, D))
        nst += len(tr)
    check("B2 q=5: trazas iguales frente al recorrido completo de subconjuntos", True)
    ri = run_implicit(code, tuple(a * 7 + l for l in code.classes[2]), (), None, 0)
    check("B5 q=5: clase 2 de la fila 7 -> el decoder retira P_0 en orden",
          ri["trace"] == [((a * 7 + l,), 5, 5) for l in code.classes[0]] and not ri["s"])
    out["q5"] = {"n": n, "checks": m, "parents": code.G, "A": code.A, "B": code.B, "counts": cf,
                 "explicit_list_built": None,
                 "traces_order_footprint_mask": {"inputs": len(inputs), "accepted_actions": nst},
                 "max_shared_qubits_between_parents": maxshare}
    note("parte B q=5: %.1f s" % (time.time() - t0))
    return out


# =====================================================================
# Parte C: F17
# =====================================================================
def part_c():
    out = {}
    code = Code(17)
    a, b = code.a, code.b
    check("C: F17 tiene (a,b,A,B,n,m)=(306,289,17,18,177157,88434)",
          (a, b, code.A, code.B, code.n, code.m) == (306, 289, 17, 18, 177157, 88434))
    check("C: cada punto esta en 18 rectas y cada recta tiene 17 puntos",
          all(len(x) == 18 for x in code.pt_lines) and all(len(x) == 17 for x in code.line_pts))
    # dos puntos distintos, una recta comun; dos rectas distintas, a lo sumo un punto
    ok = True
    for l1 in range(a):
        s1 = set(code.line_pts[l1])
        for l2 in range(l1 + 1, a):
            if len(s1.intersection(code.line_pts[l2])) > 1:
                ok = False
    check("C: dos rectas distintas comparten a lo sumo un punto (duplicados = singletons)", ok)
    cf = counts(17, 18, a, b)
    worst = 19 * (1 << 17)
    cf["all_subsets_per_parent"] = (1 << 35) - 1
    cf["implicit_groups_worst_case_per_parent"] = worst
    cf["ratio_valid_over_worst_groups_floor"] = cf["K_valid_per_parent"] // worst
    cf["ratio_all_subsets_over_worst_groups_floor"] = (1 << 35) // worst
    cf["seed_incidences"] = sum(len(x) for x in code.line_pts)
    out["counts"] = cf
    check("C1: la incidencia de la semilla tiene 5202 entradas", cf["seed_incidences"] == 5202)

    # C2 rectangulos densos con los indices fisicos de un padre real
    t0 = time.time()
    r0, j0 = 100, 37
    rid = [a * i + j0 for i in code.pt_lines[r0]]
    cid = [a * a + b * r0 + d for d in code.line_pts[j0]]
    rng = Rng(1717)
    full = (1 << 17) - 1
    pats = [("aleatorio densidad 1/2", [rng.next() & full for _ in range(18)]),
            ("todo activo", [full] * 18),
            ("aleatorio densidad 1/5",
             [sum(1 << j for j in range(17) if rng.below(5) == 0) for _ in range(18)]),
            ("nueve filas iguales", [0x155AA if i % 2 == 0 else 0 for i in range(18)])]
    dense = []
    for name, M in pats:
        EVALS[0] = 0
        wc = rect_best(rid, cid, M, 17, 0, 'col', True)
        ec = EVALS[0]
        EVALS[0] = 0
        wr = rect_best(rid, cid, M, 17, 0, 'row', True)
        er = EVALS[0]
        check("C2: rectangulo 18x17 '%s': columnas = filas" % name, wc == wr)
        dense.append({"pattern": name, "syndrome_bits": sum(m.bit_count() for m in M),
                      "winner_h": wc[0], "winner_footprint": wc[1], "winner_size": len(wc[3]),
                      "group_evaluations_by_columns": ec, "group_evaluations_by_rows": er})
    out["dense_rectangles"] = dense
    note("parte C2: %.1f s" % (time.time() - t0))

    # C3 consultas fuente-soporte
    t0 = time.time()
    P0 = code.classes[0]
    P3 = code.classes[3]
    row = 5
    VV = lambda i, l: a * i + l
    CC = lambda c, d: a * a + b * c + d
    q = {}

    def run_q(tag, E, D, Yv):
        st = {"rectangles": 0, "parents": 0, "steps": 0}
        EVALS[0] = 0
        res = run_implicit(code, E, D, Yv, 0, st)
        st["group_evaluations"] = EVALS[0]
        if res["status"] == "exit":
            outcome = "exit"
        else:
            outcome = "accept" if res["touched"] == Yv else "smaller"
        q[tag] = {"sources": len(E) + len(D), "support_size": len(Yv) + len(D), "outcome": outcome,
                  "accepted_actions": len(res["trace"]), "work": st}
        return res, outcome

    v1 = VV(5, 100)
    res, oc = run_q("Q1", (v1,), (), {v1})
    check("Q1: VV(5,100) se corrige con su singleton (score 1, huella 17) y se acepta",
          oc == "accept" and res["trace"] == [((v1,), 17, 17)] and not res["s"])
    v2 = CC(7, 200)
    res, oc = run_q("Q2", (v2,), (), {v2})
    check("Q2: CC(7,200) se corrige con su singleton (huella 18) y se acepta",
          oc == "accept" and res["trace"] == [((v2,), 18, 18)] and not res["s"])
    c3 = b * 5 + 123
    res, oc = run_q("Q3", (), (c3,), set())
    check("Q3: un check aislado no tiene candidato elegible; parada inmediata y aceptacion",
          oc == "accept" and res["trace"] == [] and res["s"] == {c3})
    q["Q3"]["best_rejected"] = {"h": res["last"][0], "footprint": res["last"][1],
                                "action": [code.name(v) for v in res["last"][3]]}
    check("Q3: el mejor candidato tiene h=1, huella 17 y no alcanza beta",
          res["last"][:2] == (1, 17) and not eligible(1, 17))
    E4 = tuple(VV(row, l) for l in P3)
    Y4 = set(E4) | set(VV(row, l) for l in P0)
    res, oc = run_q("Q4", E4, (), Y4)
    exp = [((VV(row, l),), 17, 17) for l in P0]
    check("Q4: 17 pasos VV(5,0..16) en orden, score 1, huella 17, sindrome final 0, aceptada",
          oc == "accept" and res["trace"] == exp and not res["s"])
    R4 = set(E4)
    for item in res["trace"]:
        R4.symmetric_difference_update(item[0])
    check("Q4: el residuo son las dos clases completas (34 sitios)", R4 == Y4)
    u = protected_mask(code, 51)
    check("Q4: la mascara protegida tiene peso 34 y paridad impar con el residuo: fallo silencioso",
          len(u) == 34 and len(u & R4) % 2 == 1)
    q["Q4"]["residual_weight"] = len(R4)
    q["Q4"]["protected_mask_weight"] = len(u)
    q["Q4"]["explicit_candidates_reading_initial_syndrome"] = cf["distinct_candidates_reading_full_row"]
    Y5 = set(E4)
    res, oc = run_q("Q5", E4, (), Y5)
    check("Q5: con Y=E la primera accion VV(5,0) sale de Y: rechazo por primera salida",
          oc == "exit" and res["trace"] == [] and res["exit"] == ((VV(row, 0),), 17, 17))
    E6 = tuple(VV(row, l) for l in P0)
    res, oc = run_q("Q6", E6, (), set(E6))
    check("Q6: la clase P_0 se corrige en 17 pasos y se acepta con residuo 0",
          oc == "accept" and res["trace"] == exp and not res["s"])
    out["queries"] = q
    note("parte C3: %.1f s" % (time.time() - t0))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    if os.path.exists(args.output):
        raise SystemExit("El archivo de salida ya existe.")
    t0 = time.time()
    res = {"unit": "H3_BUSQUEDA_IMPLICITA_ORACULO_20261111", "beta": "1/100",
           "order": "(-|phi_F|, mascara fisica ascendente)"}
    res["A_rectangles"] = part_a()
    res["B_affine_planes"] = part_b()
    res["C_f17"] = part_c()
    res["checks"] = len(CHECKS)
    res["status"] = "PASS"
    with open(args.output, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(res, ensure_ascii=False, indent=1) + "\n")
    note("total: %.1f s" % (time.time() - t0))
    c = res["C_f17"]["counts"]
    print("F17 candidatos deduplicados (lista explicita):", c["dedup_candidates"])
    print("F17 entradas de listas inversas:", c["inverse_list_entries"])
    print("F17 subconjuntos validos por padre:", c["K_valid_per_parent"])
    print("F17 grupos implicitos por padre, peor caso:", c["implicit_groups_worst_case_per_parent"])
    for k in sorted(res["C_f17"]["queries"]):
        v = res["C_f17"]["queries"][k]
        print(k, v["outcome"], "acciones", v["accepted_actions"], "padres", v["work"]["parents"],
              "grupos", v["work"]["group_evaluations"])
    print("PASS controles=%d" % len(CHECKS))


if __name__ == "__main__":
    main()
