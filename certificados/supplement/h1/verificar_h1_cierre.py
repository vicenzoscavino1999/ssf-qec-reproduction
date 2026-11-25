#!/usr/bin/env python3
"""H1 - cierre: test en el cociente logico y menorante mejorado.

Partes:
  A. Base explicita de 289 consultas logicas en F17 (emparejamiento identidad).
  B. Control exhaustivo en HGP([1 1 1]): n=10, k=4, d=7, todas las historias.
  C. Reverificacion por decoder de fila de la familia heredada (15,4) ASC.
  D. Antes/despues en las esquinas con acceso certificado (racionales exactos).
Solo biblioteca estandar. No ejecuta trayectorias FGL completas de F17.

Uso:  python -I -B verificar_h1_cierre.py --output <nuevo.json>
"""
import argparse
import json
import math
import os
import random
import sys
import time
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
from itertools import combinations

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

CHECKS = {}


def check(name, cond):
    CHECKS[name] = bool(cond)
    if not cond:
        raise SystemExit("FALLO DE CONTROL: " + name)


def sci(x, mode, digits=13):
    x = Fraction(x)
    if x == 0:
        return "0"
    ctx = Context(prec=digits, rounding=ROUND_FLOOR if mode == "lo" else ROUND_CEILING,
                  Emax=999999999, Emin=-999999999)
    return format(ctx.divide(Decimal(x.numerator), Decimal(x.denominator)), "E")


def both(x):
    return {"lower_reading": sci(x, "lo"), "upper_reading": sci(x, "hi")}


def log10_fr(x):
    x = Fraction(x)
    return (math.log(x.numerator) - math.log(x.denominator)) / math.log(10)


def ceil_log2(x):
    """Menor entero k con 2^k >= x (x racional positivo)."""
    x = Fraction(x)
    k = x.numerator.bit_length() - x.denominator.bit_length() - 2
    while Fraction(2) ** k < x:
        k += 1
    return k


def par(v):
    return v.bit_count() & 1


# =========================================================== presupuesto
EPS = Fraction(1, 10)
E_BUD = EPS / 8
DELTA_SB = E_BUD / 4


def budget(ell):
    eta = E_BUD * ell / 4
    return {"eta": eta, "zeta_num": eta / 4, "zeta_den": min(E_BUD, eta / 4)}


def aghp_field_bits(dim):
    """Construccion por potencias (Alon-Goldreich-Hastad-Peralta): sesgo <= (dim-1)/2^m."""
    return max(1, ceil_log2(Fraction(dim - 1) / DELTA_SB))


# ================================================================ parte A
def build_lines():
    masks = []
    for m in range(17):
        for c in range(17):
            mk = 0
            for x in range(17):
                mk |= 1 << (17 * x + (m * x + c) % 17)
            masks.append(mk)
    for c in range(17):
        mk = 0
        for y in range(17):
            mk |= 1 << (17 * c + y)
        masks.append(mk)
    return masks


LM = build_lines()
CLASSES = [list(range(17 * m, 17 * m + 17)) for m in range(17)] + [list(range(289, 306))]


def gf2_rank(rows):
    piv = {}
    for r in rows:
        while r:
            p = r.bit_length() - 1
            if p in piv:
                r ^= piv[p]
            else:
                piv[p] = r
                break
    return len(piv)


def part_a():
    # H_0: filas = puntos (289), columnas = rectas (306).
    point_rows = []
    for pt in range(289):
        row = 0
        for j in range(306):
            if LM[j] >> pt & 1:
                row |= 1 << j
        point_rows.append(row)
    rk = gf2_rank(point_rows)
    check("F17_seed_full_row_rank_289", rk == 289)

    def h0(cmask):          # sindrome de un conjunto de rectas
        s = 0
        j = 0
        while cmask:
            if cmask & 1:
                s ^= LM[j]
            cmask >>= 1
            j += 1
        return s

    inf = 0
    for j in CLASSES[17]:
        inf |= 1 << j
    cvec = []
    for m in range(17):
        c = inf
        for j in CLASSES[m]:
            c |= 1 << j
        cvec.append(c)
    check("17_kernel_vectors_of_weight_34", all(c.bit_count() == 34 and h0(c) == 0 for c in cvec))
    jpt = [17 * m for m in range(17)]            # recta y = m x, una por clase finita
    ident = all(((cvec[a] >> jpt[b]) & 1) == (1 if a == b else 0) for a in range(17) for b in range(17))
    check("evaluation_matrix_is_identity_17x17", ident)
    check("kernel_dimension_is_17", 306 - rk == 17 and gf2_rank(cvec) == 17)
    # Consultas Z: u_(a,b) = (c_a e_{j_b}^T, 0), columnas de M en ker H_0  => u en ker H_Z.
    # Logicos X:  L_(g,d) = (e_{j_g} c_d^T, 0), filas de U en ker H_0      => L en ker H_X.
    # Emparejamiento u_(a,b) . L_(g,d) = (c_a)_{j_g} (c_d)_{j_b} = delta_ag delta_db.
    pairing_identity = True
    for a in range(17):
        for b in range(17):
            for g in range(17):
                for d in range(17):
                    val = ((cvec[a] >> jpt[g]) & 1) & ((cvec[d] >> jpt[b]) & 1)
                    if val != (1 if (a == g and b == d) else 0):
                        pairing_identity = False
    check("289x289_logical_pairing_is_identity", pairing_identity)
    return {
        "k_logical": 289, "seed_rank": rk, "kernel_dimension": 17,
        "query_family": "u_(a,b) = (c_a e_{j_b}^T, 0), c_a = 1_{P_a}+1_{P_inf}, j_b = recta y=b*x",
        "query_weight_each": 34,
        "dual_family": "L_(g,d) = (e_{j_g} c_d^T, 0)",
        "pairing": "identidad 289x289: las 289 consultas son base de ker(H_Z)/fila(H_X)",
        "preprocessing": "una eliminacion binaria 289x306 (rango de H_0); las consultas son explicitas",
    }


# ================================================================ parte B
def part_b():
    n, mx = 10, 3
    P, Q = Fraction(1, 20), Fraction(1, 50)
    beta = Fraction(1, 2)
    cc = 9
    hx = [sum(1 << (3 * i + j) for j in range(3)) | 1 << cc for i in range(3)]
    hz = [sum(1 << (3 * i + j) for i in range(3)) | 1 << cc for j in range(3)]
    check("small_code_css_commutation", all(par(a & b) == 0 for a in hx for b in hz))

    def synd(v):
        return sum(par(v & hx[i]) << i for i in range(mx))

    # familia general de candidatos: subconjuntos no vacios de soportes de generadores
    cand = {}
    for g in hz:
        bits = [b for b in range(n) if g >> b & 1]
        for r in range(1, len(bits) + 1):
            for sub in combinations(bits, r):
                f = sum(1 << b for b in sub)
                phi = synd(f)
                if phi:
                    cand[f] = phi
    order = sorted(cand, key=lambda f: (f.bit_count(), f))

    def decode(s):
        corr = 0
        while True:
            best = None
            for f in order:
                phi = cand[f]
                sc = Fraction(2 * (s & phi).bit_count() - phi.bit_count(), phi.bit_count())
                if best is None or sc > best[0]:
                    best = (sc, f, phi)
            if best[0] < beta:
                return corr
            s ^= best[2]
            corr ^= best[1]

    ctab = [decode(s) for s in range(1 << mx)]
    all_v = range(1 << n)
    t_space = [v for v in all_v if synd(v) == 0]
    kerz = [v for v in all_v if all(par(v & g) == 0 for g in hz)]

    def span(vs):
        out = {0}
        for v in vs:
            out |= {x ^ v for x in out}
        return out

    s_x = span(hz)
    row_hx = span(hx)
    check("small_code_dimensions", len(t_space) == 128 and len(kerz) == 128 and len(s_x) == 8 and len(row_hx) == 8)
    # Q: base de ker H_Z (d = 7).  Lambda: base de ker H_Z modulo fila(H_X) (k = 4).
    qrows, cur = [], {0}
    for v in kerz:
        if v not in cur:
            qrows.append(v)
            cur |= {x ^ v for x in cur}
    lam, cur = [], set(row_hx)
    for v in kerz:
        if v not in cur:
            lam.append(v)
            cur |= {x ^ v for x in cur}
    check("small_code_d7_k4", len(qrows) == 7 and len(lam) == 4 and len(cur) == 128)

    def img(rows, v):
        return sum(par(v & r) << i for i, r in enumerate(rows))

    check("on_T_membership_in_S_iff_Lambda_R_zero", all((v in s_x) == (img(lam, v) == 0) for v in t_space))
    check("on_T_membership_in_S_iff_Q_R_zero", all((v in s_x) == (img(qrows, v) == 0) for v in t_space))
    check("Lambda_separates_all_16_classes", len({img(lam, v) for v in t_space}) == 16)

    # todas las historias
    b0 = (1 - P) ** n * (1 - Q) ** mx
    zp, zq = P / (1 - P), Q / (1 - Q)
    p0 = Fraction(0)
    psil = Fraction(0)
    low_sector = Fraction(0)
    mass_k = {}
    mass_d = {}
    hist = 0
    for e in all_v:
        se = synd(e)
        we = zp ** e.bit_count()
        for dn in range(1 << mx):
            hist += 1
            r = e ^ ctab[se ^ dn]
            if synd(r):
                continue
            w = b0 * we * zq ** dn.bit_count()
            p0 += w
            xk = img(lam, r)
            xd = img(qrows, r)
            mass_k[xk] = mass_k.get(xk, 0) + w
            mass_d[xd] = mass_d.get(xd, 0) + w
            if xk:
                psil += w
                if e.bit_count() + dn.bit_count() <= 1:
                    low_sector += w
    qsil = psil / p0
    check("all_8192_histories_enumerated", hist == 8192)
    check("silent_failure_is_positive", 0 < psil < p0 <= 1)
    check("zero_class_agrees_between_quotients", mass_k.get(0, 0) == mass_d.get(0, 0))

    def lam_k(s):
        return sum(w * (-1) ** par(s & x) for x, w in mass_k.items()) / p0

    def lam_d(s):
        return sum(w * (-1) ** par(s & x) for x, w in mass_d.items()) / p0

    avg_k = sum(1 - lam_k(s) for s in range(16)) / 16
    avg_d = sum(1 - lam_d(s) for s in range(128)) / 128
    check("exact_average_over_logical_quotient_equals_Qsil", avg_k == qsil)
    check("exact_average_over_full_quotient_equals_Qsil", avg_d == qsil)

    # conjunto de sesgo pequeno explicito en F_2^4: GF(2^10), x^10 + x^3 + 1
    mgf = aghp_field_bits(4)
    check("small_code_field_bits", mgf == 10)
    poly = (1 << 10) | (1 << 3) | 1

    def gmul(a, b):
        r = 0
        while b:
            if b & 1:
                r ^= a
            b >>= 1
            a <<= 1
            if a >> 10:
                a ^= poly
        return r

    mult = [0] * 16
    for x in range(1 << mgf):
        pw = [1, x, gmul(x, x)]
        pw.append(gmul(pw[2], x))
        for y in range(1 << mgf):
            s = 0
            for i in range(4):
                s |= par(pw[i] & y) << i
            mult[s] += 1
    size = 1 << (2 * mgf)
    bias = max(abs(sum(mult[s] * (-1) ** par(a & s) for s in range(16))) for a in range(1, 16))
    bias = Fraction(bias, size)
    check("explicit_multiset_bias_at_most_3_over_1024", bias <= Fraction(3, 1 << mgf))
    check("explicit_multiset_bias_within_budget", bias <= DELTA_SB)
    lam_true = [lam_k(s) for s in range(16)]
    f_true = sum(mult[s] * (1 - lam_true[s]) for s in range(16)) / size
    check("small_bias_sandwich", (1 - DELTA_SB) * qsil <= f_true <= (1 + DELTA_SB) * qsil)

    # menorantes
    wmin = min(v.bit_count() for v in t_space if v not in s_x)
    ell_w = zp ** wmin / (1 + zp ** wmin)
    ell_new = max(ell_w, low_sector)            # U_0 = 1
    check("witness_lower_bound_valid", ell_w <= qsil)
    check("sector_lower_bound_valid", 0 < low_sector <= psil)
    check("improved_lower_bound_valid_and_larger", ell_w < ell_new <= qsil)

    # presupuesto con errores adversos de acceso
    rng = random.Random(20261111)
    xi_den = p0 / b0
    xi_num = [lam_true[s] * xi_den for s in range(16)]
    worst = {}
    for label, ell in (("witness", ell_w), ("improved", ell_new)):
        bd = budget(ell)
        patterns = [([1] * 16, -1), ([-1] * 16, 1), ([1] * 16, 1), ([-1] * 16, -1)]
        for _ in range(400):
            patterns.append(([rng.choice((-1, 1)) for _ in range(16)], rng.choice((-1, 1))))
        wq = Fraction(0)
        wp = Fraction(0)
        w0 = Fraction(0)
        for sg, s0 in patterns:
            den_hat = xi_den * (1 + s0 * bd["zeta_den"])
            acc = Fraction(0)
            for s in range(1, 16):
                if mult[s]:
                    acc += mult[s] * (1 - xi_num[s] * (1 + sg[s] * bd["zeta_num"]) / den_hat)
            q_hat = min(max(acc / size, Fraction(0)), Fraction(1))      # mascara cero devuelta exacta
            p0_hat = b0 * den_hat
            wq = max(wq, abs(q_hat / qsil - 1))
            w0 = max(w0, abs(p0_hat / p0 - 1))
            wp = max(wp, abs(q_hat * p0_hat / psil - 1))
        check("budget_%s_Qsil_within_e" % label, wq <= E_BUD)
        check("budget_%s_P0_within_e" % label, w0 <= E_BUD)
        check("budget_%s_Psil_within_epsilon" % label, wp <= EPS)
        worst[label] = {"zeta_num": both(bd["zeta_num"]), "worst_relative_error_Qsil": both(wq),
                        "worst_relative_error_Psil": both(wp), "sign_patterns": len(patterns)}
    return {
        "code": "HGP de la semilla [1 1 1]: n=10, m_X=3, k=4, d=n-dim S_X=7",
        "decoder": "SSF, todos los subconjuntos de generadores con huella no vacia, beta=1/2, orden (|F|,mask)",
        "rates": {"p": str(P), "q": str(Q)}, "histories": hist, "candidates": len(order),
        "P0": both(p0), "P_sil": both(psil), "Q_sil": both(qsil),
        "queries_full_quotient_dimension": 7, "queries_logical_quotient_dimension": 4,
        "explicit_small_bias": {"field_bits": mgf, "multiset_size": size, "exact_bias": str(bias),
                                "F_over_Qsil": both(f_true / qsil)},
        "lower_bounds": {"witness_weight": wmin, "witness": both(ell_w),
                         "sector_at_most_one_fault": both(low_sector), "improved": both(ell_new)},
        "budget_with_adversarial_access_errors": worst,
    }


# ================================================================ parte C
def syndrome(line_set_mask):
    s = 0
    j = 0
    while line_set_mask:
        if line_set_mask & 1:
            s ^= LM[j]
        line_set_mask >>= 1
        j += 1
    return s


def decode_row(s):
    corr = 0
    lm = LM
    while True:
        best = -1
        bj = -1
        for j in range(306):
            h = (s & lm[j]).bit_count()
            if h > best:
                best = h
                bj = j
        if best < 9:
            return corr, s
        s ^= lm[bj]
        corr ^= 1 << bj


def part_c():
    q_lines = CLASSES[0]
    q_mask = 0
    for j in q_lines:
        q_mask |= 1 << j
    runs = 0
    silent = 0
    for ci in range(1, 18):
        cl = CLASSES[ci]
        full = 0
        for j in cl:
            full |= 1 << j
        target = full ^ q_mask                      # residual P xor Q
        for k1, k2 in combinations(cl, 2):
            e = full ^ (1 << k1) ^ (1 << k2)
            s_e = syndrome(e)
            kpts = LM[k1] | LM[k2]
            for j1, j2 in combinations(q_lines, 2):
                d = kpts & (LM[j1] | LM[j2])
                corr, s_final = decode_row(s_e ^ d)
                runs += 1
                if s_final == d and (e ^ corr) == target:
                    silent += 1
    check("rectangle_run_count", runs == 17 * 136 * 136)
    check("rectangle_15_4_all_silent_with_residual_P_xor_Q", silent == runs)
    check("rectangle_15_4_total_matches_inherited_N2", 306 * silent == 96216192)
    return {"description": "E = clase P menos 2 rectas, D = 4 cruces con 2 rectas de la primera clase; residual P xor Q",
            "row_decoder_runs": runs, "silent_per_row": silent, "silent_all_rows": 306 * silent,
            "inherited_source": "ETAPA2R1_PARTICION_REVISOR_20261104/CONTROL_RECTANGULOS.json (N_2)"}


# ================================================================ parte D
N_DATA, M_X, N_TOT, N_HOST = 177157, 88434, 265591, 11423881
THETA = Fraction(17, 361700)


def corner_row(ell, n_host):
    bd = budget(ell)
    zeta = bd["zeta_num"]
    k2 = ceil_log2(8 * n_host / zeta)
    m_trunc = -((-k2 * THETA.denominator) // THETA.numerator)       # ceil(k2/theta)
    return {"lower_bound_on_Qsil": both(ell), "log10_lower_bound": round(log10_fr(ell), 3),
            "relative_tolerance_per_partition": both(zeta),
            "bits_of_relative_precision": ceil_log2(1 / zeta),
            "KP_truncation_order_M": m_trunc}


def part_d(n154):
    nu = Fraction(306, 182451625)
    y = Fraction(1250, 74320107)
    x_e = Fraction(2850806640625, 598978202497210701)
    zs = 2 ** 14 * nu * x_e * y ** 15 / 3
    ps = zs / (1 + zs)
    qs = Fraction(51, 182451625000000000051)
    zq = qs / (1 - qs)
    check("p_star_reading", abs(ps - Fraction(1062880453202439, 10 ** 94)) < Fraction(1, 10 ** 94))
    check("q_star_odds_equal_nu_xD_over_3", zq == nu * Fraction(1, 2 * 10 ** 12) / 3)

    # ---- esquina fenomenologica
    ell_w = zs ** 34 / (1 + zs ** 34)
    absent = 1 - N_DATA * ps - M_X * qs                      # <= B_0
    check("absence_bound_positive", 0 < absent < 1)
    sectors = {"(17,0)_5202": 5202 * zs ** 17, "(16,1)_1503378": 1503378 * zs ** 16 * zq,
               "(15,4)_%d" % n154: n154 * zs ** 15 * zq ** 4}
    l_new = absent * sum(sectors.values())
    ell_new = max(ell_w, l_new)                              # U_0 = 1
    check("phenomenological_improved_exceeds_witness", ell_new > ell_w)
    phen = {"point": {"p_star": both(ps), "q_star": both(qs)},
            "before": corner_row(ell_w, N_TOT), "after": corner_row(ell_new, N_TOT),
            "sector_masses_over_B0": {k: both(v) for k, v in sectors.items()},
            "orders_gained_in_lower_bound": round(log10_fr(ell_new / ell_w), 3)}

    # ---- esquina de circuito
    t_d = zs / 55
    t_c = zq / 19
    xi = t_d / 2
    rc = xi / (1 + xi)
    rm = (t_c - xi) / (1 + t_c)
    check("circuit_corner_identities", rc / (1 - rc) == t_d / 2 and 1 / ((1 - rc) * (1 - rm)) - 1 == t_c)
    kap = (2 * rc / (3 - 2 * rc)) ** 34
    ell_l = kap / (1 + kap)
    alpha = 2 * rc / 3
    betac = 8 * rc / 15
    zw = alpha / (1 - alpha)
    zc = betac / (1 - betac)
    w0_lo = 1 - 6471288 * alpha - 3095190 * betac            # Bernoulli: (1-a)^N >= 1-Na
    one_m2t = (1 - 4 * rc / 3) ** 19 * ((15 - 16 * rc) / (15 - 8 * rc)) ** 35 * (1 - 2 * rm)
    th = (1 - one_m2t) / 2
    zmix_lo = w0_lo * (17 * zw + zc) ** 14 * (2 * zw + zc)
    rho17 = (17 * zw * th + zc / 2) / (17 * zw + zc)
    rho2 = (2 * zw * th + zc / 2) / (2 * zw + zc)
    q0 = (1 - rho17) ** 14 * (1 - rho2)
    x = rho17 / (1 - rho17)
    yy = rho2 / (1 - rho2)
    pv = 1 + 2 * x + 21 * x ** 2 + 100 * x ** 3 + 66 * x ** 4 + 220 * x ** 5 \
        + yy * (3 * x + 33 * x ** 2 + 78 * x ** 3 + 286 * x ** 4)
    lfin_lo = 306 * zmix_lo * th ** 3 * (1 - (M_X - 18) * th) * q0 * pv
    check("circuit_witness_matches_plan_3.2064e-2761",
          abs(log10_fr(ell_l) - math.log10(3.20638882424e-300) + 2461) < 1e-6)
    check("L_fin_at_least_manuscript_value", lfin_lo >= Fraction(798339561055, 10 ** 11) / Fraction(10) ** 1258)
    check("L_fin_below_next_digit", lfin_lo < Fraction(798339561057, 10 ** 11) / Fraction(10) ** 1258)
    ell_c_new = max(ell_l, lfin_lo)
    circ = {"point": {"R_c": both(rc), "R_m": both(rm)},
            "before": corner_row(ell_l, N_HOST), "after": corner_row(ell_c_new, N_HOST),
            "L_fin_lower": both(lfin_lo),
            "orders_gained_in_lower_bound": round(log10_fr(ell_c_new / ell_l), 3)}

    # ---- numero de mascaras
    d_full = N_DATA - M_X
    check("F17_quotient_dimensions", d_full == 88723 and d_full - M_X == 289)
    mb, ma = aghp_field_bits(d_full), aghp_field_bits(289)
    queries = {"delta_small_bias": str(DELTA_SB),
               "before": {"dimension": d_full, "field_bits": mb, "masks": "2^%d" % (2 * mb)},
               "after": {"dimension": 289, "field_bits": ma, "masks": "2^%d" % (2 * ma)},
               "ratio": "2^%d" % (2 * (mb - ma)),
               "asymptotic_ratio_of_dimensions_for_17_18_regular_seeds": "k/d = (1-r)^2/(1-r+r^2) = 1/307 con r=17/18",
               "construction": "potencias sobre GF(2^m), sesgo <= (dim-1)/2^m; el manuscrito cita Naor-Naor (poly)"}
    check("dimension_ratio_formula_1_over_307",
          Fraction(1, 18) ** 2 / (1 - Fraction(17, 18) + Fraction(17, 18) ** 2) == Fraction(1, 307))
    return {"phenomenological_corner": phen, "circuit_corner": circ, "queries": queries,
            "per_query_cost_form": "N * exp[O(M log max(2,Delta))] * poly(N, M, L_in), con M el orden de truncamiento",
            "KP_rate_theta": str(THETA)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    if os.path.lexists(args.output):
        raise SystemExit("La salida ya existe; use una ruta nueva.")
    t0 = time.time()
    a = part_a()
    ta = time.time()
    b = part_b()
    tb = time.time()
    c = part_c()
    tc = time.time()
    d = part_d(c["silent_all_rows"])
    td = time.time()
    result = {
        "unit": "H1_COCIENTE_Y_MENORANTE_20261111",
        "status": "PASS",
        "budget": {"epsilon": str(EPS), "e": str(E_BUD), "delta_small_bias": str(DELTA_SB),
                   "eta": "e*ell/4", "zeta_num": "e*ell/16", "zeta_den": "min(e, e*ell/16)"},
        "A_F17_logical_query_basis": a,
        "B_small_code_exhaustive_control": b,
        "C_rectangle_family_reverification": c,
        "D_before_after_at_certified_corners": d,
        "counts": {"named_checks": len(CHECKS), "F17_row_decoder_runs": c["row_decoder_runs"],
                   "small_code_histories": b["histories"], "full_FGL_trajectories_F17": 0},
        "checks": CHECKS,
    }
    with open(args.output, "x", encoding="utf-8", newline="\n") as f:
        json.dump(result, f, ensure_ascii=False, indent=1, sort_keys=True)
        f.write("\n")
    print("PASS controles=%d  A=%.1fs B=%.1fs C=%.1fs D=%.1fs" % (len(CHECKS), ta - t0, tb - ta, tc - tb, td - tc))


if __name__ == "__main__":
    main()
