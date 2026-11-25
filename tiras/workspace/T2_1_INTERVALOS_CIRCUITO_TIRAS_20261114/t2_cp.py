#!/usr/bin/env python3
"""T2-1 (BORRADOR_PLAN_T2_1.txt, U4/D-mc): intervalo de Clopper-Pearson bilateral con evaluacion conservadora.

Codigo preparado para revision (version 2: valida alfa y objetivo y declara XMAX en cp_interval, segun el
dictamen DICTAMEN_CODIGO_T2_1.txt). No se ha usado en ninguna corrida de T2-1.

Contrato:
- X ~ Bin(N, p) con N fijado antes. Con x observadas y alfa = 0.001 (5e-4 en cada cola):
    p_inf resuelve Pr_p[X >= x] = alfa/2  (p_inf = 0 si x = 0),
    p_sup resuelve Pr_p[X <= x] = alfa/2  (p_sup = 1 si x = N).
- Las colas se evaluan en aritmetica de intervalos con Decimal de PREC cifras: sumas, productos y cocientes con
  redondeo dirigido (ROUND_FLOOR para la cota inferior, ROUND_CEILING para la superior). ln y exp, que Decimal
  redondea correctamente al mas cercano (error <= medio ulp), se ensanchan dos ulps hacia fuera. Todos los
  terminos son positivos y la recursion t_(j+1) = t_j (N - j)/(j + 1) p/(1 - p) es monotona en cada factor.
- Pr_p[X >= x] crece con p (x >= 1) y Pr_p[X <= x] decrece con p (x < N). La biseccion solo avanza cuando el
  intervalo certifica el signo. Para p_inf se devuelve el mayor punto a con cota superior de Pr_a[X >= x] <= alfa/2
  (la raiz es >= a); para p_sup, el menor punto b con cota superior de Pr_b[X <= x] <= alfa/2 (la raiz es <= b).
  Si el intervalo no decide, la biseccion para y devuelve el extremo ya certificado: mas ancho, nunca mas estrecho.
- Garwood (Poisson) es solo un control de aproximacion en la autoprueba; no sustituye la cola binomial.

Uso (desde cualquier directorio; no lee ni escribe ficheros del proyecto salvo --output):
  python -I -B t2_cp.py --autoprueba [--output <nuevo.json>]
  python -I -B t2_cp.py --x <recuento> --n <historias> [--alfa 0.001]
"""
import argparse
import json
import math
import os
import sys
import time
from decimal import ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN, Context, Decimal
from fractions import Fraction as F

PREC = 60
ELIM = 999999
XMAX = 10 ** 6
ULPS = 2
DN = Context(prec=PREC, rounding=ROUND_FLOOR, Emin=-ELIM, Emax=ELIM)
UP = Context(prec=PREC, rounding=ROUND_CEILING, Emin=-ELIM, Emax=ELIM)
NE = Context(prec=PREC, rounding=ROUND_HALF_EVEN, Emin=-ELIM, Emax=ELIM)
ZERO, ONE, TWO = Decimal(0), Decimal(1), Decimal(2)


def _down(v):
    for _ in range(ULPS):
        v = NE.next_minus(v)
    return v


def _up(v):
    for _ in range(ULPS):
        v = NE.next_plus(v)
    return v


def _check_p(p):
    if not isinstance(p, Decimal) or not p.is_finite() or p < 0 or p > 1:
        raise ValueError("p debe ser un Decimal finito en [0, 1]")


def _pow_interval(base_lo, base_hi, n):
    """Cotas de b^n para 0 < base_lo <= b <= base_hi <= 1 y n entero >= 1, via exp(n ln b)."""
    l_lo = _down(NE.ln(base_lo))
    l_hi = _up(NE.ln(base_hi))
    a_lo = DN.multiply(Decimal(n), l_lo)
    a_hi = UP.multiply(Decimal(n), l_hi)
    lo = _down(NE.exp(a_lo))
    hi = _up(NE.exp(a_hi))
    return max(lo, ZERO), min(hi, ONE)


def cdf_interval(x, n, p):
    """Cotas [lo, hi] de Pr[Bin(n, p) <= x]."""
    _check_p(p)
    if n < 1 or x < 0:
        raise ValueError("n >= 1 y x >= 0")
    if x >= n:
        return ONE, ONE
    if x > XMAX:
        raise ValueError("x fuera del alcance del evaluador (XMAX)")
    if p == 0:
        return ONE, ONE
    if p == 1:
        return ZERO, ZERO
    q_lo, q_hi = DN.subtract(ONE, p), UP.subtract(ONE, p)
    t_lo, t_hi = _pow_interval(q_lo, q_hi, n)
    r_lo, r_hi = DN.divide(p, q_hi), UP.divide(p, q_lo)
    s_lo, s_hi = t_lo, t_hi
    for j in range(x):
        f_lo = DN.divide(Decimal(n - j), Decimal(j + 1))
        f_hi = UP.divide(Decimal(n - j), Decimal(j + 1))
        t_lo = DN.multiply(DN.multiply(t_lo, f_lo), r_lo)
        t_hi = UP.multiply(UP.multiply(t_hi, f_hi), r_hi)
        s_lo = DN.add(s_lo, t_lo)
        s_hi = UP.add(s_hi, t_hi)
    return max(s_lo, ZERO), min(s_hi, ONE)


def tail_interval(x, n, p):
    """Cotas [lo, hi] de Pr[Bin(n, p) >= x]."""
    _check_p(p)
    if x <= 0:
        return ONE, ONE
    if x > n:
        return ZERO, ZERO
    if p == 0:
        return ZERO, ZERO
    if p == 1:
        return ONE, ONE
    if x == n:
        return _pow_interval(p, p, n)
    c_lo, c_hi = cdf_interval(x - 1, n, p)
    return max(DN.subtract(ONE, c_hi), ZERO), min(UP.subtract(ONE, c_lo), ONE)


def _mid(a, b):
    return NE.divide(NE.add(a, b), TWO)


def _bisect(fun, a, b, below_is_a, half, target):
    """fun(p) -> (lo, hi). Invariantes: si below_is_a, la cota superior en a es <= half (a por debajo de la raiz
    de una funcion creciente) y la inferior en b es > half; si no, la inferior en a es > half y la superior en b
    es <= half (funcion decreciente). Devuelve (a, b, pasos, detenida)."""
    steps = 0
    stopped = False
    while UP.subtract(b, a) > DN.multiply(target, b):
        c = _mid(a, b)
        if not a < c < b:
            stopped = True
            break
        lo, hi = fun(c)
        steps += 1
        if below_is_a:
            if hi <= half:
                a = c
            elif lo > half:
                b = c
            else:
                stopped = True
                break
        else:
            if lo > half:
                a = c
            elif hi <= half:
                b = c
            else:
                stopped = True
                break
    return a, b, steps, stopped


def cp_lower(x, n, half, target):
    if x == 0:
        return {"valor": ZERO, "raiz_menor_que": ZERO, "pasos": 0, "detenida": False}
    fun = lambda p: tail_interval(x, n, p)  # noqa: E731
    a = ZERO                                   # Pr_0[X >= x] = 0 para x >= 1
    b = min(ONE, NE.divide(Decimal(x), Decimal(n)))
    while True:
        lo, hi = fun(b)
        if lo > half:
            break
        if hi <= half:
            a = b
        if b == ONE:
            raise RuntimeError("no se certifica Pr_1[X >= x] > alfa/2")
        b = min(ONE, NE.multiply(b, TWO))
    a, b, steps, stopped = _bisect(fun, a, b, True, half, target)
    return {"valor": a, "raiz_menor_que": b, "pasos": steps, "detenida": stopped}


def cp_upper(x, n, half, target):
    if x == n:
        return {"valor": ONE, "raiz_mayor_que": ONE, "pasos": 0, "detenida": False}
    fun = lambda p: cdf_interval(x, n, p)  # noqa: E731
    a = NE.divide(Decimal(x), Decimal(n))
    while not fun(a)[0] > half:              # Pr_0[X <= x] = 1 > alfa/2
        a = NE.divide(a, TWO) if a > Decimal("1e-300") else ZERO
    b = min(ONE, NE.multiply(NE.add(a, NE.divide(Decimal(x + 1), Decimal(n))), TWO))
    while True:
        lo, hi = fun(b)
        if hi <= half:
            break
        if lo > half:
            a = b
        if b == ONE:
            raise RuntimeError("no se certifica Pr_1[X <= x] <= alfa/2")
        b = min(ONE, NE.multiply(b, TWO))
    a, b, steps, stopped = _bisect(fun, a, b, False, half, target)
    return {"valor": b, "raiz_mayor_que": a, "pasos": steps, "detenida": stopped}


def cp_interval(x, n, alpha="0.001", target="1e-12"):
    """Intervalo de Clopper-Pearson bilateral conservador para x silenciosas en n historias."""
    if not (isinstance(x, int) and isinstance(n, int) and 0 <= x <= n and n >= 1):
        raise ValueError("se requieren enteros 0 <= x <= n, n >= 1")
    if x > XMAX and x != n:
        raise ValueError("x = %d fuera del alcance del evaluador (XMAX = %d)" % (x, XMAX))
    al, tg = Decimal(alpha), Decimal(target)
    if not (al.is_finite() and ZERO < al < ONE and tg.is_finite() and ZERO < tg < ONE):
        raise ValueError("se requieren 0 < alfa < 1 y 0 < objetivo < 1, finitos")
    half = NE.divide(al, TWO)
    lo = cp_lower(x, n, half, tg)
    hi = cp_upper(x, n, half, tg)
    return {"x": x, "N": n, "alfa": str(alpha), "p_inf": str(lo["valor"]), "p_sup": str(hi["valor"]),
            "p_inf_raiz_menor_que": str(lo["raiz_menor_que"]), "p_sup_raiz_mayor_que": str(hi["raiz_mayor_que"]),
            "pasos": [lo["pasos"], hi["pasos"]], "biseccion_detenida": [lo["detenida"], hi["detenida"]],
            "precision_cifras": PREC, "ulps_ensanche": ULPS, "anchura_relativa_objetivo": str(target)}


def intersects(cp, lo, hi):
    """D-mc: el intervalo [p_inf, p_sup] corta el intervalo certificado [lo, hi] (fracciones)?"""
    a, b = F(Decimal(cp["p_inf"])), F(Decimal(cp["p_sup"]))
    return a <= hi and lo <= b


# ------------------------------------------------------------------ autoprueba (matematica pura)
def _cdf_exact(x, n, p):
    if x >= n:
        return F(1)
    return sum(F(math.comb(n, j)) * p ** j * (1 - p) ** (n - j) for j in range(x + 1))


def _tail_exact(x, n, p):
    if x <= 0:
        return F(1)
    return 1 - _cdf_exact(x - 1, n, p)


def _poisson_cdf_interval(x, lam):
    """Cotas de Pr[Poisson(lam) <= x] (solo para el control de Garwood)."""
    if lam == 0:
        return ONE, ONE
    neg = NE.minus(lam)                      # exacto: lam tiene <= PREC cifras
    t_lo = _down(NE.exp(neg))
    t_hi = _up(NE.exp(neg))
    s_lo, s_hi = t_lo, t_hi
    for j in range(x):
        t_lo = DN.multiply(t_lo, DN.divide(lam, Decimal(j + 1)))
        t_hi = UP.multiply(t_hi, UP.divide(lam, Decimal(j + 1)))
        s_lo = DN.add(s_lo, t_lo)
        s_hi = UP.add(s_hi, t_hi)
    return max(s_lo, ZERO), min(s_hi, ONE)


def _garwood(x, half, target):
    """Extremos de Garwood para la media de Poisson (control de aproximacion)."""
    out = []
    if x == 0:
        out.append(ZERO)
    else:
        def tail(lam):
            c = _poisson_cdf_interval(x - 1, lam)
            return max(DN.subtract(ONE, c[1]), ZERO), min(UP.subtract(ONE, c[0]), ONE)
        a, b = ZERO, Decimal(x)
        while not tail(b)[0] > half:
            b = NE.multiply(b, TWO)
        out.append(_bisect(tail, a, b, True, half, target)[0])
    fun = lambda lam: _poisson_cdf_interval(x, lam)  # noqa: E731
    a, b = Decimal(x), Decimal(x + 1)
    while not fun(a)[0] > half:
        a = NE.divide(a, TWO)
    while not fun(b)[1] <= half:
        b = NE.multiply(b, TWO)
    out.append(_bisect(fun, a, b, False, half, target)[1])
    return out


def autoprueba():
    bad = []
    log = []
    t0 = time.perf_counter()
    half = Decimal("0.0005")
    target = Decimal("1e-12")

    # 1. Colas con n pequeno: la fraccion exacta cae dentro del intervalo.
    ps = [Decimal(s) for s in ("0", "0.001", "0.1", "0.123456789", "0.3333", "0.5", "0.9", "0.999", "1")]
    ncmp = 0
    for n in (1, 2, 3, 5, 10, 25, 40):
        for p in ps:
            pf = F(p)
            for x in range(0, n + 1):
                lo, hi = cdf_interval(x, n, p)
                ex = _cdf_exact(x, n, pf)
                ncmp += 1
                if not F(lo) <= ex <= F(hi):
                    bad.append(("cdf", n, str(p), x))
                lo, hi = tail_interval(x, n, p)
                ex = _tail_exact(x, n, pf)
                ncmp += 1
                if not F(lo) <= ex <= F(hi):
                    bad.append(("cola", n, str(p), x))
    log.append("1. colas exactas dentro del intervalo: %d comparaciones" % ncmp)

    # 2. Extremos con n pequeno: conservadores y con la raiz exacta encerrada por los dos puntos certificados.
    ncp = 0
    for n in (5, 10, 40):
        for x in sorted({0, 1, n // 2, n - 1, n}):
            cp = cp_interval(x, n)
            ncp += 1
            pi, pis = F(Decimal(cp["p_inf"])), F(Decimal(cp["p_inf_raiz_menor_que"]))
            ps_, pss = F(Decimal(cp["p_sup"])), F(Decimal(cp["p_sup_raiz_mayor_que"]))
            if x > 0:
                if not _tail_exact(x, n, pi) <= F(half) < _tail_exact(x, n, pis):
                    bad.append(("p_inf no encierra la raiz exacta", n, x))
                if (pis - pi) > F(2, 10 ** 12) * pis:
                    bad.append(("p_inf demasiado ancho", n, x))
            elif pi != 0:
                bad.append(("p_inf != 0 con x = 0", n))
            if x < n:
                if not _cdf_exact(x, n, ps_) <= F(half) < _cdf_exact(x, n, pss):
                    bad.append(("p_sup no encierra la raiz exacta", n, x))
                if (ps_ - pss) > F(2, 10 ** 12) * ps_:
                    bad.append(("p_sup demasiado ancho", n, x))
            elif ps_ != 1:
                bad.append(("p_sup != 1 con x = N", n))
            if any(cp["biseccion_detenida"]):
                bad.append(("biseccion detenida con n pequeno", n, x))
    log.append("2. extremos con n pequeno conservadores y ajustados: %d intervalos" % ncp)

    # 3. Formas cerradas: x = 0 (p_sup = 1 - h^(1/N)) y x = N (p_inf = h^(1/N)).
    ctx = Context(prec=90, rounding=ROUND_HALF_EVEN, Emin=-ELIM, Emax=ELIM)
    for n in (10, 1000, 7 * 10 ** 9):
        root = ctx.exp(ctx.divide(ctx.ln(half), Decimal(n)))
        cf_sup = ctx.subtract(ONE, root)
        cp = cp_interval(0, n)
        ps_ = Decimal(cp["p_sup"])
        rel = abs(F(ps_) - F(cf_sup)) / F(cf_sup)
        if not (F(ps_) >= F(cf_sup) * (1 - F(1, 10 ** 40)) and rel <= F(2, 10 ** 12)):
            bad.append(("forma cerrada x = 0", n, str(ps_), str(cf_sup)))
        log.append("3. x = 0, N = %d: p_sup = %s (forma cerrada %s)" % (n, ps_, ctx.plus(cf_sup).normalize()))
        if n <= 1000:
            cp = cp_interval(n, n)
            pi = Decimal(cp["p_inf"])
            rel = abs(F(pi) - F(root)) / F(root)
            if not (F(pi) <= F(root) * (1 + F(1, 10 ** 40)) and rel <= F(2, 10 ** 12)):
                bad.append(("forma cerrada x = N", n, str(pi), str(root)))
            log.append("3. x = N = %d: p_inf = %s (forma cerrada %s)" % (n, pi, ctx.plus(root).normalize()))

    # 4. Garwood con N grande (control de aproximacion, no sustituye la cola binomial).
    n = 7 * 10 ** 9
    for x in (1, 10, 385, 1000):
        cp = cp_interval(x, n)
        g = _garwood(x, half, target)
        lo_b, hi_b = F(Decimal(cp["p_inf"])), F(Decimal(cp["p_sup"]))
        lo_g, hi_g = F(g[0]) / n, F(g[1]) / n
        r_lo = abs(lo_b - lo_g) / lo_g
        r_hi = abs(hi_b - hi_g) / hi_g
        if r_lo > F(1, 10 ** 6) or r_hi > F(1, 10 ** 6):
            bad.append(("Garwood", x, float(r_lo), float(r_hi)))
        if any(cp["biseccion_detenida"]):
            bad.append(("biseccion detenida con N grande", x))
        log.append("4. N = 7e9, x = %d: [%.12e, %.12e]; Garwood/N [%.12e, %.12e]; dif. rel. %.2e y %.2e; pasos %s"
                   % (x, float(lo_b), float(hi_b), float(lo_g), float(hi_g), float(r_lo), float(r_hi), cp["pasos"]))

    # 5. Interseccion con un intervalo certificado.
    cp = cp_interval(385, n)
    a = F(Decimal(cp["p_inf"]))
    if not intersects(cp, a / 2, a) or intersects(cp, a / 4, a / 2):
        bad.append(("interseccion",))
    log.append("5. interseccion de intervalos comprobada en dos casos")

    # 6. Entradas fuera de contrato: se rechazan antes de buscar raices.
    nrej = 0
    for args in ((1, 10, "0"), (1, 10, "1"), (1, 10, "-0.1"), (1, 10, "NaN"), (1, 10, "Infinity"),
                 (1, 10, "0.001", "0"), (1, 10, "0.001", "1"), (1, 10, "0.001", "NaN"),
                 (XMAX + 1, 10 ** 9), (-1, 10), (11, 10), (0, 0)):
        try:
            cp_interval(*args)
            bad.append(("entrada no rechazada",) + args)
        except ValueError:
            nrej += 1
    log.append("6. entradas fuera de contrato rechazadas: %d" % nrej)
    return {"fallos": [list(map(str, b)) for b in bad], "lineas": log,
            "segundos": round(time.perf_counter() - t0, 3), "precision_cifras": PREC, "ulps_ensanche": ULPS}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--autoprueba", action="store_true")
    ap.add_argument("--x", type=int)
    ap.add_argument("--n", type=int)
    ap.add_argument("--alfa", default="0.001")
    ap.add_argument("--output")
    a = ap.parse_args()
    if a.output and os.path.exists(a.output):
        raise SystemExit("--output ya existe")
    if a.autoprueba:
        res = autoprueba()
        for line in res["lineas"]:
            print(line)
        print("FALLOS:", res["fallos"] if res["fallos"] else "ninguno")
        code = 1 if res["fallos"] else 0
    else:
        if a.x is None or a.n is None:
            raise SystemExit("faltan --x y --n")
        res = cp_interval(a.x, a.n, a.alfa)
        print(json.dumps(res, indent=1))
        code = 0
    if a.output:
        with open(a.output, "w", encoding="ascii", newline="\n") as fo:
            fo.write(json.dumps(res, ensure_ascii=True, indent=1) + "\n")
    raise SystemExit(code)


if __name__ == "__main__":
    main()
