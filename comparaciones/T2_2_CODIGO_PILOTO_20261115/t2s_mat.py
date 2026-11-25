#!/usr/bin/env python3
"""T2-2, piloto m = 6 (BORRADOR_PLAN_T2_2.txt v3): modulo matematico.

Codigo preparado para revision (version 1). No se ha ejecutado ninguna corrida: solo las autopruebas y las
pruebas de humo sinteticas de t2s_pruebas.py, sin decoder, DP, ancla ni cadenas fisicas.

Contenido (seccion del plan entre parentesis):
- representacion racional de las tasas y escalera (2.5; C2);
- leyes de sitio, soporte (V3), tabla de claves y firmas (2.3);
- tablas de logaritmos de los cocientes de pesos (2.9; C2: sin subdesbordamiento silencioso);
- histogramas por bloque y estadisticos exactos sobre los valores en coma flotante (2.11; R5);
- Bennett con medias (2.9), criterio de FF (2.8) y diagnosticos d2 a d8 con las reglas degeneradas (3.5; C1);
- cuantil de Student con precision declarada (2.9; C2);
- ancla: parada en el exito s, estimador e intervalo por inversion de la binomial negativa (2.6; R2; C3);
- estimacion, medidas de precision, decisiones D-s1 a D-s7, veredicto y regla de ampliacion (2.10, 7; C1, C2).

Tratamiento numerico (C2), fijado antes de calcular:
- Probabilidades, leyes, tasas, odds, pi_j(vacio), intervalos de referencia y decisiones: fracciones exactas.
- Constantes irracionales (ln, sqrt, exp, cuantil de Student): Decimal de DEC_PREC = 60 cifras con redondeo al
  par y ensanche de 2 ulps en el sentido conservador cuando una cota decide algo.
- Logaritmos de los cocientes de pesos: float64 correctamente redondeado desde 60 cifras; z(E) se suma con
  math.fsum (suma correctamente redondeada: el orden no influye y el resultado es reproducible).
- Medias de Bennett en dominio logaritmico (log-suma-exp con math.fsum): sin subdesbordamiento ni division por
  cero. Un termino que se anularia (exp < 2^-1022) se cuenta como subdesbordamiento: nunca es silencioso y deja
  la cantidad INDETERMINADA.
- Estadisticos del criterio de FF, de las medias por lotes y del R-hat: exactos (enteros) sobre los valores
  float de g, escalados por su maximo. Varianza nula, W = 0 o denominador nulo se detectan sin error de redondeo.
- Ningun valor no finito pasa una comparacion: las funciones lanzan ErrorNumerico o devuelven INDETERMINADO.
"""
import hashlib
import importlib.util
import json
import math
from decimal import ROUND_CEILING, ROUND_FLOOR, ROUND_HALF_EVEN, Context, Decimal
from fractions import Fraction as F

# ------------------------------------------------------------------ parametros del plan (2.12; se sellan)
ENDPOINTS = (F(1, 100), F(1, 1000), F(1, 10 ** 4), F(1, 10 ** 5), F(1, 10 ** 6))
TARGETS = ENDPOINTS[1:]
W_MIN = 3                      # f_s = 3 (V5) en lugar de d/2
RATE_DIGITS = 12               # C2: cifras significativas de las tasas intermedias (techo)
MAX_LEVELS = 400               # guarda contra bucles; no es un parametro del metodo
T_INIT = 100_000
LAMBDA_INI = 2
MAX_EXT = 4
EPS = F(1, 20)
BLOCK = 2000
L_SEEDS = 12
M_REPS = 3
S_ANCHOR = 2010
B_ANCHOR = 1_000_000
HALF_ANCHOR = "0.0125"         # alfa/2 con alfa = 0.025 (97.5 % en el ancla)
TARGET_CP = "1e-12"            # anchura relativa objetivo de la biseccion de t2_cp.py
STUDENT_P = "0.9875"
STUDENT_NU = 11
PREC_THR = F(1, 10)
RHAT_THR = F(11, 10)
BATCHES = 50
D5_K = (3, 4)
D5_Z = 3
BENNETT_ITERS = 3
SEMILLA_MAESTRA = 20261115     # fijada antes de calcular (no figura en la v2; se anade en la v3)

# ------------------------------------------------------------------ numerica (C2)
DEC_PREC = 60
_E = 999999
NE = Context(prec=DEC_PREC, rounding=ROUND_HALF_EVEN, Emin=-_E, Emax=_E)
RATE_CTX = Context(prec=RATE_DIGITS, rounding=ROUND_CEILING, Emin=-_E, Emax=_E)
RATE_MARGIN = Decimal("1e-40")
ULPS = 2
FLOAT_MIN = 2.0 ** -1022       # menor float normal; por debajo se cuenta subdesbordamiento


class ErrorNumerico(ArithmeticError):
    """Valor no finito, denominador nulo o subdesbordamiento: nunca se lee como 0 ni pasa una comparacion."""


def d_of(x):
    """Fraccion o entero -> Decimal de 60 cifras (al par; error <= medio ulp)."""
    x = F(x)
    return NE.divide(Decimal(x.numerator), Decimal(x.denominator))


def d_up(v):
    for _ in range(ULPS):
        v = NE.next_plus(v)
    return v


def d_dn(v):
    for _ in range(ULPS):
        v = NE.next_minus(v)
    return v


def up_of(x):
    """Decimal >= x (fraccion)."""
    return d_up(d_of(x))


def dn_of(x):
    """Decimal <= x (fraccion)."""
    return d_dn(d_of(x))


def fin(x, what):
    if not (isinstance(x, float) and math.isfinite(x)):
        raise ErrorNumerico("%s no finito: %r" % (what, x))
    return x


def to_float(d, what):
    """Decimal -> float correctamente redondeado; no finito o subnormal -> ErrorNumerico."""
    v = float(d)
    if not math.isfinite(v) or (v == 0.0 and d != 0) or (v != 0.0 and abs(v) < FLOAT_MIN):
        raise ErrorNumerico("%s fuera del rango normal de float64: %s" % (what, d))
    return v


def ln_frac(x, what="ln"):
    """ln de una fraccion positiva como float (desde 60 cifras)."""
    x = F(x)
    if x <= 0:
        raise ErrorNumerico("%s de un valor no positivo" % what)
    return to_float(NE.ln(d_of(x)), what)


def softplus(y):
    """log(1 + e^y) sin desbordamiento."""
    fin(y, "argumento de softplus")
    if y > 0:
        return y + math.log1p(math.exp(-y))
    return math.log1p(math.exp(y))


def dstr(x, digits=30):
    """Vista decimal (no decide nada)."""
    if isinstance(x, F):
        x = d_of(x)
    return str(Context(prec=digits, rounding=ROUND_HALF_EVEN, Emin=-_E, Emax=_E).plus(x))


def fstr(x):
    return "%d/%d" % (x.numerator, x.denominator) if x.denominator != 1 else str(x.numerator)


def canon_sha(obj):
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")
    return hashlib.sha256(raw).hexdigest()


# ------------------------------------------------------------------ tasas y escalera (2.5; C2)
LN2 = NE.ln(Decimal(2))


def rule_factor(w):
    """2^(-1/sqrt(w)) en Decimal de 60 cifras."""
    return NE.exp(NE.minus(NE.divide(LN2, NE.sqrt(d_of(w)))))


def next_rate(r, lam):
    """Siguiente tasa de la regla de BV ec. 17 con w = max(3, lambda), representada como el techo a
    RATE_DIGITS cifras significativas del valor de la regla aumentado en 1e-40 relativo. Asi la tasa
    representada es >= el valor exacto de la regla: ningun paso es mayor que el de la regla."""
    w = max(F(W_MIN), F(lam))
    v = NE.multiply(d_of(r), rule_factor(w))
    v = NE.multiply(v, NE.add(Decimal(1), RATE_MARGIN))
    return F(RATE_CTX.plus(v)), w


def build_ladder(lam_of, ends=ENDPOINTS, max_levels=MAX_LEVELS):
    """Escalera r_0 > ... > r_t por segmentos; lam_of(r) da lambda exacto (fraccion). El punto que pasaria del
    extremo inferior de un segmento se sustituye por ese extremo. t sale de la escalera efectiva."""
    ends = [F(e) for e in ends]
    if len(ends) < 2 or ends[-1] <= 0 or any(a <= b for a, b in zip(ends, ends[1:])):
        raise ValueError("los extremos deben ser positivos y estrictamente decrecientes")
    rates = [ends[0]]
    for hi, lo in zip(ends[:-1], ends[1:]):
        if rates[-1] != hi:
            raise ValueError("los extremos deben ser decrecientes y encadenados")
        r = hi
        while True:
            nxt, _ = next_rate(r, lam_of(r))
            if not 0 < nxt < r:
                raise ValueError("la regla no avanza en r = %s" % r)
            if nxt <= lo:
                nxt = F(lo)
            rates.append(nxt)
            if len(rates) > max_levels:
                raise ValueError("escalera con mas de %d niveles" % max_levels)
            if nxt == lo:
                break
            r = nxt
    return rates


def ladder_record(rates, lam_of, ends=ENDPOINTS):
    """Registro sellable de la escalera y comprobacion certificada de que ningun paso supera la regla:
    sqrt(w_j) ln(r_j/r_{j+1}) <= ln 2 con cotas dirigidas."""
    t = len(rates) - 1
    lams = [F(lam_of(r)) for r in rates]
    steps_ok = True
    for j in range(t):
        w = max(F(W_MIN), lams[j])
        lhs = d_up(NE.multiply(d_up(NE.sqrt(up_of(w))), d_up(NE.ln(up_of(rates[j] / rates[j + 1])))))
        steps_ok &= lhs <= d_dn(LN2)
    idx = {fstr(e): rates.index(F(e)) if F(e) in rates else None for e in ends}
    rec = {"tasas": [fstr(r) for r in rates], "t": t, "niveles": t + 1,
           "lambda": [fstr(x) for x in lams], "lambda_vista": [dstr(x, 12) for x in lams],
           "indices_de_los_extremos": idx, "pasos_no_mayores_que_la_regla": steps_ok,
           "representacion": {"cifras_significativas": RATE_DIGITS, "redondeo": "techo (ROUND_CEILING)",
                              "margen_relativo": str(RATE_MARGIN), "precision_decimal": DEC_PREC,
                              "w_min": W_MIN}}
    rec["sha256"] = canon_sha({k: rec[k] for k in ("tasas", "lambda", "representacion")})
    return rec


# ------------------------------------------------------------------ leyes, soporte y claves (2.3; V3)
def check_laws(laws, ncls, where):
    """laws: por sitio, lista [l(s,0), l(s,1), ..., l(s,n_s)] de fracciones. V3: todas > 0; suman 1."""
    if len(laws) != len(ncls):
        raise ValueError("%s: numero de sitios distinto" % where)
    for s, row in enumerate(laws):
        if len(row) != ncls[s] + 1:
            raise ValueError("%s: sitio %d con %d valores (se esperan %d)" % (where, s, len(row), ncls[s] + 1))
        if any((not isinstance(x, F)) or x <= 0 for x in row):
            raise ValueError("%s: V3 falla en el sitio %d (probabilidad no positiva)" % (where, s))
        if sum(row) != 1:
            raise ValueError("%s: la ley del sitio %d no suma 1" % (where, s))


def lam_from_laws(laws):
    return sum((1 - row[0] for row in laws), F(0))


class KeyTable:
    """Claves de peso: pares (s, c), c != 0, con la misma sucesion exacta de odds l_j(s,c)/l_j(s,0) en los t+1
    niveles. key_of[s][c] (c >= 1) es la clave; a[k][j] el odds comun; pi0[j] = prod_s l_j(s, 0)."""

    def __init__(self, levels_laws):
        self.t1 = len(levels_laws)
        self.NS = len(levels_laws[0])
        groups = {}
        for s in range(self.NS):
            for c in range(1, len(levels_laws[0][s])):
                vec = tuple(L[s][c] / L[s][0] for L in levels_laws)
                groups.setdefault(vec, []).append((s, c))
        self.a = sorted(groups)
        self.nkeys = len(self.a)
        self.key_of = [[None] * len(levels_laws[0][s]) for s in range(self.NS)]
        for k, vec in enumerate(self.a):
            for s, c in groups[vec]:
                self.key_of[s][c] = k
        self.pi0 = []
        for L in levels_laws:
            p = F(1)
            for s in range(self.NS):
                p *= L[s][0]
            self.pi0.append(p)

    def record(self):
        rec = {"claves": [[fstr(x) for x in vec] for vec in self.a], "clave_de_par": self.key_of,
               "numero_de_claves": self.nkeys, "sitios": self.NS, "niveles": self.t1}
        rec["sha256"] = canon_sha(rec)
        return rec

    def log_pi(self, sig, j):
        """ln pi_j(E) exacto en fraccion -> float (solo para pruebas)."""
        p = self.pi0[j]
        for k, n in sig:
            p *= self.a[k][j] ** n
        return ln_frac(p)


def sig_from_classes(cls, key_of):
    cnt = {}
    for s, c in enumerate(cls):
        if c:
            k = key_of[s][c]
            cnt[k] = cnt.get(k, 0) + 1
    return tuple(sorted(cnt.items()))


def sig_weight(sig):
    return sum(n for _, n in sig)


class LogTables:
    """Por cociente j (niveles j y j+1): z_j(E) = ln(pi_j(E)/pi_{j+1}(E)) = c0[j] + sum_k N_k d[j][k]."""

    def __init__(self, kt):
        self.t = kt.t1 - 1
        self.c0, self.d = [], []
        for j in range(self.t):
            self.c0.append(ln_frac(kt.pi0[j] / kt.pi0[j + 1], "c0"))
            self.d.append([ln_frac(kt.a[k][j] / kt.a[k][j + 1], "d") for k in range(kt.nkeys)])

    def z(self, sig, j):
        terms = [self.c0[j]]
        dj = self.d[j]
        for k, n in sig:
            terms.extend([dj[k]] * n)
        return fin(math.fsum(terms), "z")

    def record(self):
        rec = {"c0": [repr(x) for x in self.c0], "d": [[repr(x) for x in row] for row in self.d]}
        rec["sha256"] = canon_sha(rec)
        return rec


# ------------------------------------------------------------------ histogramas (2.11; R5)
def merge(blocks):
    out = {}
    for h in blocks:
        for sig, m in h.items():
            out[sig] = out.get(sig, 0) + m
    return out


def hist_n(h):
    return sum(h.values())


def side_lg(sig, lt, j, side):
    """log g(rho) con c = 1: lado '+' del nivel j usa rho_j+ = exp(z_j); lado '-' usa rho_j- = exp(-z_{j-1})."""
    if side == "+":
        return -softplus(lt.z(sig, j))
    return -softplus(-lt.z(sig, j - 1))


class Values:
    """Valores v = exp(lg - Mx) en (0, 1] de un conjunto de firmas, con Mx el maximo de lg. Cuenta los
    subdesbordamientos (v < 2^-1022): ninguno se lee como 0 sin registrarlo."""

    def __init__(self, sigs, lgfun):
        lg = {s: lgfun(s) for s in sigs}
        if not lg:
            raise ErrorNumerico("conjunto de firmas vacio")
        self.mx = max(lg.values())
        self.v, self.under = {}, 0
        for s, x in lg.items():
            v = math.exp(x - self.mx)
            if v < FLOAT_MIN:
                self.under += 1
            self.v[s] = v
        p = {}
        for s, v in self.v.items():
            if v > 0:
                num, den = v.as_integer_ratio()
                p[s] = (num, den.bit_length() - 1)
        self.E = max(e for _, e in p.values()) if p else 0
        self.iv = {s: (num << (self.E - e)) for s, (num, e) in p.items()}
        for s in self.v:
            self.iv.setdefault(s, 0)

    def sums(self, h):
        """(n, S, Q) enteros: n = sum m, S = sum m V, Q = sum m V^2, con V = v 2^E exacto."""
        n = S = Q = 0
        iv = self.iv
        for sig, m in h.items():
            V = iv[sig]
            n += m
            S += m * V
            Q += m * V * V
        return n, S, Q


def mean_frac(S, n):
    return F(S, n)


def var_frac(n, S, Q):
    """Varianza muestral exacta (n - 1 grados) en unidades V^2."""
    if n < 2:
        raise ErrorNumerico("varianza con menos de 2 muestras")
    return F(n * Q - S * S, n * (n - 1))


# ------------------------------------------------------------------ Bennett con medias (2.9)
def log_mean_g(items, lnC, sign, n):
    """log de (1/n) sum_k m_k g(exp(sign (lnC + z_k))), con items = [(z_k, m_k)]. Log-suma-exp con fsum."""
    if n <= 0 or not items:
        raise ErrorNumerico("media de Bennett sin muestras")
    terms = []
    for z, m in items:
        y = sign * (lnC + z)
        terms.append(math.log(m) - softplus(y))
    mx = max(terms)
    under = 0
    acc = []
    for x in terms:
        e = math.exp(x - mx)
        if e < FLOAT_MIN:
            under += 1
        acc.append(e)
    tot = math.fsum(acc)
    if not tot > 0:
        raise ErrorNumerico("suma nula en la media de Bennett")
    return fin(mx + math.log(tot) - math.log(n), "media de Bennett"), under


def bennett(num_items, n_num, den_items, n_den, iters=BENNETT_ITERS):
    """qhat_j(C) = C [(1/T_j) sum_a g(C rho_j+(E_a))] / [(1/T_{j+1}) sum_b g(rho_{j+1}-(E'_b)/C)], con
    rho_j+ = exp(z_j) en las muestras del nivel j y rho_{j+1}- = exp(-z_j) en las del nivel j+1. C <- 1 y tres
    veces q <- qhat(C), C <- q (FF algoritmo 2, lineas 21 a 25, con medias). Devuelve ln q final."""
    lnC = 0.0
    trace, under = [], 0
    a = b = None
    for _ in range(iters):
        a, ua = log_mean_g(num_items, lnC, 1.0, n_num)
        b, ub = log_mean_g(den_items, lnC, -1.0, n_den)
        under += ua + ub
        lnC = fin(math.fsum([lnC, a, -b]), "ln q")
        trace.append(lnC)
    return {"ln_q": lnC, "iteraciones_ln_C": [repr(x) for x in trace], "ln_media_num": a, "ln_media_den": b,
            "subdesbordamientos": under}


def bennett_items(h, lt, j):
    """[(z_j(E), m)] de un histograma; orden canonico de las firmas."""
    return [(lt.z(sig, j), h[sig]) for sig in sorted(h)]


# ------------------------------------------------------------------ criterio de FF (2.8)
def ff_sides(j, t):
    return (["+"] if j < t else []) + (["-"] if j > 0 else [])


def ff_criterion(blocks, j, t, lt, eps=EPS):
    """Criterio de FF de una cadena (c = 1) con sus bloques (sin E_0). Exacto sobre los valores float de g.
    Devuelve el registro con 'cumple' True/False; un caso no evaluable da 'cumple' False con 'indeterminado'."""
    T = sum(hist_n(h) for h in blocks)
    nb = len(blocks)
    rec = {"T": T, "bloques": nb, "lados": {}}
    if nb % 2 or T % 2:
        raise ValueError("el criterio exige un numero par de bloques iguales")
    full = merge(blocks)
    first = merge(blocks[:nb // 2])
    n1 = hist_n(first)
    if 2 * n1 != T:
        raise ValueError("la primera mitad no tiene ceil(T/2) muestras")
    sig2 = delta = None
    indet, null = [], []
    for side in ff_sides(j, t):
        try:
            vals = Values(full.keys(), lambda s: side_lg(s, lt, j, side))
            if vals.under:
                raise ErrorNumerico("%d valores de g subdesbordados" % vals.under)
            n, S, Q = vals.sums(full)
            _, S1, _ = vals.sums(first)
            if S <= 0 or S1 <= 0:
                raise ErrorNumerico("media de g nula")
            num = n * Q - S * S
            a = F(num, (n - 1) * S * S)                      # sigma_lado^2 = s^2 / (gm^2 T)
            dl = F(abs(S * n1 - S1 * n), S * n1)             # |gm - gm'| / gm
            if num == 0:
                null.append(side)
            rec["lados"][side] = {"sigma": dstr(NE.sqrt(d_of(a)), 20), "Delta": dstr(dl, 20),
                                  "varianza_nula": num == 0}
            sig2 = a if sig2 is None else max(sig2, a)
            delta = dl if delta is None else max(delta, dl)
        except ErrorNumerico as ex:
            indet.append("%s: %s" % (side, ex))
    thr_dn = d_dn(NE.divide(dn_of(eps), d_up(NE.sqrt(d_of(t + 1)))))
    rec["umbral_eps_sobre_raiz_de_t_mas_1"] = dstr(thr_dn, 20)
    if indet or sig2 is None:
        rec.update(cumple=False, indeterminado=indet or ["sin lados"], varianza_nula=null)
        return rec
    lhs_up = d_up(NE.add(d_up(NE.sqrt(up_of(sig2))), up_of(delta)))
    rec.update(sigma=dstr(NE.sqrt(d_of(sig2)), 20), Delta=dstr(delta, 20), sigma_mas_Delta_cota=dstr(lhs_up, 20),
               cumple=bool(lhs_up <= thr_dn), varianza_nula=null)
    return rec


def ff_lengths(t_init=T_INIT, lam0=LAMBDA_INI, max_ext=MAX_EXT):
    """Longitudes sucesivas: T_init, y T <- T + lambda' T_init con lambda' = 2, 4, 8, 16."""
    out, T, lam = [t_init], t_init, lam0
    for _ in range(max_ext):
        T += lam * t_init
        lam *= 2
        out.append(T)
    return out


# ------------------------------------------------------------------ diagnosticos (3.5; C1)
IND = "INDETERMINADO"


def batch_means(blocks, j, t, lt, nbatch=BATCHES):
    """d3: tau = (T/nbatch) var(medias de lote) / s^2 y ESS = T/tau por lado. s^2 = 0 o var = 0 -> INDETERMINADO."""
    nb = len(blocks)
    out = {}
    if nb % nbatch:
        return {s: {"estado": IND, "motivo": "bloques no divisibles en %d lotes" % nbatch} for s in ff_sides(j, t)}
    per = nb // nbatch
    full = merge(blocks)
    batches = [merge(blocks[i * per:(i + 1) * per]) for i in range(nbatch)]
    for side in ff_sides(j, t):
        try:
            vals = Values(full.keys(), lambda s: side_lg(s, lt, j, side))
            if vals.under:
                raise ErrorNumerico("%d valores subdesbordados" % vals.under)
            n, S, Q = vals.sums(full)
            sb = [vals.sums(b) for b in batches]
            nn = sb[0][0]
            if any(x[0] != nn for x in sb):
                raise ErrorNumerico("lotes de tamano distinto")
            s2 = var_frac(n, S, Q)
            Ssum = sum(x[1] for x in sb)
            vb = F(nbatch * sum(x[1] * x[1] for x in sb) - Ssum * Ssum, nbatch * (nbatch - 1) * nn * nn)
            if s2 == 0:
                out[side] = {"estado": IND, "motivo": "s^2 = 0"}
                continue
            if vb == 0:
                out[side] = {"estado": IND, "motivo": "varianza de las medias de lote nula (tau = 0)"}
                continue
            tau = F(nn) * vb / s2
            out[side] = {"estado": "CALCULADO", "tau": dstr(tau, 12), "ESS": dstr(F(n) / tau, 12)}
        except ErrorNumerico as ex:
            out[side] = {"estado": IND, "motivo": str(ex)}
    return out


def rhat(chains_blocks, j, t, lt, side, block=BLOCK):
    """d2: R-hat partido con las cadenas de un nivel (36 en el piloto): n = minimo de floor(T/2); de cada cadena
    sus ultimas 2n muestras en dos series de n."""
    Ts = [sum(hist_n(h) for h in bl) for bl in chains_blocks]
    n = min(T // 2 for T in Ts)
    if n % block or n < 2:
        return {"estado": IND, "motivo": "n = %d no es multiplo de B" % n}
    nbk = n // block
    series = []
    for bl in chains_blocks:
        last = bl[len(bl) - 2 * nbk:]
        series.append(merge(last[:nbk]))
        series.append(merge(last[nbk:]))
    return rhat_series(series, n, j, t, lt, side)


def rhat_series(series, n, j, t, lt, side):
    """R-hat con las series ya formadas (histogramas de n muestras cada una). Reglas C1:
    W = 0 y Bm > 0 -> R-hat infinito y alerta; W = 0 y Bm = 0 -> INDETERMINADO; no finito -> INDETERMINADO."""
    try:
        allsig = set()
        for h in series:
            allsig.update(h)
        vals = Values(allsig, lambda s: side_lg(s, lt, j, side))
        if vals.under:
            raise ErrorNumerico("%d valores subdesbordados" % vals.under)
        st = [vals.sums(h) for h in series]
        if any(x[0] != n for x in st):
            raise ErrorNumerico("series de longitud distinta")
        K = len(st)
        W = sum((var_frac(*x) for x in st), F(0)) / K
        Ssum = sum(x[1] for x in st)
        Bm = F(K * sum(x[1] * x[1] for x in st) - Ssum * Ssum, K * (K - 1) * n)
    except ErrorNumerico as ex:
        return {"estado": IND, "motivo": str(ex), "n": n}
    rec = {"n": n, "series": len(series), "W_es_cero": W == 0, "Bm_es_cero": Bm == 0}
    if W == 0:
        if Bm > 0:
            rec.update(estado="CALCULADO", R_hat="inf", alerta=True, motivo="W = 0 y Bm > 0")
        else:
            rec.update(estado=IND, motivo="W = 0 y Bm = 0")
        return rec
    V = F(n - 1, n) * W + Bm / n
    ratio = V / W
    rec.update(estado="CALCULADO", R_hat=dstr(NE.sqrt(d_of(ratio)), 12), alerta=ratio > RHAT_THR * RHAT_THR)
    return rec


def occupancy_fraction(blocks, k):
    """Fraccion exacta de muestras con |E| <= k en una cadena."""
    tot = c = 0
    for h in blocks:
        for sig, m in h.items():
            tot += m
            if sig_weight(sig) <= k:
                c += m
    return F(c, tot)


def d5_occupancy(fracs_by_seed, ref_lo, ref_hi, z=D5_Z):
    """d5: fracs_by_seed = lista (12) de listas (3) de fracciones exactas. Fraccion por semilla = media de sus
    repeticiones; ohat = media de las semillas; se^2 = var/12. Alerta si ohat queda fuera de
    [ref_lo - z se, ref_hi + z se]. Con se = 0 la banda es la referencia (no hay division: no es degenerado)."""
    f = [sum(x, F(0)) / len(x) for x in fracs_by_seed]
    L = len(f)
    o = sum(f, F(0)) / L
    var = sum(((x - o) ** 2 for x in f), F(0)) / (L - 1)
    se2 = var / L
    below = o < ref_lo and (ref_lo - o) ** 2 > z * z * se2
    above = o > ref_hi and (o - ref_hi) ** 2 > z * z * se2
    return {"ohat": dstr(o, 12), "se": dstr(NE.sqrt(d_of(se2)), 12), "referencia": [dstr(ref_lo, 12),
            dstr(ref_hi, 12)], "se_nulo": se2 == 0, "alerta": bool(below or above), "estado": "CALCULADO"}


def floor25_bracket(view):
    """Vista de 25 cifras con redondeo hacia abajo (t2_rival.show = t2_motor.dec(x, False)) -> [v, v + ulp]."""
    d = Decimal(view)
    if not d.is_finite() or d < 0:
        raise ValueError("vista no valida: %s" % view)
    if d == 0:
        return F(0), F(0)
    if len(d.as_tuple().digits) > 25:
        raise ValueError("la vista tiene mas de 25 cifras: %s" % view)
    ulp = Decimal(1).scaleb(d.adjusted() - 24)
    return F(d), F(d) + F(ulp)


# ------------------------------------------------------------------ Student (2.9; C2)
def _atan(x, ctx):
    """atan(x) para x >= 0 en el contexto ctx: reduccion x <- x/(1+sqrt(1+x^2)) y serie de Taylor."""
    k = 0
    one = Decimal(1)
    while x > Decimal("0.05"):
        x = ctx.divide(x, ctx.add(one, ctx.sqrt(ctx.add(one, ctx.multiply(x, x)))))
        k += 1
    x2 = ctx.multiply(x, x)
    term, s, n = x, x, 1
    eps = Decimal(10) ** (-(ctx.prec + 5))
    while True:
        term = ctx.minus(ctx.multiply(term, x2))
        n += 2
        add = ctx.divide(term, Decimal(n))
        if abs(add) < eps:
            break
        s = ctx.add(s, add)
    return ctx.multiply(s, Decimal(2 ** k))


_WIDE = Context(prec=DEC_PREC + 20, rounding=ROUND_HALF_EVEN, Emin=-_E, Emax=_E)
PI = NE.plus(_WIDE.multiply(Decimal(4), _atan(Decimal(1), _WIDE)))


def student_abs_cdf(t, nu):
    """P(|T_nu| <= t), nu impar, t >= 0 (Decimal), por la forma cerrada de grados impares:
    (2/pi)(theta + sin(theta)(cos + (2/3) cos^3 + ... + (2.4...(nu-3))/(3.5...(nu-2)) cos^(nu-2))),
    theta = atan(t/sqrt(nu)). Sin funciones especiales."""
    if nu < 1 or nu % 2 == 0:
        raise ValueError("nu debe ser impar >= 1")
    c = _WIDE
    t = Decimal(t)
    den = c.add(Decimal(nu), c.multiply(t, t))
    cs = c.sqrt(c.divide(Decimal(nu), den))
    sn = c.divide(t, c.sqrt(den))
    th = _atan(c.divide(t, c.sqrt(Decimal(nu))), c)
    s, coef, pw, c2 = Decimal(0), Decimal(1), cs, c.multiply(cs, cs)
    for k in range((nu - 1) // 2):
        s = c.add(s, c.multiply(coef, pw))
        coef = c.multiply(coef, c.divide(Decimal(2 * k + 2), Decimal(2 * k + 3)))
        pw = c.multiply(pw, c2)
    pi_w = c.multiply(Decimal(4), _atan(Decimal(1), c))
    return NE.plus(c.divide(c.multiply(Decimal(2), c.add(th, c.multiply(sn, s))), pi_w))


def student_quantile(p=STUDENT_P, nu=STUDENT_NU, digits=25):
    """Cuantil t_{nu; p} (p > 1/2) por biseccion sobre la forma cerrada hasta 1e-45 relativo, y techo a
    `digits` cifras tras sumar 1e-40 relativo: el valor devuelto es >= el cuantil exacto (conservador)."""
    p = Decimal(p)
    target = NE.subtract(NE.multiply(Decimal(2), p), Decimal(1))       # P(|T| <= t) = 2p - 1
    lo, hi = Decimal(0), Decimal(10)
    while student_abs_cdf(hi, nu) < target:
        hi = NE.multiply(hi, Decimal(2))
        if hi > Decimal(10) ** 12:
            raise ErrorNumerico("cuantil de Student fuera de rango")
    for _ in range(400):
        if NE.subtract(hi, lo) <= NE.multiply(hi, Decimal("1e-45")):
            break
        mid = NE.divide(NE.add(lo, hi), Decimal(2))
        if student_abs_cdf(mid, nu) < target:
            lo = mid
        else:
            hi = mid
    v = NE.multiply(hi, NE.add(Decimal(1), Decimal("1e-40")))
    return Context(prec=digits, rounding=ROUND_CEILING, Emin=-_E, Emax=_E).plus(v)


# ------------------------------------------------------------------ ancla (2.6; R2; C3)
def load_cp(path, sha):
    """Importa t2_cp.py (T2-1, sin cambios) tras comprobar su SHA-256."""
    with open(path, "rb") as fh:
        got = hashlib.sha256(fh.read()).hexdigest()
    if got != sha:
        raise ValueError("SHA-256 distinto del registrado: " + path)
    spec = importlib.util.spec_from_file_location("t2_cp_t2_2", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def anchor_prefix(blocks, s, Ba):
    """Parada exacta en el exito s sobre el prefijo contiguo de bloques terminados en orden logico.
    blocks: dict b -> {"historias": n, "posiciones": [posiciones 1..n de las silenciosas, crecientes]}.
    Devuelve T (indice logico de la silenciosa numero s) o None si el prefijo aun no tiene s silenciosas."""
    cum, b = 0, 0
    while b in blocks:
        rec = blocks[b]
        pos = rec["posiciones"]
        if rec["historias"] != Ba or any(not 1 <= p <= Ba for p in pos) or any(
                x >= y for x, y in zip(pos, pos[1:])):
            raise ValueError("bloque %d del prefijo incompleto o mal formado" % b)
        if cum + len(pos) >= s:
            T = b * Ba + pos[s - cum - 1]
            return {"T": T, "s": s, "bloque_de_parada": b, "bloques_del_prefijo": b + 1,
                    "silenciosas_del_prefijo_hasta_T": s,
                    "silenciosas_del_bloque_de_parada_tras_T": len(pos) - (s - cum),
                    "historias_del_bloque_de_parada_tras_T": Ba - pos[s - cum - 1]}
        cum += len(pos)
        b += 1
    return {"T": None, "silenciosas_del_prefijo": cum, "bloques_del_prefijo": b}


def anchor_seeds(blocks, nseeds, Ba):
    """Indices logicos y (bloque, posicion) de las silenciosas 1..nseeds en orden logico."""
    out = []
    b = 0
    while len(out) < nseeds and b in blocks:
        for p in blocks[b]["posiciones"]:
            out.append({"indice_logico": b * Ba + p, "bloque": b, "posicion": p})
            if len(out) == nseeds:
                break
        b += 1
    if len(out) < nseeds:
        raise ValueError("menos de %d silenciosas en el prefijo" % nseeds)
    return out


def anchor_interval(s, T, cp, half=HALF_ANCHOR, target=TARGET_CP):
    """Intervalo al 97.5 % por inversion conservadora de la binomial negativa (R2):
    l = cp_lower(s, T) y u = cp_upper(s - 1, T - 1) de t2_cp.py con alfa/2 = 0.0125 en cada cola; u = 1 si
    T = s. Se toma el campo 'valor' (lado conservador). Alcance (C3): describe el procedimiento sin truncacion
    con historias Bernoulli ideales (V6); no es una cobertura condicionada a terminar bajo el techo."""
    if not (isinstance(s, int) and isinstance(T, int)) or isinstance(s, bool) or isinstance(T, bool):
        raise ValueError("s y T deben ser enteros")
    if s < 2 or T < s:
        raise ValueError("se requieren s >= 2 y T >= s")
    if s - 1 > cp.XMAX:
        raise ValueError("s - 1 = %d fuera del alcance de t2_cp (XMAX = %d)" % (s - 1, cp.XMAX))
    h, tg = Decimal(half), Decimal(target)
    if not (h.is_finite() and 0 < h < Decimal("0.5") and tg.is_finite() and 0 < tg < 1):
        raise ValueError("alfa/2 u objetivo no validos")
    lo = cp.cp_lower(s, T, h, tg)
    if T == s:
        up = {"valor": Decimal(1), "raiz_mayor_que": Decimal(1), "pasos": 0, "detenida": False}
    else:
        up = cp.cp_upper(s - 1, T - 1, h, tg)
    l, u = F(lo["valor"]), F(up["valor"])
    if not 0 <= l <= u <= 1:
        raise ErrorNumerico("intervalo del ancla mal ordenado")
    return {"l": l, "u": u, "vista": [str(lo["valor"]), str(up["valor"])],
            "l_raiz_menor_que": str(lo["raiz_menor_que"]), "u_raiz_mayor_que": str(up["raiz_mayor_que"]),
            "pasos": [lo["pasos"], up["pasos"]], "biseccion_detenida": [lo["detenida"], up["detenida"]],
            "alfa_por_cola": half, "objetivo": target, "precision_cifras": cp.PREC, "ulps_ensanche": cp.ULPS}


# ------------------------------------------------------------------ estimacion (2.9, 2.10)
def estimate_rate(lnq_by_rep, jt, anchor_iv, phat0, tq, L=L_SEEDS, M=M_REPS):
    """lnq_by_rep[(l, i)] = [ln q_0, ..., ln q_{t-1}] (l = 1..L, i = 1..M); jt = indice de la tasa objetivo.
    Qhat_i = exp(sum_{j < jt} ln q_j); Qbar_l = media de M; Qbar = media de L; s_Q con L - 1 grados.
    [Q_lo, Q_hi] = Qbar -/+ t s_Q/sqrt(L) (semiancho con techo dirigido); I = [l max(Q_lo, 0), u Q_hi]."""
    keys = sorted(lnq_by_rep)
    if keys != [(l, i) for l in range(1, L + 1) for i in range(1, M + 1)]:
        raise ValueError("faltan repeticiones")
    qh = {}
    for key in keys:
        v = math.exp(fin(math.fsum(lnq_by_rep[key][:jt]), "ln Qhat"))
        if not (math.isfinite(v) and v >= FLOAT_MIN):
            raise ErrorNumerico("Qhat fuera del rango normal")
        qh[key] = F(v)
    qbl = [sum((qh[(l, i)] for i in range(1, M + 1)), F(0)) / M for l in range(1, L + 1)]
    qbar = sum(qbl, F(0)) / L
    s2 = sum(((x - qbar) ** 2 for x in qbl), F(0)) / (L - 1)
    if s2 == 0:
        hw = F(0)                    # exacto: sin ensanche de 2 ulps alrededor de 0 (fracciones gigantes)
    else:
        hw = F(d_up(NE.multiply(d_up(Decimal(tq)), d_up(NE.sqrt(up_of(s2 / L))))))
    qlo, qhi = qbar - hw, qbar + hw
    I = (anchor_iv["l"] * max(qlo, F(0)), anchor_iv["u"] * qhi)
    x = phat0 * qbar
    rec = {"Qbar": qbar, "s_Q2": s2, "Q_lo": qlo, "Q_hi": qhi, "x": x, "I": I, "s_Q_nula": s2 == 0,
           "Qbar_por_semilla": qbl, "Qhat_por_repeticion": {"%d.%d" % k: qh[k] for k in keys}}
    return rec


def precision_measures(x, I):
    lo, hi = I
    out = {}
    if lo > 0:
        out["E_nom"] = max(abs(x - lo) / lo, abs(x - hi) / hi)
    else:
        out["E_nom"] = None
    if lo + hi > 0:
        out["c_h"] = 2 * lo * hi / (lo + hi)
        out["rho_h"] = (hi - lo) / (hi + lo)
    return out


# ------------------------------------------------------------------ decisiones (7.2 a 7.4; C1, C2)
def e_comp_of(J):
    L, U = J
    return (U - L) / (U + L)


def decide_rate(est, J, flags, dmc=None):
    """Decisiones de una tasa objetivo con fracciones exactas. flags: ff_cumple, a1, a2, indeterminados (lista
    de motivos de C1). dmc: [p_inf, p_sup] exactos (solo 1e-3)."""
    x, (Ilo, Ihi) = est["x"], est["I"]
    L, U = J
    pm = precision_measures(x, (Ilo, Ihi))
    D = {}
    Enom = pm["E_nom"]
    D["D_s1"] = "PRECISA_NOMINAL" if (Ilo > 0 and Enom <= PREC_THR) else "IMPRECISA_NOMINAL"
    D["D_s2"] = "CRITERIO_FF_CUMPLE" if flags["ff_cumple"] else "CRITERIO_FF_NO_CUMPLE"
    if Ilo <= L and U <= Ihi:
        D["D_s3"] = "CUBRE_EN_ESTA_REALIZACION"
    elif Ihi < L or U < Ilo:
        D["D_s3"] = "DISJUNTO"
    else:
        D["D_s3"] = "COBERTURA_INDETERMINADA"
    D["compatible"] = not (Ihi < L or U < Ilo)
    Emax = Emin = None
    if x > 0 and L > 0:
        Emax = max(abs(x - L) / L, abs(x - U) / U)
        Emin = F(0) if L <= x <= U else min(abs(x - L) / L, abs(x - U) / U)
        D["D_s4"] = ("PRECISION_REAL_ACREDITADA" if Emax <= PREC_THR else
                     "PRECISION_REAL_DESCARTADA" if Emin > PREC_THR else "PRECISION_REAL_INDETERMINADA")
    else:
        D["D_s4"] = "PRECISION_REAL_INDETERMINADA"
    if flags["a1"] or flags["a2"]:
        D["D_s5"] = "ALERTA_MEZCLA"
    elif flags["indeterminados"]:
        D["D_s5"] = "DIAGNOSTICO_INDETERMINADO"
    else:
        D["D_s5"] = "SIN_ALERTA"
    if dmc is not None:
        D["D_s6"] = "COMPATIBLE_DMC" if (Ilo <= dmc[1] and dmc[0] <= Ihi) else "DISJUNTO_DMC"
    if Emax is not None and Enom is not None:
        D["D_s7"] = ("REAL_NO_SUPERA_DECLARADO" if Emax <= Enom else
                     "REAL_SUPERA_DECLARADO" if Emin > Enom else "INDETERMINADO")
    else:
        D["D_s7"] = "INDETERMINADO"
    D["estimacion_solo_diagnostico"] = not flags["ff_cumple"]
    D["complementa_nominal"] = (e_comp_of(J) > PREC_THR and D["D_s1"] == "PRECISA_NOMINAL"
                                and D["D_s2"] == "CRITERIO_FF_CUMPLE" and D["D_s5"] == "SIN_ALERTA"
                                and D["D_s3"] != "DISJUNTO")
    num = {"E_nom": Enom, "E_max": Emax, "E_min": Emin, "c_h": pm.get("c_h"), "rho_h": pm.get("rho_h"),
           "e_comp_4": e_comp_of(J)}
    return D, num


def verdict(per_rate, incomplete_reason=None):
    """Veredicto operativo (7.3): el primero que se cumpla. per_rate: dict tasa -> D (o None si no evaluada)."""
    if incomplete_reason or any(v is None for v in per_rate.values()):
        return "PILOTO_INCOMPLETO"
    if any(v["D_s3"] == "DISJUNTO" for v in per_rate.values()):
        return "PILOTO_CON_DISJUNCION"
    if any(v["D_s2"] == "CRITERIO_FF_NO_CUMPLE" or v["D_s5"] == "ALERTA_MEZCLA" for v in per_rate.values()):
        return "PILOTO_CON_ALERTA"
    if any(v["D_s5"] == "DIAGNOSTICO_INDETERMINADO" for v in per_rate.values()):
        return "PILOTO_CON_DIAGNOSTICO_INDETERMINADO"
    return "PILOTO_OPERATIVO_SIN_ALERTAS"


def extension_rule(ver, per_rate):
    hi = [fstr(TARGETS[0]), fstr(TARGETS[1])]
    comp = [k for k, v in per_rate.items() if v is not None and v["complementa_nominal"]]
    if ver == "PILOTO_OPERATIVO_SIN_ALERTAS" and any(k in hi for k in comp):
        rec = "RECOMENDAR_AMPLIACION"
        why = "veredicto operativo y COMPLEMENTA_NOMINAL en " + ", ".join(k for k in comp if k in hi)
    else:
        rec = "NO_RECOMENDAR_AMPLIACION"
        why = ("veredicto %s" % ver if ver != "PILOTO_OPERATIVO_SIN_ALERTAS" else
               "sin COMPLEMENTA_NOMINAL en 1e-3 ni en 1e-4")
    return {"regla": rec, "motivo": why, "complementa_nominal_en": comp,
            "lectura_global": "SIN_VENTAJA_ACREDITADA" if not comp else "COMPLEMENTA_NOMINAL",
            "VENTAJA_TEMPORAL": "NO_ACREDITADA"}
