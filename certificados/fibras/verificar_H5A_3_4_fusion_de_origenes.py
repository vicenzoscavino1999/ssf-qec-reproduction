"""HGP (3,4). Fusiones de orígenes distintos con varias llegadas del mismo origen."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import permutations, product
from math import ceil, prod
from pathlib import Path
from runpy import run_path

HERE = Path(__file__).resolve().parent
SHORT = run_path(str(HERE / "verificar_H5A_3_4_reuniones_cortas.py"))
BLOCK = SHORT["BLOCK"]
Z = SHORT["Z"]
CELLS = SHORT["cells"]
TERMINALS = SHORT["TERMINALS"]
COL_PERMS = tuple(permutations(range(3)))
UNIT = Q(1, 10**6)


def dot(left, right):
    return sum(x*y for x, y in zip(left, right))


def apply_matrix(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def nonempty_paths(matrix):
    return [[x-int(i == j) for j, x in enumerate(row)] for i, row in enumerate(matrix)]


def round_up(value):
    return Q(ceil(value * 10**12), 10**12)


def check_mass(vector):
    return 12 * dot((1, 3, 6), vector)


@lru_cache(maxsize=None)
def canonical(mask):
    """Tipo bajo permutaciones independientes de las cuatro filas y las tres columnas."""
    images = []
    for perm in COL_PERMS:
        rows = sorted(sum(1 << perm[c] for c in range(3) if mask & (1 << (3*r+c)))
                      for r in range(4))
        images.append(sum(row << (3*i) for i, row in enumerate(rows)))
    return min(images)


def bipartitions(mask):
    """Particiones no ordenadas en dos bloques no vacíos, con el menor bit a la izquierda."""
    rest = mask ^ (mask & -mask)
    sub = rest
    while sub:
        yield mask ^ sub, sub
        sub = (sub-1) & rest


@lru_cache(maxsize=None)
def partitions(mask):
    """(número de bits frescos, tipos de bloques heredados), con su multiplicidad."""
    if not mask:
        return Counter({(0, ()): 1})
    bit = mask & -mask
    rest = mask ^ bit
    out = Counter()
    for (fresh, groups), count in partitions(rest).items():
        out[fresh+1, groups] += count
    sub = rest
    while True:
        group = canonical(bit | sub)
        for (fresh, groups), count in partitions(rest ^ sub).items():
            out[fresh, tuple(sorted((*groups, group)))] += count
        if not sub:
            return out
        sub = (sub-1) & rest


@lru_cache(maxsize=None)
def incidence(mask):
    """Estados representados cuya salida contiene todas las marcas."""
    if mask.bit_count() == 1:
        return (12, 36, 72)
    if mask.bit_count() == 2:
        a, b = CELLS(mask)
        if a//3 != b//3:
            return (0, 12, 0) if a%3 == b%3 else (0, 0, 12)
    return (0, 0, 0)


def direct_origins(mask, hc, hv):
    kind = incidence(mask)
    if kind == (12, 36, 72):
        return 21*9*hc + 10*28*hv
    if kind == (0, 12, 0):
        return 15*9*hc + 8*28*hv
    if kind == (0, 0, 12):
        return 15*6*hc + 8*16*hv
    return 0


class OpenBranches:
    """Salidas de un mismo tronco que llegan por varios brazos a un lector fijado."""

    def __init__(self, forward, forward_star, reverse_star, endpoint):
        self.one = forward[0]
        self.endpoint = endpoint
        remainder = nonempty_paths(forward_star)
        self.long = tuple(map(round_up, BLOCK["multiply"](forward, remainder)[0]))
        self.split = round_up(sum(remainder[0][1:]))
        self.reverse = nonempty_paths(reverse_star)

    @lru_cache(maxsize=None)
    def extra_forward(self, marks):
        total = 0
        for left, right in bipartitions(marks):
            a = self.maximum(canonical(left))
            b = self.maximum(canonical(right))
            # Dos brazos directos ya forman una sola emisión.
            total += a * (right.bit_count() == 1) + b * (left.bit_count() == 1) + a*b
        return round_up(2*self.split*total)

    @lru_cache(maxsize=None)
    def forward(self, source, marks):
        types = [i for i, count in enumerate(incidence(marks)) if count]
        closed = 0
        if types:
            once = min(self.endpoint(source, CELLS(marks)), sum(self.one[i] for i in types))
            closed = once + sum(self.long[i] for i in types)
        return round_up(closed + self.extra_forward(canonical(marks)))

    @lru_cache(maxsize=None)
    def maximum(self, marks):
        # 4 y 3 representan una coordenada exterior en cada eje.
        return max(self.forward(site, marks) for site in product(range(5), range(4)))

    @lru_cache(maxsize=None)
    def open_fork(self, marks):
        """Masa inversa hasta la salida bifurcante, excluida su producción."""
        bound = [0, Q(0), Q(0)]
        if marks.bit_count() <= 1:
            return tuple(bound)
        for x in CELLS(marks):
            row, col = divmod(x, 3)
            rest = marks ^ (1 << x)
            row_sets = SHORT["other_axis"](4, row)
            col_sets = SHORT["other_axis"](3, col)
            bound[1] += 12 * max(sum(self.forward((r, col), rest) for r in rr) for rr in row_sets)
            bound[2] += 12 * max(sum(self.forward((r, c), rest) for r, c in product(rr, cc))
                                 for rr, cc in product(row_sets, col_sets))
        for left, right in bipartitions(marks):
            a, b = canonical(left), canonical(right)
            back_a, back_b = self.positive(a), self.positive(b)
            for i in (1, 2):
                bound[i] += 2*min(back_a[i]*self.maximum(b), back_b[i]*self.maximum(a))
        return tuple(map(round_up, bound))

    @lru_cache(maxsize=None)
    def positive(self, marks):
        """Rama con al menos un paso nuevo antes de su terminal o primera bifurcación abierta."""
        heads = tuple(x+y for x, y in zip(incidence(marks), self.open_fork(marks)))
        return tuple(round_up(dot(heads, col)) for col in zip(*self.reverse))


def reader_polynomials(winners, rows, cols):
    result = []
    # Número de particiones del conjunto de luces y un elemento que distingue el conjunto fresco.
    bell = (1, 2, 5, 15, 52, 203, 877)
    for terminal in TERMINALS:
        terms = Counter()
        for lit, foot in winners:
            if lit ^ foot != terminal:
                continue
            mult = 3 if foot in rows else 4 if foot in cols else 1
            split = partitions(lit)
            assert sum(split.values()) == bell[lit.bit_count()]
            for (fresh, groups), count in split.items():
                if groups:
                    terms[groups] += Q(count, mult) * Z**fresh
        result.append(terms)
    return result


class FusionBound:
    def __init__(self, winners, rows, cols):
        prior = run_path(str(HERE / "verificar_H5A_3_4_reunion_de_ramas.py"))
        forward = [row[1:] for row in prior["causal_rows"](Z)[0]]
        reverse = BLOCK["reverse_rows"](winners, rows, cols)
        relaxed = [BLOCK["reverse_rows"](winners, rows, cols, i) for i in (1, 2)]
        k01, k11 = SHORT["short_rows"](winners, rows, cols)
        self.paths = SHORT["nested_paths"](reverse, k01, k11, relaxed, forward, prior, Q(1))
        arms = tuple(map(Q, (".009746", ".000000623", ".0003633")))
        kernel = [row[1:] for row in prior["fork_rows"](Z, Q(1), arms)]
        all_forward = BLOCK["path_sum"](kernel)
        remainder = nonempty_paths(all_forward)
        assert all(x <= cap for x, cap in zip(remainder[0], arms))
        long_arm = BLOCK["multiply"](forward, remainder)[0]
        self.closed_fork = SHORT["add_long_arms"](SHORT["add"](k01, k11), relaxed, long_arm, self.paths)
        endpoint, _ = SHORT["endpoint_prices"](winners, rows, cols)
        self.branches = OpenBranches(forward, all_forward, self.paths, endpoint)
        self.rows = reader_polynomials(winners, rows, cols)
        self.multi = [{js: w for js, w in row.items() if len(js) >= 2} for row in self.rows]
        self.marks = {j for row in self.multi for groups in row for j in groups}
        assert len(self.marks) == 21 and tuple(map(len, self.multi)) == (123, 14, 39)

    def single(self, hc, hv):
        source = 27*hc + 40*hv
        forks = apply_matrix(self.closed_fork, (source,)*3)
        first = [v+sum(co*direct_origins(js[0], hc, hv) for js, co in row.items() if len(js) == 1)
                 for v, row in zip(forks, self.rows)]
        return apply_matrix(self.paths, first)

    def merged(self, values, hc, hv, arity=None):
        free_source = 27*hc + 40*hv
        fork_source = tuple(free_source+x for x in values)
        prices = {j: direct_origins(j, hc, hv) + dot(incidence(j), values)
                     + dot(self.branches.open_fork(j), fork_source) for j in self.marks}
        joined = [sum(co*prod(prices[j] for j in groups) for groups, co in row.items()
                      if arity is None or len(groups) == arity) for row in self.multi]
        return apply_matrix(self.paths, joined)

    def exact_pair_coefficient(self, hc, hv):
        return check_mass(self.merged(self.single(hc, hv), hc, hv, arity=2))

    def verify(self, hc, hv, coordinates):
        trial = tuple(Q(x)*UNIT for x in coordinates)
        leaves = self.single(hc, hv)
        unions = self.merged(trial, hc, hv)
        assert all(a+b < bound for a, b, bound in zip(leaves, unions, trial))
        return check_mass(trial)


def contact_exponents(hc, hv, mass, tilt, covers):
    # Campos escalados simultáneamente por rho entre 0 y 1.
    query, cc, vv = covers
    lam = Q(tilt)
    rate = Z*lam*(1 + lam/(2*(1-lam/3)))
    assert 0 < lam < 3
    bq = query[2] + query[3]
    cost = (2*cc[2]+cc[3])*hc + (2*vv[2]+vv[3])*hv + (bq+query[2]-1)*mass
    assert lam-bq*rate >= cost
    out = []
    for v, c in ((1, 3), (2, 2)):
        def base(cover):
            return v*cover[0] + c*cover[1] + (2*v+3*c)*cover[3]
        out.append(base(cc)*hc + base(vv)*hv + (mass+rate)*base(query))
    return out


def main():
    if not __debug__:
        raise RuntimeError("Activar los assert de Python para verificar las desigualdades.")
    winners, rows, cols = run_path(str(HERE / "verificar_H5A_3_4_cadenas_un_predecesor.py"))["catalogo"]()
    SHORT["check_symmetry"](winners)
    assert max(h.bit_count() for h, f in winners if h != f) == 6
    model = FusionBound(winners, rows, cols)
    single_cc = check_mass(model.single(1, 0))
    single_vv = check_mass(model.single(0, 1))
    assert Q("520.6014") < single_cc < Q("520.6015")
    assert Q("774.3965") < single_vv < Q("774.3966")
    cc = model.exact_pair_coefficient(1, 0)
    vv = model.exact_pair_coefficient(0, 1)
    mixed = model.exact_pair_coefficient(1, 1) - cc - vv
    pair_caps = (131065526, 392550876, 293958731)
    assert all(cap-1 < value < cap for value, cap in zip((cc, mixed, vv), pair_caps))
    print("Un origen, coeficientes CC y VV:", float(single_cc), float(single_vv))
    print("Dos orígenes, techos CC+CC, CC+VV, VV+VV:", ceil(cc), ceil(mixed), ceil(vv))
    mono = model.verify(UNIT, 0, ("4.7793", "1.3620", "8.1568"))
    double_field = model.verify(2*UNIT, 0, ("15.630", "4.636", "26.365"))
    both = model.verify(UNIT/2, UNIT/2, ("6.514", "1.890", "11.082"))
    assert mono == Q(".0006936732") and double_field == Q(".002252736")
    assert both == Q(".000944112")
    tail = mono - single_cc*UNIT - cc*UNIT**2
    assert 0 < tail < Q(".000042007")
    print("Por check, h_CC=10^-6: total <=", float(mono), "y desde tres orígenes <", 0.000042007)
    print("A campo CC doble:", float(double_field), "y en la esquina mixta (5e-7,5e-7):", float(both))
    covers = tuple(BLOCK["check_cover"](winners, rows, cols, i, j) for i, j in ((2, 2), (5, 4), (4, 5)))
    assert covers == ((210, 252, 90, 114), (16641, 19665, 7353, 8748), (24832, 29584, 11008, 13104))
    for hc, hv, value, tilt, bounds in (
        (UNIT, 0, mono, ".319", ("2.713309", "2.522949")),
        (UNIT/2, UNIT/2, both, ".45", ("3.826855", "3.558220")),
    ):
        exponents = contact_exponents(hc, hv, value, tilt, covers)
        assert all(Q(b)-UNIT < x < Q(b) for x, b in zip(exponents, bounds))
        print("Exponentes hasta el contacto (1V+3C, 2V+2C):", *map(float, exponents))
    print("Los campos ponderan datos marcados. El resto global de G3 se evalúa por separado.")


if __name__ == "__main__":
    main()
