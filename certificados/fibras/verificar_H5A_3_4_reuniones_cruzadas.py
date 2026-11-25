"""HGP (3,4). Dos grupos pueden alimentar repetidamente dos reuniones simultáneas."""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from math import factorial
from pathlib import Path
from runpy import run_path

HERE = Path(__file__).resolve().parent
FUSION = run_path(str(HERE / "verificar_H5A_3_4_fusion_de_origenes.py"))
SHORT, BLOCK = FUSION["SHORT"], FUSION["BLOCK"]
Z, UNIT = FUSION["Z"], FUSION["UNIT"]
TERMINALS, cells = FUSION["TERMINALS"], FUSION["CELLS"]
dot, matvec = FUSION["dot"], FUSION["apply_matrix"]
mass, up = FUSION["check_mass"], FUSION["round_up"]
multiply, add = BLOCK["multiply"], SHORT["add"]
TYPES = tuple(combinations_with_replacement(range(3), 2))
MONOMIALS = tuple(combinations_with_replacement(range(4), 2))


def tighten_paths(model, winners, rows, cols):
    """La primera ocurrencia del brazo principal se paga desde una sola celda física."""
    reverse = BLOCK["reverse_rows"](winners, rows, cols)
    relaxed = [BLOCK["reverse_rows"](winners, rows, cols, j) for j in (1, 2)]
    k01, k11 = SHORT["short_rows"](winners, rows, cols)
    near = add(k01, k11)
    first = tuple(row[0]/12 for row in reverse)
    side = (model.branches.long[0], sum(model.branches.long[1:]))

    def long_forks(paths):
        positive = matvec(paths, first)
        result = []
        for i in range(3):
            value = sum(weight*(matrix[i][0]/12 + dot(matrix[i], positive))
                        for weight, matrix in zip(side, relaxed))
            result.append((0, 36*value, 72*value))
        return result

    paths = [[Q(i == j) for j in range(3)] for i in range(3)]
    for _ in range(100):
        forks = add(near, long_forks(paths))
        mixed = multiply(add(reverse, forks), paths)
        updated = [[up(x+int(i == j)) for j, x in enumerate(row)] for i, row in enumerate(mixed)]
        assert all(x <= y for old, new in zip(paths, updated) for x, y in zip(old, new))
        if updated == paths:
            assert all(x+int(i == j) <= paths[i][j] for i, row in enumerate(mixed)
                       for j, x in enumerate(row))
            model.paths, model.closed_fork = paths, forks
            return first
        paths = updated
    raise ArithmeticError("No se obtuvo una supersolución para los tramos de dos terminales.")


@lru_cache(maxsize=None)
def axes(size, origin):
    """Coincidencias de los demás vecinos con un eje del lector de llegada."""
    inside = [x for x in range(size) if x != origin]
    return [choices+(size,)*(size-1-k) for k in range(size) for choices in combinations(inside, k)]


def unordered_product(left, right, indices):
    return tuple(left[i]*right[j] + (left[j]*right[i] if i != j else 0) for i, j in indices)


class Crossings:
    def __init__(self, model, first, winners, rows, cols):
        self.paths = model.paths
        self.arms = model.branches   # Mayorantes de abanicos ya verificados en v26.
        self.first = first
        self.pairs = [{js: w for js, w in row.items() if len(js) == 2} for row in model.multi]
        seeds, repeat = [], []
        for i, j, equal in [(*pair, False) for pair in TYPES] + [(0, 0, True)]:
            origins, earlier = self.layer(i, j, equal)
            seeds.append(origins)
            repeat.append([*earlier, Q(0)])
        self.repeat = repeat
        shape = tuple(map(Q, ("1", ".00017461", ".027404", ".00000003162", ".00000463",
                              ".00073873", "1")))
        self.factor = max(x/y for x, y in zip(matvec(repeat, shape), shape))
        assert self.factor < Q(".03037")
        star = BLOCK["path_sum"](repeat)
        assert all(x >= 0 for row in star for x in row)
        closing = self.closing_rows(winners, rows, cols)
        # Coeficientes de h_C, h_V, W_columna y W_diagonal, por grado dos.
        full = multiply(self.paths, multiply(closing, multiply(star, seeds)))
        self.coefficients = [tuple(map(up, row)) for row in full]

    @lru_cache(maxsize=None)
    def delivery(self, marks):
        front = [x+y for x, y in zip(FUSION["incidence"](marks), self.arms.open_fork(marks))]
        return tuple(up(dot(front, col)) for col in zip(*self.paths))

    @lru_cache(maxsize=None)
    def physical_start(self, marks):
        return up(dot(self.delivery(marks), self.first))

    @lru_cache(maxsize=None)
    def direct_partner(self, marks, same_row):
        values = [Q(0), Q(0)]
        if same_row:
            assert marks.bit_count() == 1
        mark_row = cells(marks)[0]//3
        for source_r, source_c in product(range(5), range(4)):
            if same_row and source_r != mark_row:
                continue
            rr, cc = axes(4, source_r), axes(3, source_c)
            vertical = max(sum(self.arms.forward((i, source_c), marks) for i in a) for a in rr)
            diagonal = max(sum(self.arms.forward(site, marks) for site in product(a, b))
                           for a, b in product(rr, cc))
            values[0] = max(values[0], 12*vertical)
            values[1] = max(values[1], 12*diagonal)
        return tuple(map(up, values))

    @lru_cache(maxsize=None)
    def fork_prices(self, a, b, same_row):
        """Un solo bifurcador abastece a y b en dos lectores, cada cual por su propia celda."""
        direct = a.bit_count() == b.bit_count() == 1 and not same_row
        both = min(self.physical_start(a)*self.arms.maximum(b),
                   self.physical_start(b)*self.arms.maximum(a))
        branches = []
        for i, through in enumerate((36, 72)):
            price = through*both
            if a.bit_count() == 1:
                price += self.direct_partner(b, same_row)[i]
            if b.bit_count() == 1:
                price += self.direct_partner(a, same_row)[i]
            branches.append(up(price))
        # Con ambos brazos directos el dato se deduplica por sus posiciones posibles.
        source = (27*sum(branches), 40*sum(branches), *branches)
        direct_source = (135, 224, 12, 12) if direct else (0, 0, 0, 0)
        origins = tuple(x+y for x, y in zip(source, direct_source))
        # En una capa repetida se remonta hasta los dos estados anteriores.
        forks = tuple(x+12*int(direct) for x in branches)
        earlier = tuple(up(dot(forks, col)) for col in zip(*self.paths[1:]))
        return origins, earlier

    def layer(self, i, j, same_row):
        costs = [Q(0)]*len(MONOMIALS)
        repeat = [Q(0)]*len(TYPES)
        singleton = FUSION["canonical"](1)
        for (a, b), left in self.pairs[i].items():
            for (c, d), right in self.pairs[j].items():
                price = left*right
                terms = [(False, price)]
                if same_row and a == b == c == d == singleton:
                    # Dos huellas horizontales, cada una de precio 1/3, en la misma fila.
                    assert price >= Q(1, 9)
                    terms = [(False, price-Q(1, 9)), (True, Q(1, 9))]
                for equal, weight in terms:
                    for p, q in ((c, d), (d, c)):
                        f1, r1 = self.fork_prices(a, p, equal)
                        f2, r2 = self.fork_prices(b, q, equal)
                        for n, value in enumerate(unordered_product(f1, f2, MONOMIALS)):
                            costs[n] += weight*value
                        for n, value in enumerate(unordered_product(r1, r2, TYPES)):
                            repeat[n] += weight*value
        return list(map(up, costs)), list(map(up, repeat))

    def closing_rows(self, winners, rows, cols):
        result = []
        for i, terms in enumerate(self.pairs):
            totals = [Q(0)]*len(TYPES)
            for (a, b), weight in terms.items():
                choices = unordered_product(self.delivery(a), self.delivery(b), TYPES)
                for j, value in enumerate(choices):
                    totals[j] += weight*value
            same = Q(0)
            for lit, foot in winners:
                if foot ^ lit != TERMINALS[i]:
                    continue
                pairs = sum(a//3 == b//3 for a, b in combinations(cells(lit), 2))
                mult = 3 if foot in rows else 4 if foot in cols else 1
                same += Q(144*pairs, mult)*Z**(lit.bit_count()-2)
            totals[0] -= same
            assert min(totals) >= 0
            result.append(tuple(map(up, (*totals, same))))
        return result

    def value(self, state, hc, hv):
        inputs = (hc, hv, state[1], state[2])
        return matvec(self.coefficients, tuple(inputs[i]*inputs[j] for i, j in MONOMIALS))



def disjoint_guard(widths, geometry):
    """Débito mínimo por primeras extinciones en franjas disjuntas de paridad uno."""
    @lru_cache(maxsize=None)
    def debit(ranks):
        if not any(ranks):
            return 0
        choices = []
        for lit, foot, in_row, in_col, increment in geometry:
            if any(width*lit < rank*foot for width, rank in zip(widths, ranks)):
                continue
            cost = 138*(2*lit-foot)-increment
            assert cost >= 0
            allowed = [range(min(r, in_row if w == 3 else in_col)+1) for w, r in zip(widths, ranks)]
            for counts in product(*allowed):
                if 0 < sum(counts) <= lit:
                    remaining = tuple(r-a for r, a in zip(ranks, counts))
                    choices.append(cost+debit(remaining))
        return min(choices)
    debits = {r: debit(r) for r in product(*(range(w+1) for w in widths))}
    worst = max([0]+[138*(2*sum(r)-sum(widths))-cost for r, cost in debits.items()])
    return debits, worst


def coverage_discounts(verifier, geometry, old):
    guard = verifier["field_coefficients"]()[0]
    def upper_factor(saving):
        return 1/sum((guard*saving)**j/factorial(j) for j in range(9))
    caps = {}
    for width, expected in ((3, ([0, 0, 72, 414], 66)), (4, ([0, 0, 0, 222, 552], 54))):
        debits, cap = disjoint_guard((width,), geometry)
        assert ([debits[r,] for r in range(width+1)], cap) == expected
        caps[width] = cap
    out = {}
    ordinary_savings = {}
    for (v, c), inherited in old.items():
        n = v+c
        cap = caps[3+c] if n == 1 else verifier["FAR"]["group_budget"](v, c, geometry, n, guard=True)
        saving = max(0, 246*v+340*c-cap)
        ordinary_savings[v, c] = saving
        out[v, c] = upper_factor(saving)
        assert 0 < out[v, c] <= inherited
    pair_factors = {}
    for (i, j), expected, extra in (((0, 0), 342, 138), ((0, 1), 342, 126), ((1, 1), 552, 54)):
        widths = (3+i, 3+j)
        _, cap = disjoint_guard(widths, geometry)
        assert max(caps[w] for w in widths) <= cap == expected
        v, c = 2-i-j, i+j
        assert 246*v+340*c-cap-ordinary_savings[v, c] == extra
        pair_factors[i, j] = upper_factor(extra)
    return out, pair_factors


def overlap_floor(verifier, old_case):
    _, u, h, physical, visible, root = old_case[-1]
    g = 40*Q(physical[0])+27*Q(physical[1])+Q(root)
    fields = (2*Q(u[0]), 2*Q(u[1]), 2*g, Q(h[0]), Q(h[1]), Q(visible))
    return (*verifier["field_coefficients"](), *(Q(19, 25)*x for x in fields))


GLOBAL_CASES = (
    ("Bernoulli", "5.26e-7", None, ".567", ".685",
     (("9.914079608117197e-7", "1.156118735194398e-6"),
      ("6.956966700197279e-7", "7.460821132121261e-7"),
      ("3.864716960986921e-8", "1.37106438383504e-7"),
      ("4.73952176187161e-9", "1.920877433061565e-8"),
      "5.246059226595564e-7", "1.792172590084889e-6")),
    ("VV primero", "9.90e-8", 4, ".558", ".750",
     (("1.416150337583454e-6", "1.622832009421659e-6"),
      ("4.570443026847518e-7", "8.573869079158714e-7"),
      ("5.186001912185859e-8", "2.681534479620213e-7"),
      ("6.791381167922173e-9", "3.764134809756153e-8"),
      "7.158887035545343e-7", "1.455817082256635e-6")),
    ("Cualquier orden de columnas", "9.02e-8", 7, ".555", ".780",
     (("1.291577552031522e-6", "1.480504385152409e-6"),
      ("6.56328432804184e-7", "7.884804388764516e-7"),
      ("8.42896694550738e-8", "2.930816650779389e-7"),
      ("9.106698319377819e-9", "4.054389127367144e-8"),
      "8.294006460971507e-7", "1.329862191796146e-6")),
)


def check_global():
    verifier = run_path(str(HERE / "verificar_H5A_3_4_cola_por_lectores.py"))
    assert verifier["F5"]["RM"] == Q(".001")
    verifier["check_local"]()
    data, circuit = verifier["v21"]["previous_tables"]()
    assert data[2] == circuit[3]
    discounts, pair = coverage_discounts(verifier, circuit[4], data[2])
    data = (*data[:2], discounts, data[3])
    circuit = (*circuit[:3], discounts, *circuit[4:])
    for case, prior in zip(GLOBAL_CASES, verifier["CASES"]):
        assert case[2] == prior[2] and Q(case[1]) > Q(prior[1])
        verifier["verify_model"](case, data, circuit, floor_fields=overlap_floor(verifier, prior),
                                 disjoint_pair_factors=pair)
    floor = overlap_floor(verifier, verifier["CASES"][0])
    for a, s, lower in ((".567", ".685", "1.0009189"), (".569", ".673", "1.0181109")):
        assert verifier["obstruction"](data, Q(a), Q(s), floor_fields=floor) > Q(lower)


def main():
    if not __debug__:
        raise RuntimeError("El verificador requiere los assert de Python.")
    assert Z == Q(".0012") and UNIT == Q("1e-6")
    catalog = run_path(str(HERE / "verificar_H5A_3_4_cadenas_un_predecesor.py"))["catalogo"]
    winners, rows, cols = catalog()
    SHORT["check_symmetry"](winners)
    assert all(foot in rows for lit, foot in winners if lit.bit_count() == 2 and foot != lit)
    assert all(not any((foot^lit) & row == (foot^lit) for row in rows)
               for lit, foot in winners if (foot^lit).bit_count() == 2)
    model = FUSION["FusionBound"](winners, rows, cols)
    previous = mass(model.single(1, 0))
    first = tighten_paths(model, winners, rows, cols)
    cross = Crossings(model, first, winners, rows, cols)
    print("Factor por cada capa sucesiva, en norma ponderada:", float(cross.factor))
    mono = mass(model.single(1, 0))
    assert Q("520.60") < previous < Q("520.61") and Q("161.9148") < mono < Q("161.9149")
    print("Coeficiente lineal CC por check:", float(previous), "->", float(mono))
    two = model.exact_pair_coefficient(1, 0)
    two_cross = mass(cross.value(model.single(1, 0), 1, 0))
    assert 56433752 < two < 56433753 and 23383761 < two_cross < 23383762
    print("Grado dos CC, fusión simple y bloque de cruces:", float(two), float(two_cross))

    covers = tuple(BLOCK["check_cover"](winners, rows, cols, i, j) for i, j in ((2, 2), (5, 4), (4, 5)))
    cases = (
        (UNIT, 0, ("1.8650", ".2570", "3.1803"), ".0002606136", ".136", ("1.138537", "1.058837")),
        (UNIT/2, UNIT/2, ("2.6975", ".4288", "4.4934"), ".0003713316", ".193", ("1.606332", "1.493734")),
        (2*UNIT, 0, ("6.4345", "1.2706", "10.3895"), ".0008709996", None, None),
        (UNIT, UNIT, ("12.0700", "2.5693", "18.9740"), ".0016034628", None, None),
    )
    for hc, hv, coordinates, bound, tilt, limits in cases:
        state = tuple(Q(x)*UNIT for x in coordinates)
        leaf, join = model.single(hc, hv), model.merged(state, hc, hv)
        block = cross.value(state, hc, hv)
        assert all(a+b+c < limit for a, b, c, limit in zip(leaf, join, block, state))
        assert mass(state) == Q(bound)
        print(f"h_C={float(hc):g}, h_V={float(hv):g}. Masa por check <= {bound}.")
        if tilt:
            actual = FUSION["contact_exponents"](hc, hv, mass(state), tilt, covers)
            assert all(Q(bound)-UNIT < value < Q(bound) for value, bound in zip(actual, limits))
            print("  Exponentes para 1V+3C y 2V+2C:", *map(float, actual))
    tail = Q(cases[0][3]) - mono*UNIT - (two+two_cross)*UNIT**2
    assert 0 < tail < Q(".000018882")
    print("En h_C=10^-6, h_V=0, desde tres orígenes <", 0.000018882)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--global", dest="include_global", action="store_true",
                        help="reconstruye las tablas y verifica las tres filas físicas completas")
    args = parser.parse_args()
    main()
    if args.include_global:
        check_global()
