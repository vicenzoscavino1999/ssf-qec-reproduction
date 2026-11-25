"""Reuniones (0,1) y (1,1) con celdas de llegada en un mismo lector de HGP (3,4)."""

from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import ceil
from pathlib import Path
from runpy import run_path


HERE = Path(__file__).resolve().parent
BLOCK = run_path(str(HERE / "verificar_H5A_3_4_balance_por_pieza.py"))
Z = BLOCK["Z"]
TERMINALS = (1, 9, 17)  # Una celda, pareja vertical, pareja diagonal.
THROUGH = (1, 3, 6)


def cells(mask):
    return tuple(i for i in range(12) if mask & (1 << i))


def check_symmetry(winners):
    # Las transposiciones adyacentes generan todas las permutaciones de cada eje.
    assert len(winners) == 771
    for axis, length in ((0, 4), (1, 3)):
        for start in range(length - 1):
            new_bits = []
            for site in range(12):
                rc = list(divmod(site, 3))
                if rc[axis] in (start, start + 1):
                    rc[axis] = 2 * start + 1 - rc[axis]
                new_bits.append(1 << (3 * rc[0] + rc[1]))
            def moved(mask):
                return sum(new_bits[k] for k in cells(mask))
            assert {(moved(h), moved(f)) for h, f in winners} == set(winners)


def subsets(sites, sizes):
    return [s for n in sizes for s in combinations(sites, n)]


@lru_cache(maxsize=None)
def other_axis(size, origin):
    """Vecinos del otro lector. El índice size representa un vecino exterior."""
    others = tuple(i for i in range(size) if i != origin)
    return [inside + (size,) * (size - 1 - len(inside))
            for inside in subsets(others, range(size))]


@lru_cache(maxsize=None)
def local_mask(source, target):
    """Llevar la celda de entrada a (0,0), conservando las igualdades de ambos ejes."""
    r, c = source
    rows, cols = {r: 0}, {c: 0}
    mask = 0
    for site in target:
        a, b = divmod(site, 3)
        i = rows.setdefault(a, len(rows))
        j = cols.setdefault(b, len(cols))
        mask |= 1 << (3 * i + j)
    return mask


def endpoint_prices(winners, rows, cols):
    """Un paso de un único padre. Cada marca de salida deberá estar en H de la reunión."""
    per_input = {}
    all_inputs = {}
    for lit, foot in winners:
        birth = foot ^ lit
        mult = 3 if foot in rows else 4 if foot in cols else 1
        price = Z ** (lit.bit_count() - 1) / mult
        for marked in subsets(cells(birth), (1, 2)):
            mask = sum(1 << i for i in marked)
            if lit & 1:
                per_input[mask] = per_input.get(mask, 0) + price
            # El lector que contiene las marcas está fijado y se elige el padre en H.
            all_inputs[marked] = all_inputs.get(marked, 0) + lit.bit_count() * price

    @lru_cache(maxsize=None)
    def fixed_input(source, target):
        # Hasta doce lectores contienen simultáneamente la entrada y las salidas marcadas.
        return 12 * per_input.get(local_mask(source, target), 0)

    @lru_cache(maxsize=None)
    def free_input(target):
        return 12 * all_inputs.get(target, 0)

    return fixed_input, free_input


def short_rows(winners, rows, cols):
    price, free = endpoint_prices(winners, rows, cols)

    @lru_cache(maxsize=None)
    def direct(lit, x):
        # La bifurcación aporta x directamente y la otra celda sólo al paso intermedio.
        sites = tuple(i for i in cells(lit) if i != x)
        marks = subsets(sites, (1, 2))
        head_r, head_c = divmod(x, 3)
        def value(r, c):
            if r < 4 and c < 3 and lit & (1 << (3 * r + c)):
                return 0
            return sum(price((r, c), mark) * Z ** (len(sites) - len(mark)) for mark in marks)
        weights = {(r, c): value(r, c) for r in range(5) if r != head_r for c in range(4)}
        vertical = max(sum(weights[r, head_c] for r in rr) for rr in other_axis(4, head_r))
        diagonal = max(sum(weights[r, c] for r, c in product(rr, cc))
                       for rr in other_axis(4, head_r) for cc in other_axis(3, head_c))
        return (12 * vertical, 12 * diagonal)

    @lru_cache(maxsize=None)
    def arbitrary_input(target):
        return max(price((r, c), target) for r in range(5) for c in range(4))

    @lru_cache(maxsize=None)
    def two_arms(lit):
        marks = subsets(cells(lit), (1, 2))
        cost = Q(0)
        for a, b in combinations(marks, 2):
            if not set(a).isdisjoint(b):
                continue
            # Una orientación basta. Se ancla el primer interior en a y el otro en b.
            weight = min(free(a) * arbitrary_input(b), free(b) * arbitrary_input(a))
            cost += weight * Z ** (lit.bit_count() - len(a) - len(b))
        return (12 * THROUGH[1] * cost, 12 * THROUGH[2] * cost)

    k01, k11 = [], []
    for terminal in TERMINALS:
        direct_row = [Q(0)] * 3
        twin_row = [Q(0)] * 3
        for lit, foot in winners:
            if foot ^ lit != terminal:
                continue
            mult = 3 if foot in rows else 4 if foot in cols else 1
            for x in cells(lit):
                for i, weight in enumerate(direct(lit, x), 1):
                    direct_row[i] += weight / mult
            for i, weight in enumerate(two_arms(lit), 1):
                twin_row[i] += weight / mult
        k01.append(direct_row)
        k11.append(twin_row)
    return k01, k11


def add(left, right):
    return [[x + y for x, y in zip(r, s)] for r, s in zip(left, right)]


def scale(factor, matrix):
    return [[factor * value for value in row] for row in matrix]


def add_long_arms(base, relaxed, side, main_paths):
    """En el origen de la bifurcación tiene que haber dos celdas de salida."""
    join = [[2 * (side[0] * relaxed[0][i][j] + sum(side[1:]) * relaxed[1][i][j])
             for j in range(3)] for i in range(3)]
    arm = BLOCK["multiply"](join, main_paths)
    return add(base, [[0, row[1], row[2]] for row in arm])


def cc_bound(paths):
    return 12 * 27 * sum(n * (sum(row) - 1) for n, row in zip(THROUGH, paths))


def nested_paths(reverse, k01, k11, relaxed, forward, prior, tilt):
    """Mayorantes hacia delante y hacia atrás para brazos que también pueden bifurcarse."""
    star, mul = BLOCK["path_sum"], BLOCK["multiply"]
    if tilt == 1:
        arm_limits = (".009746", ".000000623", ".0003633")
        witness = (("1.00973", ".08450", ".14570"),
                   (".00000271", "1.01339", ".02233"),
                   (".00012696", ".09400", "1.20631"))
    else:
        assert tilt == Q(6, 5)
        arm_limits = (".011742", ".000000834", ".0004639")
        witness = (("1.0117221", ".209302", ".4053589"),
                   (".0000119", "1.0478638", ".090665"),
                   (".0002086", ".2913711", "1.6483681"))
    arms = tuple(map(Q, arm_limits))
    step = [row[1:] for row in prior["fork_rows"](Z, tilt, arms)]
    forward_paths = star(step)
    nonempty = [[x - int(i == j) for j, x in enumerate(row)]
                for i, row in enumerate(forward_paths)]
    assert all(x <= y for x, y in zip(nonempty[0], arms))
    # Desde una sola celda el primer paso es lineal. Se exige al menos uno más.
    side = mul(scale(tilt, forward), nonempty)[0]
    base = add(scale(tilt, reverse), add(scale(tilt**2, k01), scale(tilt**3, k11)))
    paths = [list(map(Q, values)) for values in witness]
    # La reunión añade una ocurrencia. La bifurcación se cobra en su propio prefijo.
    kernel = add_long_arms(base, relaxed, [tilt * x for x in side], paths)
    iterate = mul(kernel, paths)
    assert all(int(i == j) + iterate[i][j] < paths[i][j]
               for i, j in product(range(3), repeat=2))
    if tilt == 1:
        assert side[0] < Q(".000132503") and sum(side[1:]) < Q(".000050786")
        assert max(map(sum, kernel)) < Q(".247")
    return paths


def contact_bounds(incoming, lam):
    """Aplicar la ley de las consultas de la pieza al mismo campo h_CC de v24."""
    tau = Q(lam, 10**6)
    assert tau < 1
    rate = Z * lam * (1 + tau / (2 * (1 - tau / 3)))
    by = 23454 + 293 * incoming
    assert lam - 204 * rate >= by
    out = []
    for v, c in ((1, 3), (2, 2)):
        u = 2 * v + 3 * c
        qbase = 210 * v + 252 * c + 114 * u
        dbase = 16641 * v + 19665 * c + 8748 * u
        ey = Z * qbase / (1 - 204 * Z)
        mean = dbase + incoming * qbase + by * ey
        slope = dbase + (incoming + rate) * qbase
        assert mean < ceil(mean) and slope < ceil(slope)
        out.append((ceil(mean), ceil(slope)))
    return out


def main():
    if not __debug__:
        raise RuntimeError("Ejecutar con las comprobaciones de Python activadas.")
    winners, rows, cols = run_path(str(HERE / "verificar_H5A_3_4_cadenas_un_predecesor.py"))["catalogo"]()
    check_symmetry(winners)
    prior = run_path(str(HERE / "verificar_H5A_3_4_reunion_de_ramas.py"))
    forward = [row[1:] for row in prior["causal_rows"](Z)[0]]
    reverse = BLOCK["reverse_rows"](winners, rows, cols)
    relaxed = [BLOCK["reverse_rows"](winners, rows, cols, a) for a in (1, 2)]
    star, mul = BLOCK["path_sum"], BLOCK["multiply"]
    side = mul(mul(forward, forward), star(forward))[0]
    kernel = add_long_arms(reverse, relaxed, side, star(reverse))
    k01, k11 = short_rows(winners, rows, cols)
    for label, short, ceiling, count in (
        ("Largas y (0,1)", k01, ".160509", 421),
        ("Largas, (0,1) y (1,1)", k11, ".160627", 422),
    ):
        kernel = add(kernel, short)
        largest = max(map(sum, kernel))
        incoming = cc_bound(star(kernel))
        assert Q(ceiling) - Q("0.000001") < largest < Q(ceiling)
        assert count - 1 < incoming < count
        print(f"{label}: máxima suma de fila {float(largest):.12f} < {ceiling}, orígenes CC < {count}.")
    assert kernel[2][2] < Q(".110")
    paths = nested_paths(reverse, k01, k11, relaxed, forward, prior, Q(1))
    nested = cc_bound(paths)
    assert 696 < nested < 697
    print(f"Anidaciones de cualquier profundidad: supersolución racional, orígenes CC < {ceil(nested)}.")
    tilted = cc_bound(nested_paths(reverse, k01, k11, relaxed, forward, prior, Q(6, 5)))
    assert 2164 < tilted < 2165
    print("Por check de salida, cola con al menos n ocurrencias cobradas: < 2165 (5/6)^n.")
    for label, count, lam, expected in (
        ("serie", 422, 205000, [(1627606, 1714907), (1513539, 1594705)]),
        ("serie-paralelo", 697, 322000, [(2522338, 2731725), (2345397, 2540071)]),
    ):
        values = contact_bounds(count, lam)
        assert values == expected
        for group, (mean, slope) in zip(("1V+3C", "2V+2C"), values):
            print(f"{label}, {group}: media < {mean}, pendiente < {slope}.")
    print("El campo de contactos satisface 0 <= h_CC <= 10^-6. No se evalúa una nueva esquina global.")


if __name__ == "__main__":
    main()

