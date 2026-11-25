"""HGP (3,4). Balance adaptativo de una pieza y contactos con brazos entrantes.
Cálculo exacto con la biblioteca estándar y el catálogo físico del proyecto.
"""

from fractions import Fraction as Q
from math import ceil
from pathlib import Path
from runpy import run_path


HERE = Path(__file__).resolve().parent
Z = Q(12, 10000)
THROUGH = (1, 3, 6)


def axis_split(degree, radius):
    """Bolas de la paridad final desde los vecinos de una raíz del árbol."""
    adj = [set()]
    level = [0]
    for depth in range(radius + 1):
        following = []
        for vertex in level:
            valence = degree if depth % 2 == 0 else 7 - degree
            for _ in range(valence - len(adj[vertex])):
                child = len(adj)
                adj.append({vertex})
                adj[vertex].add(child)
                following.append(child)
        level = following
    balls = []
    for root in adj[0]:
        seen = {root}
        for _ in range(radius):
            seen = {u for v in seen for u in adj[v]}
        balls.append(seen)
    common = set.intersection(*balls)
    private = [ball - common for ball in balls]
    assert len(set.union(*private)) == sum(map(len, private))
    private_sizes = {len(part) for part in private}
    assert len(private_sizes) == 1
    return len(common), private_sizes.pop()


def check_cover(winners, rows, cols, row_radius, col_radius):
    ar, br = axis_split(4, row_radius)
    ac, bc = axis_split(3, col_radius)
    def size(mask):
        if mask == 0:
            return 0
        nr = sum(bool(row & mask) for row in rows)
        nc = sum(bool(col & mask) for col in cols)
        return ar * ac + br * ac * nr + ar * bc * nc + br * bc * mask.bit_count()
    costs = [Q(size(f) - size(h), 2 * h.bit_count() - f.bit_count()) for h, f in winners]
    return size(rows[0]), size(cols[0]), size(1), max(costs)


def reverse_rows(winners, rows, cols, extra=0):
    """Mayorante de los padres de una ocurrencia por lector y máscara terminal."""
    singles = tuple(1 << j for j in range(12))
    pairs = [s for s in range(4096) if s.bit_count() == 2]
    vertical = tuple(s for s in pairs if any(s & col == s for col in cols))
    diagonal = tuple(s for s in pairs if not any(s & line == s for line in rows + cols))
    states = (singles, vertical, diagonal)
    assert tuple(map(len, states)) == (12, 18, 36)
    assert all(sum(bool(s & bit) for s in group) == k
               for group, k in zip(states, THROUGH) for bit in singles)
    assert {f ^ h for h, f in winners} == {0, *singles, *vertical, *diagonal}

    result = []
    for terminals in states:
        signatures = set()
        for terminal in terminals:
            row = [Q(0) for _ in states]
            for h, f in winners:
                if f ^ h != terminal:
                    continue
                n = h.bit_count()
                mult = 3 if f in rows else 4 if f in cols else 1
                one = Z ** max(0, n - 1 - extra)
                two = Z ** max(0, n - 2 - extra)
                for i, parents in enumerate(states):
                    mass = THROUGH[i] * n * one
                    if i:
                        doubles = sum(s & h == s for s in parents)
                        mass += doubles * max(0, two - 2 * one)
                    row[i] += 12 * mass / mult
            signatures.add(tuple(row))
        assert len(signatures) == 1
        result.append(signatures.pop())
    return result


def check_incoming(rows):
    caps = tuple(map(Q, ("1.1141", "1.00040", "1.0934")))
    tilted = tuple(map(Q, ("1.191", "1.001", "1.156")))
    for tilt, values in ((Q(1), caps), (Q(8, 5), tilted)):
        assert all(1 + tilt * sum(r * g for r, g in zip(row, values)) <= limit
                   for row, limit in zip(rows, values))
    # Un lector fijado toca franjas de hasta 3(1+4(3-1)) originales CC.
    root_choices = 3 * (1 + 4 * (3 - 1))
    per_check = 12 * root_choices * sum(n * (g - 1) for n, g in zip(THROUGH, caps))
    assert Q(218) < per_check < 219
    return ceil(per_check)


def multiply(left, right):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*right)]
            for row in left]


def path_sum(matrix):
    """Resuelve B=I+MB y verifica la positividad para dominar todas las longitudes."""
    size = len(matrix)
    identity = [[Q(i == j) for j in range(size)] for i in range(size)]
    aug = [[identity[i][j] - matrix[i][j] for j in range(size)] + identity[i]
           for i in range(size)]
    for i in range(size):
        pivot = aug[i][i]
        assert pivot
        aug[i] = [x / pivot for x in aug[i]]
        for j in range(size):
            if j != i:
                coeff = aug[j][i]
                aug[j] = [x - coeff * y for x, y in zip(aug[j], aug[i])]
    bound = [row[size:] for row in aug]
    assert all(x >= 0 for row in bound for x in row)
    product = multiply(matrix, bound)
    assert all(bound[i][j] == identity[i][j] + product[i][j]
               for i in range(size) for j in range(size))
    return bound


def reconvergences(winners, rows, cols, reverse):
    prior = run_path(str(HERE / "verificar_H5A_3_4_reunion_de_ramas.py"))
    unary, _ = prior["causal_rows"](Z)
    forward = [row[1:] for row in unary]
    all_forward = path_sum(forward)
    all_reverse = path_sum(reverse)
    joins = [reverse_rows(winners, rows, cols, extra) for extra in (1, 2)]
    power = forward
    for minimum in (1, 2):
        side = multiply(power, all_forward)[0]
        prices = (side[0], side[1] + side[2])
        join = [[2 * sum(w * mat[i][j] for w, mat in zip(prices, joins))
                 for j in range(3)] for i in range(3)]
        main_arm = multiply(join, all_reverse)
        kernel = [[reverse[i][j] + (main_arm[i][j] if j in (1, 2) else 0)
                   for j in range(3)] for i in range(3)]
        if minimum == 1:
            # Este mayorante relaja la posición final del brazo secundario.
            assert kernel[2][2] > Q("1.58")
            power = multiply(power, forward)
            continue
        assert prices[0] < Q(".000099872") and prices[1] < Q(".000022719")
        caps = tuple(map(Q, ("1.16470", "1.01336", "1.17105")))
        assert all(1 + sum(x * g for x, g in zip(row, caps)) <= cap
                   for row, cap in zip(kernel, caps))
        per_check = 12 * 27 * sum(n * (g - 1) for n, g in zip(THROUGH, caps))
        assert 398 < per_check < 399
        return ceil(per_check)


def main():
    if not __debug__:
        raise RuntimeError("Se necesitan las comprobaciones de Python activadas.")
    filename = HERE / "verificar_H5A_3_4_cadenas_un_predecesor.py"
    winners, rows, cols = run_path(str(filename))["catalogo"]()
    assert len(winners) == 771
    assert all(h & f == h and 2 * h.bit_count() > f.bit_count() for h, f in winners)
    qv, qc, qsite, qcost = check_cover(winners, rows, cols, 2, 2)
    dv, dc, dsite, dcost = check_cover(winners, rows, cols, 5, 4)
    assert (qv, qc, qsite, qcost) == (210, 252, 90, 114)
    assert (dv, dc, dsite, dcost) == (16641, 19665, 7353, 8748)
    reverse = reverse_rows(winners, rows, cols)
    linear = check_incoming(reverse)
    joined = reconvergences(winners, rows, cols, reverse)
    # Sólo para diagnosticar la relajación que regala el segundo padre.
    for bit in (1 << j for j in range(12)):
        three = [(h, f) for h, f in winners if f in rows and f ^ h == bit]
        assert len(three) == 1
        lit, _ = three[0]
        assert Q(lit.bit_count() * 12, 3) == 8
    bq = qsite + qcost
    assert bq == 204 and bq * Z < 1
    assert Q(".00115") < Z
    rc = Q(1, 10**6)
    assert Q(".001") + 5 * 2 * rc / 3 + 7 * 4 * rc / (15 - 8 * rc) < Z
    print("771 pares. Incrementos por unidad de descenso:", qcost, dcost)
    print("Orígenes CC por check terminal, recorridos no vacíos. Lineal:", linear,
          "con reuniones en serie:", joined)
    print(" Tipo        V  C | crédito  E|Q|   zonas CC  E[CC con brazos]  pendiente exponencial")
    for name, incoming, lam, expected in (
        ("lineal", linear, 124000, ((967131, 1035140), (899476, 962706))),
        ("reuniones", joined, 204000, ((1552774, 1740378), (1443965, 1618386))),
    ):
        by = dsite + dcost + dsite + incoming * (bq + qsite - 1)
        hmax = Q(1, 10**6)
        rate = Z * lam / (1 - lam * hmax)
        assert lam * hmax < 1 and lam - bq * rate >= by
        for (v, c), limits in zip(((1, 3), (2, 2)), expected):
            credit = 2 * v + 3 * c
            qbase = qv * v + qc * c + qcost * credit
            dbase = dv * v + dc * c + dcost * credit
            mean_ones = Z * qbase / (1 - bq * Z)
            base = dbase + incoming * qbase
            values = (qbase + bq * mean_ones, dbase + (dsite + dcost) * mean_ones,
                      base + by * mean_ones, base + rate * qbase)
            upper = tuple(map(ceil, values))
            assert upper[:2] == ((2940, 228661) if v == 1 else (2734, 212898))
            assert upper[2:] == limits and all(x < y for x, y in zip(values, upper))
            print(f" {name:9s} {v}  {c} | {credit:7d} {upper[0]:5d} {upper[1]:10d}"
                  f" {upper[2]:17d} {upper[3]:21d}")
    print("Campo CC entre 0 y 10^-6. La tabla cubre la clase del lema, no la fila completa de G3.")


if __name__ == "__main__":
    main()

