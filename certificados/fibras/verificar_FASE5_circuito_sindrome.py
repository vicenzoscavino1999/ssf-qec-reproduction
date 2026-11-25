"""HGP (3,4). Certificado racional por los checks que un fallo de circuito puede encender."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement, permutations, product
from math import ceil, comb, factorial, floor
from pathlib import Path
from runpy import run_path

HERE = Path(__file__).resolve().parent
RM = Q(".001")
B_PASS, THETA, ALPHA, ETA, TILT = map(Q, (".001145", ".00115", ".002", ".0001", ".015"))
GRID = 10**26
AXIS_PAIRS = ((0, 1), (0, 2), (1, 2))
OVERLAP = ((9, 9), (16, 8))
CANONICAL = (("1.068897711e-6", "1.223239294e-6"), ("3.442391362e-7", "6.435991854e-7"),
             ("4.133824791e-8", "2.340104217e-7"), ("3.627668276e-9", "2.374346905e-8"),
             "7.293164044e-7", "1.179440849e-6")
INTERLEAVED = (("9.988591961e-7", "1.143426562e-6"), ("5.072115042e-7", "6.089217595e-7"),
               ("6.864893412e-8", "2.746819638e-7"), ("5.621632945e-9", "2.990595467e-8"),
               "6.544026537e-7", "1.057394828e-6")


def previous(name):
    return run_path(str(HERE / f"verificar_H5A_3_4_{name}.py"))


def rounded(x):
    return Q(ceil(x*GRID), GRID)


def exp_upper(x):
    # El resto geométrico del Taylor de grado 24 mayoriza la exponencial.
    return rounded(taylor(x))


taylor = previous("medidas_y_un_dato")["exp_upper"]


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def evaluate(poly, w):
    return sum(x*w[0]**v*w[1]**c for (v, c), x in poly.items())


def laws(kind, rate, first):
    n = 3+kind
    wait, ctrl = 2*rate/(3-2*rate), 8*rate/(15-8*rate)
    one = (12-n)*wait+n*ctrl
    total = (1+wait)**(12-n)*(1+ctrl)**n-1
    many = 1-one/total
    assert 0 < many < Q("3e-7")
    # Las dosis de un solo control, por número de celdas omitidas en el observado.
    masses = [first*wait+ctrl/2] + [ctrl]*(n-1) + [(12-n-first)*wait+ctrl/2]
    assert sum(masses) == one
    # La h-ésima arista procesada tiene a lo sumo (first-1)+(h-1) colores prohibidos.
    # Probar también la dominación acumulada de las esperas para cualquier secuencia admisible.
    for gates in combinations_with_replacement(range(1, 11), n):
        if len(set(gates)) < n or any(g > first+j for j, g in enumerate(gates)):
            continue
        waits = [0]*(n+1)
        for t in range(12):
            if t not in gates:
                waits[sum(g <= t for g in gates)] += 1
        for h in range(n):
            assert sum(waits[:h+1]) <= first

    regular = (2, 3, 6) if kind == 0 else (2, 2, 4, 6)
    isolated = (2, 3, 6) if kind == 0 else (2, 2, 6, 6)
    capacity, general, disjoint, weak = (Counter() for _ in range(4))
    for b in range(n+1):
        prob = comb(n, b)*THETA**b*(1-THETA)**(n-b)
        for h, dose in enumerate(masses):
            visible = n-h
            limit = min(n, visible+b)
            weight = dose*prob/total
            capacity[limit] += weight
            general[6*visible-sum(regular[:visible])+4*b] += weight
            if 2*limit >= n:
                disjoint[6*(2*limit-n)-sum(regular[:limit])] += weight
                weak[6*(2*limit-n)-sum(isolated[:limit])] += weight
        general[6*n-sum(regular)+4*b] += many*prob
        disjoint[6*n-sum(regular)] += many*prob
        weak[6*n-sum(isolated)] += many*prob
    assert sum(capacity.values()) == 1-many and sum(general.values()) == 1
    return total, many, capacity, general, disjoint, weak


def overlap_counts(trees, neighbors):
    """Un árbol enraizado con una arista de solapamiento distinguida."""
    result = [Counter(), Counter()]
    for n in range(2, 11):
        for i in (0, 1):
            attach = Counter()
            for j in (0, 1):
                for key, mu in trees[j].items():
                    attach[key] += OVERLAP[i][j]*mu
                for key, mu in result[j].items():
                    attach[key] += neighbors[i][j]*mu
            for (v, c), a in trees[i].items():
                for (vv, cc), b in attach.items():
                    if v+c+vv+cc == n:
                        result[i][v+vv, c+cc] += a*b
    return result


def triple_counts(isolated, profiles, stable, rooted_trees):
    counts, overlaps = Counter(), Counter()
    trees = tuple(rooted_trees(3))

    @lru_cache(None)
    def cap(types, relations):
        return isolated["triple_budget"](types, relations, profiles, stable)

    for types in product((0, 1), repeat=3):
        options = [(0,)+tuple(isolated["positions"](types[i], types[j])) for i, j in AXIS_PAIRS]
        for edges in product(*options):
            if sum(bool(e) for e in edges) < 2:
                continue
            if len(set(types)) == 1 and sum(e >= 3 for e in edges) == 2 and 0 not in edges:
                continue
            relation = dict(zip(AXIS_PAIRS, edges))
            key = min((tuple(types[i] for i in order),
                       tuple(relation[tuple(sorted((order[i], order[j])))] for i, j in AXIS_PAIRS))
                      for order in permutations(range(3)))
            # VV→CC y CC→VV tienen multiplicidades distintas. Orientar cada
            # árbol desde el dato fijado antes de contar las dos ubicaciones.
            bounds = []
            for tree in trees:
                count = 1
                for parent, child in tree:
                    rel = relation[tuple(sorted((parent, child)))]
                    if rel == 0:
                        break
                    count *= isolated["positions"](types[parent], types[child])[rel]
                else:
                    bounds.append(count)
            weight = Q(min(bounds), 2)
            address = types[0], types.count(0), types.count(1), cap(*key)
            counts[address] += weight
            if 4 in edges:
                overlaps[address] += weight
    assert cap.cache_info().currsize == 132
    return counts, overlaps


@lru_cache(None)
def rank_caps(nv, nc, profiles, donors):
    """Saldo condicional para cualquier vector de capacidades de los lectores."""
    widths = (3,)*nv+(4,)*nc

    def sorted_types(vec):
        return tuple(sorted(vec[:nv]))+tuple(sorted(vec[nv:]))

    @lru_cache(None)
    def debit(a, rank):
        if sum(bool(k) for k in a) > donors:
            return 10**6
        count = sum(a)
        row, col = max(a[:nv], default=0), max(a[nv:], default=0)
        return min((2*h-v for h, v, hr, hc, _ in profiles
                    if rank*v <= 12*h and count <= h and row <= hr and col <= hc
                    and (rank < 12 or (v-count <= 3 and 2*count >= v))), default=10**6)

    @lru_cache(None)
    def cost(r):
        if not any(r):
            return 0
        rank = max(12*k//n for k, n in zip(r, widths))
        best = 10**6
        for a in product(*(range(k+1) for k in r)):
            if not any(a) or any(j != nv and r[j] == r[j-1] and a[j] < a[j-1]
                                 for j in range(1, len(r))):
                continue
            d = debit(sorted_types(a), rank)
            if d < best:
                best = min(best, d+cost(sorted_types(tuple(k-j for k, j in zip(r, a)))))
        return best

    caps = {}
    for vv, cc in product(combinations_with_replacement(range(4), nv),
                          combinations_with_replacement(range(5), nc)):
        r = vv+cc
        total, width = sum(r), sum(widths)
        debit_r = cost(r)
        general = max(0, total-(width-total)%2-debit_r)
        disjoint = max(0, 2*total-width-debit_r) if all(2*k >= n for k, n in zip(r, widths)) else -1
        for j, k in enumerate(r):
            if k:
                old = caps[sorted_types(tuple(x-(i == j) for i, x in enumerate(r)))]
                general, disjoint = max(general, old[0]), max(disjoint, old[1])
        caps[r] = general, disjoint
    return caps


def convolve_hist(left, right, merge):
    new = Counter()
    for key, x in left.items():
        for other, y in right.items():
            new[merge(key, other)] += x*y
    return {k: rounded(x) for k, x in new.items()}


def polys(hist, trees, crowded, discounts, geometry, triple, laws_by_type, neighbors, anchor, ball, Y):
    overlap = overlap_counts(trees, neighbors)
    ordinary, marked = [Counter(), Counter()], [Counter(), Counter()]
    general, disjoint = {(0, 0): {0: Q(1)}}, {(0, 0): {0: Q(1)}}
    tuples = [[{(): Q(1)}] for _ in (0, 1)]
    for i in (0, 1):
        for _ in range(6):
            tuples[i].append(convolve_hist(tuples[i][-1], laws_by_type[i][2],
                                           lambda x, k: tuple(sorted(x+(k,)))))
    for n in range(1, 11):
        for v in range(n+1):
            c = n-v
            key = v, c
            axis = 0 if v else 1
            sub = (v-1, c) if v else (v, c-1)
            for dist, pos in ((general, 3), (disjoint, 4)):
                dist[key] = convolve_hist(dist[sub], laws_by_type[axis][pos], int.__add__)
            by_donor = [rank_caps(v, c, geometry, j) for j in range(2, n+1)] if 2 <= n <= 6 else []
            for i in (0, 1):
                if key not in hist[i]:
                    continue
                chi = discounts.get(key, Q(1))
                if n == 1:
                    ordinary[i][key] = chi*Y*sum(laws_by_type[i][4].values())
                    continue
                geom = hist[i][key]
                weights, limits = [], []
                for k in range(max(geom)+1):
                    total = sum(x for h, x in geom.items() if h >= k)
                    with_edge = min(total, overlap[i][key])
                    if n == 3:
                        totals = [sum(x for (ix, vv, cc, h), x in table.items()
                                      if (ix, vv, cc) == (i, v, c) and h >= k) for table in triple]
                        total = min(total, totals[0])
                        with_edge = min(total, with_edge, totals[1])
                    def event(dist):
                        return min(Q(1), sum(x for h, x in dist.items() if k == 0 or h >= 6*k))
                    pg, pa = event(disjoint[key]), event(general[key])
                    weights.append((total-with_edge)*pg+with_edge*max(pa, pg))
                    volumes = [total, total, total]+[
                        min(total, crowded[j][i].get(key, 0)) for j in range(3, n+1) if by_donor]+[0]
                    limits.append((total, with_edge, volumes))
                if by_donor:
                    prob_many = 1-(1-laws_by_type[0][1])**v*(1-laws_by_type[1][1])**c
                    conditional = [t*prob_many for t, _, _ in limits]
                    for vv, pv in tuples[0][v].items():
                        for cc, pc in tuples[1][c].items():
                            row = [tab[vv+cc] for tab in by_donor]
                            for k, (_, edge, volumes) in enumerate(limits):
                                ja = next((j for j, (a, _) in enumerate(row, 2) if a >= k), n+1)
                                jg = next((j for j, (_, g) in enumerate(row, 2) if g >= k), n+1)
                                conditional[k] += pv*pc*min(volumes[ja], volumes[jg]+edge)
                    weights = [min(x, y) for x, y in zip(weights, conditional)]
                assert all(x >= y >= 0 for x, y in zip(weights, weights[1:]))
                ordinary[i][key] = chi*(weights[0]+sum(x*(Y**(-k)-Y**(1-k))
                                                        for k, x in enumerate(weights[1:], 1)))
                marked[i][key] = chi*sum(x*(anchor(ball, i, k)-anchor(ball, i, k-1))
                                          for k, x in enumerate(weights[1:], 1))
    return ordinary, marked


def bulk_tail(low, high, weak_edges, strong_edges):
    """Supersolución racional positiva y extracción de los diez grados en los mismos pesos."""
    def rows(d, e):
        result = []
        for i in (0, 1):
            x, z = dot(weak_edges[i], d), dot(strong_edges[i], e)
            result.append((rounded(exp_upper(x)*(low[i]+high[i]*(exp_upper(z)-1))),
                           rounded(high[i]*exp_upper(x+z))))
        return tuple(zip(*result))

    d = e = (Q(0), Q(0))
    for _ in range(600):
        nd, ne = rows(d, e)
        if (nd, ne) == (d, e):
            break
        assert all(a <= b for a, b in zip(d+e, nd+ne)) and max(ne) < Q(".01")
        d, e = nd, ne
    else:
        raise AssertionError("La supersolución de la cola no estabilizó.")
    ds, es = [[Q(0)], [Q(0)]], [[Q(0)], [Q(0)]]
    wf, ff = [[Q(1)], [Q(1)]], [[Q(1)], [Q(1)]]
    for k in range(1, 11):
        for i in (0, 1):
            ds[i].append(low[i]*wf[i][k-1]+high[i]*(ff[i][k-1]-wf[i][k-1]))
            es[i].append(high[i]*ff[i][k-1])
        for i in (0, 1):
            def coeff(j, full):
                parents = weak_edges[i][0]*ds[0][j]+weak_edges[i][1]*ds[1][j]
                if full:
                    parents += strong_edges[i][0]*es[0][j]+strong_edges[i][1]*es[1][j]
                return j*parents
            for forest, full in ((wf, False), (ff, True)):
                forest[i].append(sum(coeff(j, full)*forest[i][k-j] for j in range(1, k+1))/k)
    remainder = [d[i]-sum(ds[i]) for i in (0, 1)]
    assert min(remainder) >= 0
    return remainder


def certificate(law, poly, pure, halos, far, ball, rate, a, s, witness):
    # m, U, h, H, v y la marca de cobertura son valores prefijados, todos racionales.
    m, u, h, H = [tuple(map(Q, row)) for row in witness[:4]]
    visible, root = map(Q, witness[4:])
    Y = a**6
    target = 4*rate/(15-8*rate)
    odds_check = (1-2*rate/3)**(-5)/(1-RM)/(1-target)**7-1
    y = odds_check*exp_upper(TILT)
    z = [v[0]*exp_upper(TILT+ALPHA) for v in law]
    assert y/(1+y) <= THETA
    R = max(dot(row, u) for row in far["N"])
    g = 40*H[0]+27*H[1]+root
    I = 280*u[0]+189*u[1]+588*g
    ell = 89*B_PASS+I
    eps = exp_upper(ell)*(exp_upper(ell)-1)
    assert eps < 1 and R < ALPHA and max(u) < ETA

    def upper_log(t):
        return t+t*t/(2*(1-t))

    collision = Q(89, 2)*upper_log(THETA**2)
    guard = upper_log(THETA)-collision
    ordinary = B_PASS+collision-guard
    assert guard >= Q(".001") and ordinary > 0
    ratios = [y*(1+B_PASS)*(1+eps)/B_PASS,
              y*(1-THETA)*exp_upper(99*guard+90*ordinary+2*I)/THETA]
    dd, dg, gd, gg = far["kernels"](ball, Y, s)
    w, depart = [], []
    for i in (0, 1):
        exposure = 2*dot([row[i] for row in halos[1:]], (u[0], u[1], g))
        rhs = z[i]*exp_upper(ETA+2*R+guard*(504, 628)[i]+ordinary*halos[0][i]+exposure)
        ratios.append(rhs/m[i])
        w.append(rounded(m[i]*exp_upper(dot([k[0] for k in dd[i]], h)+dg[i][0]*visible)))
        depart.append(dot([k[1] for k in dd[i]], h)+dg[i][1]*visible)

    low = [rounded(w[i]*sum(p*a**(-k) for k, p in law[i][5].items())) for i in (0, 1)]
    high = [rounded(w[i]*sum(p*a**(-k) for k, p in law[i][3].items())) for i in (0, 1)]
    assert all(a <= b for a, b in zip(low, high))
    tails = bulk_tail(low, high, far["W"], far["B"])
    trees = [mi*Q("1.014") for mi in m]
    init = (210, 252), (640, 784), (441, 513), (1344, 1596)
    removed = []
    for i in (0, 1):
        group = evaluate(poly[i], w)
        naked = evaluate(poly[i], m)
        removed.append(naked)
        pointed = sum(val*w[0]**v*w[1]**c*(v*depart[0]+c*depart[1]) for (v, c), val in poly[i].items())
        connected = exp_upper(dot(far["N"][i], trees))
        assert m[i]*connected <= trees[i]
        fixed = min(group-naked, pointed)+tails[i]+m[i]*(connected-1)
        start = dot((B_PASS, 2*u[0], 2*u[1], 2*g), [row[i] for row in init])
        simple = z[i]*exp_upper(ETA+2*R+start)*(1+y)**(3+i)*sum(law[i][4].values())
        ratios.extend(((group+tails[i])/h[i], fixed/H[i], (simple+H[i])/u[i]))
    J = dot([k[0] for k in gd], h)+gg[0]*visible
    one = dot([k[1] for k in gd], h)+gg[1]*visible
    raw = exp_upper(J)*one-dot([k[1] for k in gd], removed)
    assert raw > 0
    outer = Y/s**2
    largest = max(m)
    ratio = exp_upper(Q(1))*1297*largest*Q(12, 11)**9
    assert ratio < 1
    direct_tail = 2*2000*11**9*Q(11**9, factorial(10))*largest*(1297*largest)**10/(1-ratio)
    distant = outer*(outer*raw+(exp_upper(J)-1)*one/2)+direct_tail+sum(evaluate(p, m) for p in pure)
    ratios.extend((outer*(exp_upper(J)-1-J)/visible, distant/root))
    for i, ratio in enumerate(ratios):
        assert 0 < ratio < 1, (i, float(ratio))
    return R, max(ratios), tails, ratios[0], ratios[1]


def early_controls_at_target(far):
    # Condición necesaria de estas filas con a=.56 y el precio de cobertura >=.001.
    def lower_exp(x):
        return sum(x**j/factorial(j) for j in range(7))

    def down(x):
        return Q(floor(x*10**16), 10**16)

    rate = Q("3.75e-7")
    base = [rate*Q(c, 15)*lower_exp(Q(gamma, 1000))
            for c, gamma in ((44, 504), (74, 628))]
    a = Q(".56")
    lo = [base[0]/a**7, base[1]/a**8]
    hi = [base[0]/a**7, base[1]/a**10]
    d = e = (Q(0), Q(0))
    for step in range(1, 51):
        x = [dot(line, d) for line in far["W"]]
        z = [dot(line, e) for line in far["B"]]
        d = [down(lower_exp(x[i])*(lo[i]+hi[i]*(lower_exp(z[i])-1))) for i in (0, 1)]
        e = [down(hi[i]*lower_exp(x[i]+z[i])) for i in (0, 1)]
        if e[1] >= Q(".1"):
            assert hi[1]*lower_exp(Q("12.9")) > Q(".1")
            return step
    raise AssertionError("No se alcanzó la obstrucción del mayorante.")


def main():
    far, history, contacts, stripes, base, isolated, prior = map(previous, (
        "propagacion_lejana", "medidas_por_historia", "contacto_de_checks", "franjas",
        "medidas_y_un_dato", "medidas_aisladas", "exceso_sin_compensar"))
    profiles, stable, reps = contacts["geometry"]()
    wins = far["winners"](reps)
    assert len(wins) == 771
    assert all(max(((foot ^ lit) & row).bit_count() for row in far["ROWS"]) <= 1
               and (foot ^ lit).bit_count() <= 2 for lit, foot in wins)
    history["far_credits"](profiles, stable)
    assert all(0 <= 6-credit <= 4 for cr in ((2, 3, 6), (2, 2, 4, 6)) for credit in cr)
    _, cap = stripes["cost_tables"](stripes["light_profiles"](), previous("fondo_pasivo")["step_costs"],
                                    previous("pares_adyacentes")["stripe_profiles"])
    small = history["data_coefficients"](contacts["positions"](profiles, stable, cap, stripes["coefficients"](4)))
    triples = isolated["count_triples"](profiles, stable, contacts["rooted_trees"])
    trees, six = stripes["coefficients"](10), stripes["coefficients"](6)
    hist = prior["cumulative"](trees, small, triples, isolated["group_series"](10), far["pair_counts"](wins, stable))
    prior["refine_tails"](hist, six, reps)
    profiles_f = tuple(sorted(far["cover_profiles"](reps)))
    discounts = {}
    for n in range(1, 7):
        for v in range(n+1):
            c = n-v
            saving = max(0, 246*v+340*c-far["group_budget"](v, c, profiles_f, n, guard=True))
            discounts[v, c] = 1/sum(Q(saving, 1000)**j/factorial(j) for j in range(6))
    crowded = far["crowded_coefficients"](six)
    trip = triple_counts(isolated, profiles, stable, contacts["rooted_trees"])
    for i, table in enumerate(triples):
        for (v, c, k), count in table.items():
            address = i, v, c, k
            assert 0 <= trip[1][address] <= trip[0][address] <= count
            if v and c:
                assert trip[0][address] == count
    assert all((v, c, k) in triples[i] for i, v, c, k in trip[0])
    halos = contacts["exposure"](reps, base["budgets"])
    for rate, a, s, first, witness in ((Q("7.5e-8"), Q(".56"), Q(".76"), 4, CANONICAL),
                                       (Q("7e-8"), Q(".555"), Q(".78"), 7, INTERLEAVED)):
        laws_by_type = [laws(0, rate, first), laws(1, rate, 7)]
        regular, marked = polys(hist, trees, crowded, discounts, profiles_f, trip, laws_by_type, far["N"],
                                prior["anchor"], history["ball"], a**6)
        row, highest, tails, passive, active = certificate(laws_by_type, regular, marked, halos, far,
                                                          history["ball"], rate, a, s, witness)
        pv, pc = [sum(x[4].values()) for x in laws_by_type]
        if first == 4:
            assert row < Q(".000601") and highest < Q(".99883")
            assert tails[0] < Q("6e-11") and tails[1] < Q("6e-10")
            assert pv < Q(".457") and pc < Q(".805")
        else:
            assert row < Q(".000711") and highest < Q(".99894")
        order = "VV primero" if first == 4 else "cualquier orden de los datos"
        print(f"{order}, rc<={float(rate):g}, rm<=1e-3. Fila {float(row):.12g} sobre {float(ALPHA):g}.")
        print(f"Mayor cociente {float(highest):.12f}. Pasivo {float(passive):.12f}, activo {float(active):.12f}.")
        print(f"Restos desde once por raíz VV y CC {float(tails[0]):.12g}, {float(tails[1]):.12g}.")
        print(f"Factores para una franja disjunta VV y CC {float(pv):.9g}, {float(pc):.9g}.")
    print(f"En rc=3.75e-7, los controles de franja completa por sí solos excluyen estas filas "
          f"tras {early_controls_at_target(far)} iteraciones inferiores.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para verificar las desigualdades.")
    main()

