"""HGP (3,4). Posiciones conjuntas de las salidas de un grupo de datos."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from math import ceil, comb, factorial
from pathlib import Path
from runpy import run_path


HERE = Path(__file__).resolve().parent
old = run_path(str(HERE / "verificar_H5A_3_4_reunion_de_ramas.py"))
prior = old["prior"]
RADII = (((4, 6), (5, 5), (3, 4)), ((5, 5), (6, 4), (4, 3)))
DISTANCES = (((2, 4), (3, 3)), ((3, 3), (4, 2)))
DIRECT = ((640, 441, 30), (784, 513, 36))


def displacements(i, j):
    m, n = DISTANCES[i][j]
    return [(a, b) for a in range(m % 2, m + 1, 2)
            for b in range(n % 2, n + 1, 2) if a or b]


class SeedTree:
    """Dos centros a distancia dada y las bolas cortas que se usan en sus franjas."""
    def __init__(self, first, distance):
        self.neighbors = [[q for q in (j - 1, j + 1) if 0 <= q <= distance]
                          for j in range(distance + 1)]
        self.types = [first ^ (j % 2) for j in range(distance + 1)]
        self.left, self.right = 0, distance

    def adjacent(self, node):
        neighbors = self.neighbors[node]
        while len(neighbors) < 3 + self.types[node]:
            other = len(self.types)
            neighbors.append(other)
            self.types.append(1 - self.types[node])
            self.neighbors.append([node])
        return tuple(neighbors)

    def ball(self, node, radius):
        level = {node}
        seen = set(level)
        selected = set() if radius % 2 else set(level)
        for d in range(1, radius + 1):
            level = {y for x in level for y in self.adjacent(x)} - seen
            seen.update(level)
            if d % 2 == radius % 2:
                selected.update(level)
        return selected


@lru_cache(None)
def initial_overlap(i, j, d1, d2):
    axes = (SeedTree(i, d1), SeedTree(i, d2))
    def stripe(kind, centers):
        a, b = centers
        if kind == 0:
            return set(product((a,), axes[1].adjacent(b)))
        return set(product(axes[0].adjacent(a), (b,)))
    left = stripe(i, (0, 0))
    right = stripe(j, (d1, d2))
    def covers(sites):
        volumes = []
        for r1, r2 in ((2, 2), (2, 3), (3, 2), (3, 3)):
            cover = set()
            guards = Counter()
            for a, b in sites:
                for y in product(axes[0].ball(a, r1), axes[1].ball(b, r2)):
                    cover.add(y)
                    if (r1, r2) == (2, 2):
                        guards[y] += 1 + (a == y[0])
            if (r1, r2) == (2, 2):
                volumes.append(sum(min(2, n) for n in guards.values()))
            volumes.append(len(cover))
        return tuple(volumes)
    one, two, both = covers(left), covers(right), covers(left | right)
    expected = ((258, 210, 640, 441, 1344), (288, 252, 784, 513, 1596))
    assert one == expected[i] and two == expected[j]
    return tuple(a + b - c for a, b, c in zip(one, two, both))


@lru_cache(None)
def ball_intersection(i, distance, first, second):
    """Puntos del tipo de destino en dos bolas de un árbol (3,4)."""
    assert (first - second - distance) % 2 == 0
    total = 0
    for pos in range(distance + 1):
        kind = i ^ (pos % 2)
        length = min(first - pos, second - (distance - pos))
        parity = (first - pos) % 2
        if length >= 0 and not parity:
            total += 1
        off_path = 3 + kind - int(pos > 0) - int(pos < distance)
        count = off_path
        for depth in range(1, length + 1):
            if depth % 2 == parity:
                total += count
            count *= 2 + (kind ^ (depth % 2))
    return total


@lru_cache(None)
def overlap_layers(i, j, a, b, target, last=24):
    first = RADII[i][target]
    second = RADII[j][target]
    return tuple(ball_intersection(i, a, first[0] + 2*l, second[0] + 2*l)
                 * ball_intersection(i, b, first[1] + 2*l, second[1] + 2*l)
                 for l in range(last + 1))


@lru_cache(None)
def link_overlap(i, j, a, b, target, Y, s):
    outer, rho = Y / s**2, Y * s
    ends = 1 if target == 2 else 2
    assert 0 < outer < 1 and 0 < 36 * rho < 1 and 0 < s < 1
    weight = Q(1)
    saved = 0
    for l, area in enumerate(overlap_layers(i, j, a, b, target)):
        new = weight * (outer if l < ends else rho)
        saved += (weight - new) * area
        weight = new
    return saved


def tree_minimum(v, c, weights):
    if v == 0:
        return (c - 1) * weights[1, 1]
    if c == 0:
        return (v - 1) * weights[0, 0]
    cross = weights[0, 1]
    return cross + (v - 1) * min(weights[0, 0], cross) + (c - 1) * min(weights[1, 1], cross)


def savings(Y, s, exposure, children):
    geom, both = {}, {}
    for i, j in ((0, 0), (0, 1), (1, 1)):
        positions = displacements(i, j)
        initial = {k: sum(x*y for x, y in zip(initial_overlap(i, j, *k), exposure)) for k in positions}
        extra = {k: sum(x*link_overlap(i, j, *k, t, Y, s) for t, x in enumerate(children)) for k in positions}
        geom[i, j] = min(initial.values())
        both[i, j] = min(initial[k]+extra[k] for k in positions)
    return geom, both


def previous_tables():
    """Reconstruir las entradas de v19 desde las definiciones del proyecto."""
    v15 = prior["previous"]("verificar_H5A_3_4_contacto_oculto.py")
    fase = prior["previous"]("verificar_FASE5_circuito_sindrome.py")
    ball, halo, series, geometry, hist, _ = v15["v11"]["data_tables"]()
    discounts = v15["v11"]["geometry_discounts"](geometry)
    base, trees, crowded, disc, geom, trip, hc_halo = prior["circuit_data"](fase)
    limits = prior["triple_limits"]()
    hp = prior["refine"](hist, limits)
    hc = prior["refine"](base, limits, circuit=True)
    caps_p, caps_c = {}, {}
    for i in (0, 1):
        cp, cc = prior["aligned_limits"](hp, i)
        caps_p.update(cp)
        caps_c.update(cc)
    hp, hc = prior["apply_caps"](hp, caps_p), prior["apply_caps"](hc, caps_c)
    vol = old["triple_reader"](trees, 8)
    cp, cc = old["block_caps"](hp, vol, 8, hc)
    other_p, other_c = old["separated_caps"](hp, hc, 8)
    for caps, extra in ((cp, other_p), (cc, other_c)):
        for key, value in extra.items():
            caps[key] = min(value, caps.get(key, value))
    hp = prior["apply_caps"](hp, cp)
    hc, hr = old["circuit_histograms"](hc, hp, cc)
    return (halo, series, discounts, hp), (hc_halo, trees, crowded, disc, geom, trip, hc, hr)


def singleton_cc(rate, first, theta, Y):
    """Un aislado CC con dos luces iniciales reserva dos unidades cercanas."""
    wait, ctrl = 2*rate/(3-2*rate), 8*rate/(15-8*rate)
    mass = (1+wait)**8*(1+ctrl)**4-1
    multiple = mass-8*wait-4*ctrl
    costs = [min(Y, Y**2+12*theta*Y)]
    for h in range(1, 5):
        cost = Q(0)
        for x, y in product(range(5-h), range(h+1)):
            lit = 4-h-x+y
            prob = comb(4-h,x)*comb(h,y)*theta**(x+y)*(1-theta)**(4-x-y)
            if lit >= 2:
                cost += prob*(Y**2 if lit == 2 else Y)
        costs.append(cost)
    gates = ctrl*sum((a+b)/2 for a,b in zip(costs, costs[1:]))
    candidates = []
    for used in combinations(range(1,11), 4):
        if all(t <= first+j for j,t in enumerate(used)):
            waits = sum(costs[sum(g<=t for g in used)] for t in range(12) if t not in used)
            candidates.append((wait*waits+gates+multiple*Y)/mass)
    return max(candidates)


def pair_mark(marked, anchor, ball, discounts, circuit):
    """Raíz lejana, un solo grupo y dos datos en la misma línea de síndrome."""
    for i, (count, readers) in enumerate(((9,40), (12,27))):
        key = (2-2*i, 2*i)
        full = anchor(ball, i, 2)
        fiber = anchor(ball, i, 2, fiber=True)
        assert full > fiber > readers and anchor(ball, i, 1) == 0
        cap = discounts[key]*count*(fiber-readers)
        if circuit:
            # Toda la marca de este monomio (incluidos los controles múltiples)
            # proviene del único escalón positivo, k=2.
            cap = min(cap, marked[i][key]*Q(fiber-readers, full))
        marked[i][key] = min(marked[i][key], cap)


def upper_factor(loss):
    assert loss >= 0
    loss = Q((loss*10**12)//1, 10**12)
    den = sum(loss**j/factorial(j) for j in range(9))
    return Q(-((-10**12)//den), 10**12)


def joint_rows(regular, marked, b, theta, Y, s, u, h, physical, visible, root):
    """Filas evaluadas exclusivamente en los campos de esta supersolución."""
    log_bound = lambda t: t+t*t/(2*(1-t))
    collision = Q(89,2)*log_bound(theta**2)
    guard = log_bound(theta)-collision
    g = 40*physical[0]+27*physical[1]+root
    exposure = (guard, b+collision-guard, 2*u[0], 2*u[1], 2*g)
    children = (*h, visible)
    assert min(*exposure, *children) > 0
    initial, both = savings(Y, s, exposure, children)
    for i in (0,1):
        forbidden = sum(x*y for x,y in zip(DIRECT[i], children))
        for key in regular[i]:
            assert 1 <= sum(key) <= 10
            regular[i][key] *= upper_factor(forbidden + tree_minimum(*key, both))
        for key in marked[i]:
            assert 1 <= sum(key) <= 10
            marked[i][key] *= upper_factor(tree_minimum(*key, initial))


PHEN = dict(name="Salidas de un mismo grupo",
    rates=("5.15e-7",".001145",".004",".000912","1.01"),
    m=("9.72527903382381e-7","1.13458547975817e-6"),
    u=("6.82433427961456e-7","7.3285074143111e-7"),
    h=("3.84860828482341e-8","1.2521729434091e-7"),
    physical=("5.12743124788144e-9","2.02361963077159e-8"),
    outer=(".00115","5.67213467676418e-7","2.04376855371464e-6"),
    d=("6.77373586312163e-5",".000167676688496055"),
    e=("6.77373586312163e-5",".000451974950076228"))

CIRCUIT = (
    ("VV primero", "9.25e-8", 4, ".562", ".725",
     (("1.32168312754132e-6","1.51394347770133e-6"),
      ("4.26097369224994e-7","7.95599903956131e-7"),
      ("4.30445586918816e-8","2.02191351901184e-7"),
      ("5.57961037487384e-9","3.00985490626542e-8"),
      "6.31335859286648e-7","1.46417160304293e-6")),
    ("Cualquier orden de columnas", "8.45e-8", 7, ".56", ".74",
     (("1.20935099343879e-6","1.38584274547163e-6"),
      ("6.14034662870051e-7","7.33681137896696e-7"),
      ("6.95724636009173e-8","2.10393905326822e-7"),
      ("7.8071593508529e-9","3.36213604427016e-8"),
      "7.3562347635933e-7","1.47723059792533e-6")))


def verify_rows(data, circuit):
    v15 = prior["previous"]("verificar_H5A_3_4_contacto_oculto.py")
    fase = prior["previous"]("verificar_FASE5_circuito_sindrome.py")
    far = v15["v11"]["prop"]
    ball = v15["v11"]["previous"]("medidas_por_historia")["ball"]
    anchor = v15["v10"]["anchor"]
    halo, series, discounts, hp = data
    a, s = Q(".569"), Q(".673")
    w = PHEN
    theta, visible, root = map(Q, w["outer"])
    m, u, h, phys, d, e = (tuple(map(Q,w[k])) for k in ("m","u","h","physical","d","e"))
    reg, marked = prior["physical_prices"](v15, a, theta, ball, discounts, hp)
    pair_mark(marked, anchor, ball, discounts, False)
    joint_rows(reg, marked, Q(w["rates"][1]), theta, a**6, s, u,h,phys,visible,root)
    marked[0][0,0] = v15["single_tail"](m)
    kernels = v15["corrected_kernel"](far["kernels"](ball,a**6,s), w, reg, a**6/s**2)
    eff = v15["v11"]["effective"]
    eff["certify"](a**6, s, halo, reg, marked, kernels, lambda x,i:eff["full_tail"](series,a,x,i,d,e), w)
    halo, trees, crowded, discounts, geometry, trip, hc, hr = circuit
    theta = fase["THETA"]
    for label, rate, first, a, s, w in CIRCUIT:
        rate, a, s = map(Q,(rate,a,s))
        law = [old["circuit_law"](i, rate, first if i==0 else 7, theta) for i in (0,1)]
        reg, marked = fase["polys"](hr, trees, crowded, discounts, geometry, trip,
                                    law, far["N"], anchor, ball, a**6)
        extra, extra_marked = old["surcharge"](hc, hr, law, discounts, a**6, anchor, ball)
        for i in (0,1):
            reg[i].update(extra[i])
            marked[i].update(extra_marked[i])
            key = (1-i, i)
            cost = old["singleton_cost"](i, rate, first if i==0 else 7, theta, a**6)
            if i:
                improved = singleton_cc(rate, 7, theta, a**6)
                assert 0 < improved < cost
                cost = improved
            reg[i][key] = min(reg[i][key], discounts[key]*cost)
            assert max(hc[i][2-2*i,2*i]) <= 2 and max(hr[i][2-2*i,2*i]) <= 2
        pair_mark(marked, anchor, ball, discounts, True)
        u,h,phys = (tuple(map(Q,row)) for row in w[1:4])
        visible, root = map(Q,w[4:])
        joint_rows(reg, marked, fase["B_PASS"],theta,a**6,s,u,h,phys,visible,root)
        row, highest, tails, passive, active = fase["certificate"](law,reg,marked,halo,far,ball,rate,a,s,w)
        print(f"{label}, rc={float(rate):g}, rm=.001. Fila {float(row):.12g} sobre .002.")
        print(f"  Máximo {float(highest):.12f}. Pasivo {float(passive):.12f}, activo {float(active):.12f}.",flush=True)


def bulk_degrees(weights, a, weak, strong, maximum):
    """Coeficientes positivos por grado de la serie de árboles usada desde once datos."""
    low = [w/a**(7+i) for i,w in enumerate(weights)]
    high = [w/a**(7+3*i) for i,w in enumerate(weights)]
    d = [[Q(0)] for _ in range(2)]
    e = [[Q(0)] for _ in range(2)]
    x = [[],[]]
    z = [[],[]]
    ex = [[Q(1)], [Q(1)]]
    ez = [[Q(1)], [Q(1)]]
    for degree in range(1, maximum+1):
        n = degree-1
        for i in (0,1):
            x[i].append(sum(k*row[n] for k,row in zip(weak[i],d)))
            z[i].append(sum(k*row[n] for k,row in zip(strong[i],e)))
            if n:
                for series, exponential in ((x[i],ex[i]),(z[i],ez[i])):
                    exponential.append(sum(j*series[j]*exponential[n-j] for j in range(1,n+1))/n)
            mixed = sum(ex[i][n-j]*ez[i][j] for j in range(1,n+1))
            d[i].append(low[i]*ex[i][n]+high[i]*mixed)
            e[i].append(high[i]*(ex[i][n]+mixed))
    return d


def check_common_outputs():
    # Cada una de las r celdas iniciales de la franja sigue encendida en la
    # primera salida. Las otras k-1 huellas que emiten apagan al menos dos.
    def outside(width, k):
        bounds = []
        for r in range(k,width):
            first = max(2, ceil(Q(r,width-r)))
            bounds.append(max(k-1, 2*(k-1)+first-r))
        return min(bounds)
    assert (outside(3,2),outside(4,2),outside(4,3)) == (2,2,4)
    u, z = Q(".00000101"), Q(".0012")
    all_groups = 64*1297*u**2/(1-8*3*1297*u)
    vv = u*comb(207,2)*z**2+all_groups
    cc = u*(3*comb(248,2)*z**2+comb(248,4)*z**4)+all_groups
    assert vv < Q("1.185e-7") and cc < Q("2.214e-7")
    print("Subconjuntos de dos o más salidas autónomas desde un dato fijado en G3: VV <1.185e-7, CC <2.214e-7.")
    for i,j in ((0,0),(0,1),(1,1)):
        for d1,d2 in displacements(i,j):
            for target in range(3):
                for axis, distance in enumerate((d1,d2)):
                    tree = SeedTree(i,distance)
                    for layer in (0,1):
                        a = RADII[i][target][axis]+2*layer
                        b = RADII[j][target][axis]+2*layer
                        direct = tree.ball(tree.left,a) & tree.ball(tree.right,b)
                        assert ball_intersection(i,distance,a,b) == len(direct)
        mins = [min(overlap_layers(i,j,*k,t)[0] for k in displacements(i,j)) for t in range(3)]
        assert mins == [1792,1197,90]
    # Las columnas corresponden a todos los datos próximos y los centros
    # que leen la franja del dato de raíz.
    for i in (0,1):
        tree = SeedTree(i,0)
        near = [len(tree.ball(0,r))*len(tree.ball(0,s)) for r,s in DISTANCES[i]]
        reader_radii = (1,2) if i==0 else (2,1)
        near.append(len(tree.ball(0,reader_radii[0]))*len(tree.ball(0,reader_radii[1])))
        assert tuple(near)==DIRECT[i]


def check_g3_obstruction(data):
    """Condición necesaria del mayorante efectivamente usado, con b y theta fijos."""
    halos, series, _, _ = data
    v11 = prior["previous"]("verificar_H5A_3_4_contacto_oculto.py")["v11"]
    far = v11["prop"]
    b, theta, p, q = Q(".001145"), Q(".00115"), Q("1e-6"), Q(".001")
    def lower_exp(x):
        return sum(x**k/factorial(k) for k in range(8))
    def upper_log(x):
        return x+x*x/(2*(1-x))
    collision = Q(89,2)*upper_log(theta*theta)
    guard = upper_log(theta)-collision
    other = b+collision-guard
    mass = (Q("1.84e-6"),Q("2.14e-6"))
    for i,f in enumerate((504,628)):
        assert mass[i] < p*lower_exp(Q(".0191")+f*guard+halos[0][i]*other)
    for a, stop, minimum in ((Q(".569"),11,Q("1.0035")),(Q(".575"),14,Q("1.00138"))):
        levels=bulk_degrees(mass,a,far["W"],far["B"],stop)
        for i in (0,1):
            ten = sum(n*mass[0]**nv*mass[1]**nc/a**nu for (nv,nc,nu),n in series[i].items())
            assert sum(levels[i][:11])==ten
        h = [sum(row[11:stop+1]) for row in levels]
        I = 588*(40*h[0]+27*h[1])
        phi = lower_exp(89*b+I)
        ratio = q/(1-2*q)*lower_exp(Q(".015"))*(1+b)*(1+phi*(phi-1))/b
        assert ratio > minimum
        print(f"G3, a≤{float(a)}, grados 11 a {stop}. Cociente pasivo necesario >{float(minimum):g}.")
        if stop==11:
            ball = v11["previous"]("medidas_por_historia")["ball"]
            links = far["kernels"](ball,a**6,Q(".673"))[0]
            S = min(sum(edge[0]*v for edge,v in zip(row,h)) for row in links)
            marked_only = sum((11*S)**k/factorial(k) for k in range(1,13))
            assert S>Q(".32424") and marked_only > Q("34.39")
            print("Misma cola con al menos un hijo y enlaces originales, s=.673: factor de grado once >34.39.")


def main():
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para comprobar las desigualdades.")
    check_common_outputs()
    data, circuit = previous_tables()
    verify_rows(data, circuit)
    check_g3_obstruction(data)


if __name__ == "__main__":
    main()
