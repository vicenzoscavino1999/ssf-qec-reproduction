"""Dos conmutaciones del mismo check y distancia conjunta en HGP (3,4)."""
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement, product
from math import factorial
from pathlib import Path
from runpy import run_path


v12 = run_path(str(Path(__file__).with_name("verificar_H5A_3_4_episodios_compartidos.py")))
v11 = v12["v11"]
v10 = v11["prior"]
dot, exp_upper = v11["dot"], v11["upper"]


def nonnegative_sqrt2(a, b):
    """Decide a + b sqrt(2) >= 0 usando enteros."""
    if a < 0:
        return b > 0 and 2*b*b >= a*a
    return b >= 0 or a*a >= 2*b*b


def check_hidden_contacts():
    _, _, reps = v11["previous"]("contacto_de_checks")["geometry"]()
    winners = v11["prop"]["winners"](reps)
    distance = [[int(x//3 != y//3) + int(x%3 != y%3) for y in range(12)] for x in range(12)]
    tight = False
    for lit, footprint in winners:
        start = [x for x in range(12) if lit & (1 << x)]
        births = [x for x in range(12) if (lit ^ footprint) & (1 << x)]
        for heights in product(range(3), repeat=len(births)):
            if any(abs(i-j)>distance[x][y] for x,i in zip(births,heights) for y,j in zip(births,heights)):
                continue
            floor = [max([0] + [h-distance[x][y] for y,h in zip(births,heights)]) for x in start]
            diff = [floor.count(k)-heights.count(k) for k in range(3)]
            assert nonnegative_sqrt2(diff[0]+2*diff[2], diff[1])
            tight |= diff == [2, 0, -1]
    assert len(winners) == 771 and tight
    print("Distancia conjunta: 771 ganadores, factor sqrt(2) por paso de un eje.", flush=True)


def check_marked_trees():
    checked = 0
    n = 6
    for parents in product(*(range(i) for i in range(1,n))):
        neighbors = [set() for _ in range(n)]
        for i,p in enumerate(parents,1):
            neighbors[i].add(p)
            neighbors[p].add(i)
        for size in (2,3,4):
            for terminals in combinations(range(n),size):
                core = set(range(n))
                while True:
                    leaves = {i for i in core.difference(terminals) if len(core & neighbors[i])<=1}
                    if not leaves:
                        break
                    core -= leaves
                routes = {i:() for i in core}
                while len(routes)<n:
                    for i in set(range(n))-routes.keys():
                        known = neighbors[i] & routes.keys()
                        if known:
                            p, = known
                            routes[i] = (i,)+routes[p]
                for a,b in combinations_with_replacement(set(range(n)).difference(terminals),2):
                    if a!=b and a in core and b in core:
                        assert min(len(neighbors[i]&core) for i in (a,b))>=2
                    elif a!=b:
                        first,second = sorted((a,b),key=lambda i:len(routes[i]))
                        assert second not in core.union(routes[first])
                    else:
                        assert (core.union(routes[a]) & {a}) == {a}
                    checked += 1
    print(f"Elección de las dos conmutaciones comprobada en {checked} marcas de árboles.", flush=True)


def root_totals(kernel, witness, regular, group=Q(1)):
    """Costes de los hijos de una raíz generador en la supersolución congelada."""
    _, _, gd, gg = kernel
    h = tuple(map(Q, witness["h"]))
    m = tuple(map(Q, witness["m"]))
    visible = Q(witness["visible"]) if "visible" in witness else Q(witness["outer"][1])
    exp_j = exp_upper(dot([link[0] for link in gd], h) + gg[0]*visible)
    marked = [link[1] for link in gd]
    one = dot(marked, h) + gg[1]*visible
    removed = [group*v11["effective"]["evaluate"](p,m) for p in regular]
    base = dot(marked, removed)
    assert one > base >= 0
    return exp_j*one-base, (exp_j-1)*one


def corrected_kernel(kernel, witness, regular, outside):
    raw, at_least_two = root_totals(kernel,witness,regular)
    # División sólo en esta supersolución. No es un volumen geométrico de enlace.
    factor = (outside*raw + at_least_two/2)/raw
    assert 0 < factor < 1
    dd,dg,gd,gg = kernel
    return dd,dg,tuple((volume,marked*factor) for volume,marked in gd),(gg[0],gg[1]*factor)


def single_tail(m):
    largest = max(m)
    ratio = exp_upper(Q(1))*1297*largest*Q(12,11)**9
    assert ratio < 1
    return 2*2000*11**9*Q(11**9,factorial(10))*largest*(1297*largest)**10/(1-ratio)


RECTANGLE = dict(name="Contacto lejano marcado dos veces", a=".555",s=".7984",
    rates=("3.75e-7",".001145",".004",".000671","1.01"),
    m=("7.043188734393e-7","8.205884160728e-7"),u=("4.961693700309e-7","5.491758377318e-7"),
    h=("4.448763414933e-8","3.131839159086e-7"),physical=("4.403715698574e-9","3.198901990376e-8"),
    outer=(".001147193821302","7.315596589479e-7","1.047649197369e-6"),
    d=("5.645577281966e-5",".00014132235091"),e=("5.645577281966e-5",".0004029237191351"))


def check_rectangle(tables, discounts):
    w = RECTANGLE
    ball,halos,series,_,hist,_ = tables
    a,s = Q(w["a"]),Q(w["s"])
    y = a**6
    reg,pure = v12["prices"](y,Q(1),Q(w["outer"][0]),ball,discounts,hist)
    pure[0][0,0] = single_tail(tuple(map(Q,w["m"])))
    original = v11["prop"]["kernels"](ball,y,s)
    kernel = corrected_kernel(original,w,reg,y/s**2)
    d,e = (tuple(map(Q,w[t])) for t in ("d","e"))
    v11["effective"]["certify"](y,s,halos,reg,pure,kernel,
        lambda x,i:v11["effective"]["full_tail"](series,a,x,i,d,e),w)


def check_g3(tables,discounts):
    ball,_,_,_,hist,trees = tables
    ct = v11["previous"]("contacto_de_checks")
    _,stable,reps = ct["geometry"]()
    effective = v11["effective"]
    winners = v11["prop"]["winners"](reps)
    pairs = effective["pair_classes"](winners,stable,effective["effective_budgets"](winners,stable))
    w = v10["WITNESS"]
    y,s,b = Q(43,2000),Q(".9999"),Q(".00115")
    reg,pure = v10["sector"](y,ball,hist,pairs)
    m,u,H = (tuple(map(Q,w[t])) for t in ("m","u","physical"))
    d = tuple(map(Q,w["d"]))
    v10["old"]["certify"](w,y,s,ball,tables[1],reg,pure,
        lambda x,i:v10["tail"](trees(10),y,x,i,d),scale=1,group=Q(1),
        single_group_factor=1/y,tree_factor=Q("1.006"),b=b)
    outside = y/s**2
    kernel = v11["prop"]["kernels"](ball,y,s)
    raw,fork = root_totals(kernel,w,reg)
    direct = sum(effective["evaluate"](tab,m) for tab in pure) + single_tail(m)/y
    rooted = outside*(outside*raw+fork/2) + direct
    assert 0 < rooted < Q(w["rooted"])
    fixed_g = 40*H[0] + 27*H[1] + rooted
    factors,_ = v10["marked_sensitivity"](b,v11["Q_RATE"],u,H,fixed_g)
    fixed = [H[i]+(u[i]-H[i])*factors[i] for i in (0,1)]
    broad = max(dot(row,fixed) for row in v11["N"])/(v10["BONUS"]*y)
    # Las otras tres cotas disjuntas de v10 y las dos normas históricas de v12.
    other = Q(".001723672")+Q(".000111820")+Q(".000156022")
    sparse = v12["projected"](v12["SPARSE"],tables,discounts)
    angular = v12["projected"](v12["AXES"],tables,discounts)
    floor = v12["AXES"]["group"]**4 * Q(v12["AXES"]["xi"])**2
    base = other+broad+sparse+angular/floor
    plus = base+angular
    print(f"Raíz física G3 para la marca de v10 {float(rooted):.12g}. "
          f"Fila con L>=6Y {float(broad):.12g}.")
    print(f"Subtotal de v13 {float(base):.12g}. Con la misma clase de v14 "
          f"{float(plus):.12g}. Margen {float(Q('.004')-plus):.12g}.")
    assert broad<Q(".000549140") and base<Q(".003046603") and plus<Q(".003214411")
    assert base+Q(".00095339")<Q(".004") and plus+Q(".00078558")<Q(".004")
    assert other+broad+sparse+angular/Q(".120")<Q(".004")
    return broad,base,plus


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para verificar las desigualdades.")
    check_hidden_contacts()
    check_marked_trees()
    tables = v11["data_tables"]()
    discounts = v11["geometry_discounts"](tables[3])
    check_rectangle(tables,discounts)
    check_g3(tables,discounts)
