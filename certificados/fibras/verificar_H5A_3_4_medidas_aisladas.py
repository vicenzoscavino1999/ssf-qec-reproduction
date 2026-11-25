"""H5-A. Medidas aisladas y grupos de datos. Python 3.10+, biblioteca estándar."""

from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations, product
from math import ceil, floor, factorial
from pathlib import Path
from runpy import run_path

W = ((486, 324), (576, 384))
B = ((154, 117), (208, 129))
NEIGHBORS = ((640, 441), (784, 513))
ROWS = tuple(7 << (3*i) for i in range(4))
COLS = tuple(sum(1 << (3*i+j) for i in range(4)) for j in range(3))


def guard_cost(representatives, budgets):
    # Bolas de radio dos en los factores del rectángulo 4×3 del recubrimiento.
    rows = [set(range(4)) | set(range(4+6*i, 10+6*i)) for i in range(4)]
    cols = [set(range(3)) | set(range(3+6*j, 9+6*j)) for j in range(3)]
    cells = [{(r,c): 1+int(r == i) for r in rows[i] for c in cols[j]}
             for i,j in product(range(4), range(3))]

    @lru_cache(None)
    def cover(mask):
        votes = Counter()
        for i, cell in enumerate(cells):
            if mask >> i & 1:
                votes.update(cell)
        return sum(min(n,2) for n in votes.values())

    for bits in range(1,128):
        phi = 0
        for j, stripe in enumerate(ROWS+COLS):
            if bits >> j & 1:
                phi ^= stripe
        if phi and 2*phi.bit_count() >= 3*bits.bit_count():
            assert phi.bit_count() >= 4 or phi in ROWS
            for i,j in combinations(range(12),2):
                if (phi >> i & 1) and (phi >> j & 1) and i//3 != j//3:
                    assert 100*(4-phi.bit_count()) < phi.bit_count()
    costs = {(0, phi.bit_count(), lit.bit_count(), cover(phi)-cover(lit))
             for lit,phi in representatives}
    ranks = budgets(costs,lambda _,increment: increment)
    assert ranks == ((138,108,0),(138,138,64,0))
    initial = cover(ROWS[0]), cover(COLS[0])
    assert initial == (258,288) and cover(1) == 99
    return tuple(initial[i]+sum(ranks[i]) for i in range(2))


# Clase 0 fuera de G, 1 débil, 2 fuerte restante,
# 3 igual coordenada fija sin solapamiento, 4 huellas con intersección.
def positions(i,j):
    if i != j:
        return {1:324,2:117,4:9} if i == 0 else {1:576,2:208,4:16}
    return {1:486,2:90,3:64,4:9} if i == 0 else {1:384,2:72,3:57,4:8}


def triple_budget(types, edges, profiles, stable, *, with_debits=False):
    pairs = ((0,1),(0,2),(1,2))
    relation = dict(zip(pairs,edges))
    widths = tuple(3+i for i in types)
    states = list(product(*(range(n+1) for n in widths)))
    index = {r:i for i,r in enumerate(states)}
    rank = [max(12*k//n for k,n in zip(r,widths)) for r in states]
    allowed = {k:{} for k in set(rank)}

    for lit,h,v,options in profiles:
        regular, perfect = set(), set()

        def assign(i, used, previous, counts):
            if i == 3:
                if used:
                    regular.add(index[counts])
                    if h == v and stable[lit & ~used]:
                        perfect.add(index[counts])
                return
            for mask, num, blocked in options[types[i]]:
                if mask & used:
                    continue
                axis = (ROWS,COLS)[types[i]]
                coordinate = next((j for j,s in enumerate(axis) if s & mask), -1)
                for j,(other,coord,hit) in enumerate(previous):
                    if not (mask and other):
                        continue
                    cl = relation[j,i]
                    if cl == 0 or (cl == 1 and mask & hit):
                        break
                    if types[j] == types[i] and ((coord == coordinate) != (cl >= 3)):
                        break
                else:
                    assign(i+1, used|mask, previous+((mask,coordinate,blocked),), counts+(num,))

        assign(0,0,(),())
        for k, moves in allowed.items():
            if k*v <= 12*h:
                for move in perfect if k == 12 else regular:
                    cost = 2*h-v
                    moves[move] = min(moves.get(move,cost),cost)

    labels = list(range(3))
    for (i,j),cl in relation.items():
        if cl == 4:
            a,b = labels[i],labels[j]
            labels = [a if c == b else c for c in labels]
    blocks = [[i for i,c in enumerate(labels) if c == b] for b in set(labels)]
    discount = [0]
    best = 0
    for i,r in enumerate(states):
        if i:
            moves = allowed[rank[i]]
            possible = (index[a] for a in product(*(range(z+1) for z in r)))
            discount.append(min(moves[m]+discount[i-m] for m in possible if m in moves))
        balance = 0
        for group in blocks:
            n = sum(widths[j] for j in group)
            s = sum(r[j] for j in group)
            if len(group) == 1:
                balance += 2*s-n
            elif s > n-2:
                break
            else:
                balance += s-(n-s)%2
        else:
            best = max(best,balance-discount[i])
    if with_debits:
        full_moves = {states[j]: cost for j,cost in allowed[12].items()}
        return best, dict(zip(states,discount)), full_moves
    return best


def count_triples(profiles, stable, rooted_trees):
    rows = [Counter(), Counter()]
    pairs = ((0,1),(0,2),(1,2))
    trees = list(rooted_trees(3))
    cached = lru_cache(None)(lambda t,e: triple_budget(t,e,profiles,stable))
    for types in product(range(2),repeat=3):
        for edges in product(*[(0,)+tuple(positions(types[i],types[j])) for i,j in pairs]):
            relation = dict(zip(pairs,edges))
            bounds = []
            for tree in trees:
                count = 1
                for i,j in tree:
                    cl = relation[min(i,j),max(i,j)]
                    if cl == 0:
                        break
                    count *= positions(types[i],types[j])[cl]
                else:
                    bounds.append(count)
            if bounds:
                key = min((tuple(types[j] for j in p),
                           tuple(relation[min(p[i],p[j]),max(p[i],p[j])] for i,j in pairs))
                          for p in permutations(range(3)))
                k = cached(*key)
                rows[types[0]][types.count(0),types.count(1),k] += Q(min(bounds),2)
    assert rows[1][0,3,4] == 52032 and rows[0][3,0,3] == 2048
    return rows


def group_series(limit):
    # Coeficientes por número de VV, de CC y por saldo en sextos.
    d,e = [{},{}],[{},{}]
    weak_forest = [{(0,0,0):Q(1)} for _ in range(2)]
    full_forest = [{(0,0,0):Q(1)} for _ in range(2)]
    weak_log, full_log = [{},{}],[{},{}]

    def add_root(i, forest, dest, power, sign=1):
        for (v,c,k),weight in forest.items():
            dest[v+1-i,c+i,k+power] += sign*weight

    def exponential(term, forest, degree):
        result = defaultdict(Q)
        for (a,b,k),val in term.items():
            for (c,d2,l),factor in forest.items():
                if a+b+c+d2 == degree:
                    result[a+c,b+d2,k+l] += (a+b)*val*factor/degree
        return result

    for n in range(1,limit+1):
        for i in (0,1):
            dn,en = defaultdict(Q),defaultdict(Q)
            weak = {k:v for k,v in weak_forest[i].items() if sum(k[:2]) == n-1}
            full = {k:v for k,v in full_forest[i].items() if sum(k[:2]) == n-1}
            add_root(i,weak,dn,7+i)
            add_root(i,weak,dn,7+3*i,-1)
            add_root(i,full,dn,7+3*i)
            add_root(i,full,en,7+3*i)
            assert all(x >= 0 for x in dn.values())
            d[i].update({k:v for k,v in dn.items() if v})
            e[i].update(en)
        keys = set(d[0]) | set(d[1]) | set(e[0]) | set(e[1])
        for i in (0,1):
            for key in keys:
                if sum(key[:2]) == n:
                    weak_log[i][key] = sum(W[i][j]*d[j].get(key,0) for j in (0,1))
                    full_log[i][key] = weak_log[i][key]+sum(B[i][j]*e[j].get(key,0) for j in (0,1))
            weak_forest[i].update(exponential(weak_log[i],weak_forest[i],n))
            full_forest[i].update(exponential(full_log[i],full_forest[i],n))
    return d

def pair_rows(Y):
    return ((576+64/Y+9/Y**2,441/Y+9/Y**2),
            (784/Y+16/Y**2,384+72/Y+65/Y**2))


def certify(tables, triples, series, halos, guard_sizes, geometric_rows, exp_upper):
    p,q,b,alpha,eta = map(Q,("2.6e-7",".001",".00116",".001",".0001"))
    a = Q(".540626883624")
    Y = a**6
    assert 0 < 36*Y < 1 and guard_sizes == (504,628)

    def upper(x):
        return Q(ceil(exp_upper(x)*10**24),10**24)

    def cols(A,x):
        return sum(k*v for k,v in zip(A,x))

    m = tuple(map(Q,("5.300793027156959e-7","6.326453178005993e-7")))
    u = tuple(map(Q,("3.6832451357961466e-7","4.623636566733828e-7")))
    h = tuple(map(Q,("6.351484662926536e-8","3.20512316990277e-7")))
    physical = tuple(map(Q,("1.3847600325766217e-8","8.65692678381971e-8")))
    d = tuple(map(Q,("4.724587973873814e-5",".0001222116205515124")))
    e = tuple(map(Q,("4.7245879738738135e-5",".00036841248946804854")))
    visible, rooted = Q("3.528397180234843e-7"),Q("1.3915309628676695e-5")
    theta = Q(".00116874036")
    g = 40*physical[0]+27*physical[1]+rooted
    R = max(cols(line,u) for line in NEIGHBORS)
    I = 280*u[0]+189*u[1]+588*g
    y = q*upper(Q(".015"))/(1-2*q)
    z = p*upper(Q(".015")+alpha)/(1-2*p)
    ell = 89*b+I
    eps = upper(ell)*(upper(ell)-1)
    assert eps < 1 and R < Q(".00052596") < alpha and max(u) < eta
    assert R*Q(6,13)**3 < Q(".00005172")
    prices = [y*(1+b)*(1+eps)/b]

    def neglog_upper(t):
        return t+t*t/(2*(1-t))

    collision = Q(89,2)*neglog_upper(theta*theta)
    log_zero = neglog_upper(theta)
    guard = log_zero-collision
    ordinary = b+2*collision-log_zero
    assert guard > 0 and ordinary > 0
    prices.append(y*(1-theta)*upper(99*guard+90*ordinary+2*I)/theta)
    for i in (0,1):
        correction = 2*(halos[1][i]*u[0]+halos[2][i]*u[1]+halos[3][i]*g)
        prices.append(z*upper(eta+2*R+guard*guard_sizes[i]+ordinary*halos[0][i]+correction)/m[i])

    dd,dg,gd,gg,pointed = geometric_rows(Y)
    w = tuple(m[i]*upper(cols(dd[i],h)+dg[i]*visible) for i in (0,1))
    pair = pair_rows(Y)
    initial = [(210,252),(640,784),(441,513),(1344,1596)]
    tree = tuple(mi*Q("1.003") for mi in m)
    tails = []
    for i in (0,1):
        # Comparación de la recurrencia de series con la enumeración etiquetada
        # heredada. En esos cuatro grados coinciden todos los coeficientes.
        reference = Counter()
        for (nv,nc,k,small),n in tables[i].items():
            reference[nv,nc,k] += n
        assert reference == {k:n for k,n in series[i].items() if sum(k[:2]) <= 4}

        small = Counter()
        for (nv,nc,k,short),n in tables[i].items():
            if nv+nc in (3,4):
                small[nv,nc] += n/Y**short
        alternative = Counter()
        for (nv,nc,k),n in triples[i].items():
            alternative[nv,nc] += n/Y**k
        assert Q(".08618") < alternative[1,2]/small[1,2] < Q(".08619")
        for key,new in alternative.items():
            small[key] = min(small[key],new)
        for (nv,nc,k),n in series[i].items():
            if nv+nc >= 5:
                small[nv,nc] += n/Y**(k//6)

        def first_ten(values):
            poly = sum(n*values[0]**nv*values[1]**nc for (nv,nc),n in small.items())
            return values[i]*(Y+cols(pair[i],values))+poly

        X,Z = cols(W[i],d),cols(B[i],e)
        low,high = w[i]/a**(7+i),w[i]/a**(7+3*i)
        full = upper(X)*(low+high*(upper(Z)-1))
        prices.extend((full/d[i],high*upper(X+Z)/e[i]))
        unrounded = sum(n*w[0]**nv*w[1]**nc/a**k for (nv,nc,k),n in series[i].items())
        tail = full-unrounded
        # Coeficiente de todos los árboles débiles de grado 11 del propio tipo.
        assert tail >= Q(11**9,factorial(10))*W[i][i]**10*low**11
        tails.append(tail)
        branches = first_ten(w)+tail
        connected = upper(cols(NEIGHBORS[i],tree))
        assert m[i]*connected <= tree[i]
        physical_row = branches-first_ten(m)+m[i]*(connected-1)
        factors = b,2*u[0],2*u[1],2*g
        simple = z*upper(eta+2*R+sum(v*s[i] for v,s in zip(factors,initial)))*(1+y)**(3+i)
        prices.extend((branches/h[i], physical_row/physical[i], (simple+physical[i])/u[i]))
    J = cols(gd,h)+gg*visible
    marked = cols(pointed,h+(visible,))
    prices.extend((Y*(upper(J)-1-J)/visible,Y*upper(J)*marked/rooted))
    for j,ratio in enumerate(prices):
        assert 0 < ratio < 1,(j,float(ratio))
    print(f"p=2.6e-7, q=.001. Fila uniforme {float(R):.12g} < .00052596.")
    print(f"Respuesta pasiva {float(prices[0]):.12f}; medidas activas {float(prices[1]):.12f}.")
    print(f"Mayor cociente {float(max(prices)):.12f}; colas desde once por raíz VV y CC "
          f"{float(tails[0]):.10g}, {float(tails[1]):.10g}.")


def joint_obstruction(geometric_rows):
    # Sólo se usan cotas inferiores de las filas que definen este mayorante.
    # Los extremos del intervalo de Y se eligen en la dirección favorable.
    p = Q("9e-7")
    def down(x):
        return Q(floor(x*10**18),10**18)
    def exp_low(x):
        x = min(x,Q(20))
        return sum(x**k/factorial(k) for k in range(9))

    m = [down(p*exp_low(Q(k,1000))) for k in (504,628)]
    grid = [Q(0)]+[Q(j,1000) for j in range(19,28)]+[Q(1,36)]
    steps = []
    for lo,hi in zip(grid,grid[1:]):
        # Una raíz sexta inferior por búsqueda entera.
        left,right = 1,30000
        while left < right:
            mid = (left+right+1)//2
            if hi*Q(mid,10000)**6 <= 1:
                left = mid
            else:
                right = mid-1
        base = Q(left,10000)
        assert hi*base**6 <= 1 < hi*(base+Q(1,10000))**6
        L = geometric_rows(lo)[0]
        pairs = pair_rows(hi)
        d,e,h = [Q(0)]*2,[Q(0)]*2,[Q(0)]*2
        for t in range(1,151):
            x = [down(m[i]*exp_low(sum(k*u for k,u in zip(L[i],h)))) for i in (0,1)]
            h = [down(x[i]*(lo+sum(k*u for k,u in zip(pairs[i],x)))) for i in (0,1)]
            X = [sum(k*u for k,u in zip(row,d)) for row in W]
            Z = [sum(k*u for k,u in zip(row,e)) for row in B]
            e = [down(x[i]*base**(7+3*i)*exp_low(X[i]+Z[i])) for i in (0,1)]
            d = [down(x[i]*exp_low(X[i])*(base**(7+i)+base**(7+3*i)*(exp_low(Z[i])-1)))
                 for i in (0,1)]
            if e[1] >= Q(".1"):
                break
        assert e[1] >= Q(".1")
        assert m[1]*base**10*exp_low(Q("12.9")) > Q(".1")
        steps.append(t)
    print("Obstrucción de las filas de parejas y grupos grandes desde p=9e-7 "
          f"para todo 0<Y<1/36. Iteraciones por intervalo {steps}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para comprobar las desigualdades.")
    folder = Path(__file__).resolve().parent
    def previous(name):
        return run_path(str(folder/f"verificar_H5A_3_4_{name}.py"))
    history,contacts,striped = map(previous,("medidas_por_historia","contacto_de_checks","franjas"))
    base = previous("medidas_y_un_dato")
    profiles,stable,representatives = contacts["geometry"]()
    _,caps = striped["cost_tables"](striped["light_profiles"](),previous("fondo_pasivo")["step_costs"],
                                    previous("pares_adyacentes")["stripe_profiles"])
    hist = contacts["positions"](profiles,stable,caps,striped["coefficients"](4))
    series = group_series(10)
    checktrees = striped["coefficients"](10)
    for i in (0,1):
        geometry_only = Counter()
        for (nv,nc,k),n in series[i].items():
            geometry_only[nv,nc] += n
        assert geometry_only == {key:n for key,n in checktrees[i].items() if n}
    certify(history["data_coefficients"](hist),count_triples(profiles,stable,contacts["rooted_trees"]),
            series,contacts["exposure"](representatives,base["budgets"]),
            guard_cost(representatives,base["budgets"]),history["geometric_rows"],base["exp_upper"])
    joint_obstruction(history["geometric_rows"])


