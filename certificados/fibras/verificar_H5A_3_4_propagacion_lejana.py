"""HGP (3,4). Presupuesto de pasos lejanos y certificado uniforme; Python 3.10+."""

from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from heapq import heappop, heappush
from itertools import combinations, combinations_with_replacement, permutations, product
from math import ceil, factorial, floor
from pathlib import Path
from runpy import run_path

ROWS = tuple(7 << (3*i) for i in range(4))
COLS = tuple(sum(1 << (3*i+j) for i in range(4)) for j in range(3))
N = ((640,441),(784,513))
W = ((486,324),(576,384))
B = ((154,117),(208,129))
LIMIT = 6


def dot(row, values):
    return sum(k*x for k,x in zip(row,values))


def previous(stem):
    return run_path(str(Path(__file__).with_name(f"verificar_H5A_3_4_{stem}.py")))


def winners(reps):
    def permute(mask, bits):
        return sum(bit for j,bit in enumerate(bits) if mask >> j & 1)
    found = set()
    for rows,cols in product(permutations(range(4)),permutations(range(3))):
        bits = tuple(1 << (3*rows[j//3]+cols[j%3]) for j in range(12))
        found.update((permute(h,bits),permute(p,bits)) for h,p in reps)
    assert len(found) == 771 and max((p ^ h).bit_count() for h,p in found) == 2
    return found


def projected_pair(first, second, max_product, local, stable):
    """Distancia mínima de descenso en los bits observados de las dos franjas."""
    sites = sorted(set(first)|set(second))
    index = {site:1 << i for i,site in enumerate(sites)}
    stripes = [sum(index[s] for s in stripe) for stripe in (first,second)]
    widths = [len(first),len(second)]
    def subsets(items, limit):
        return (p for k in range(min(limit,len(items))+1) for p in combinations(items,k))
    moves = {}
    for rows,cols in product(subsets(sorted({r for r,c in sites}),4),
                             subsets(sorted({c for r,c in sites}),3)):
        bits = [index.get((r,c),0) for r in rows+(None,)*(4-len(rows))
                for c in cols+(None,)*(3-len(cols))]
        covered = sum(bits)
        if not covered or (covered&stripes[0]).bit_count()*(covered&stripes[1]).bit_count() > max_product:
            continue
        raw = sum(1 << i for i,val in enumerate(bits) if val)
        @lru_cache(None)
        def project(mask):
            return sum(bit for i,bit in enumerate(bits) if mask >> i & 1)
        for lit,phi in local:
            hit = project(phi)
            if not hit:
                continue
            h,v = lit.bit_count(),phi.bit_count()
            key = hit,project(lit),12*h//v
            costs = moves.setdefault(key,[1000,1000])
            delta = 2*h-v
            costs[0] = min(costs[0],delta)
            if h == v and stable[phi & ~raw]:
                costs[1] = min(costs[1],delta)
    size = 1 << len(sites)
    graph = [[] for _ in range(2*size)]
    for state in range(size):
        rank = max(12*(state&stripe).bit_count()//n for stripe,n in zip(stripes,widths))
        if rank < 12:
            graph[size+state].append((state,0))
        for (hit,lit,top),(cost,perfect) in moves.items():
            if top >= rank and state & hit == lit:
                graph[state].append((state ^ hit,cost))
                if perfect < 1000:
                    graph[size+state].append((size+(state ^ hit),perfect))
    cap = 0
    for final in range(size):
        if any(2*(final&stripe).bit_count() > n for stripe,n in zip(stripes,widths)):
            continue
        start = final ^ stripes[0] ^ stripes[1]
        queue = [(0,size+start)]
        distance = {size+start:0}
        while queue:
            cost,state = heappop(queue)
            if cost != distance[state]:
                continue
            if state % size == final:
                cap = max(cap,start.bit_count()-final.bit_count()-cost)
                break
            for other,step in graph[state]:
                new = cost+step
                if new < distance.get(other,10**6):
                    distance[other] = new
                    heappush(queue,(new,other))
    return cap


def pair_counts(local,stable):
    column = [(j,0) for j in range(4)]
    disjoint = [(j,0) for j in range(4,8)]
    assert [projected_pair(column,disjoint,k,local,stable) for k in (1,3,4)] == [0,1,2]
    row = [(0,j) for j in range(3)]
    second = [(j,3) for j in range(1,5)]
    assert projected_pair(row,second,1,local,stable) == 0
    # Coeficientes constantes, de saldo uno y de saldo dos. Para la marca
    # física se verifican por separado las comparaciones de multiplicidad.
    return (((576,54,9),(0,279,9)),((0,496,16),(384,72,20)))


def cover_profiles(reps):
    rows = [set(range(4)) | set(range(4+6*i,10+6*i)) for i in range(4)]
    cols = [set(range(3)) | set(range(3+6*j,9+6*j)) for j in range(3)]
    sites = [{(r,c):1+(r == i) for r in rows[i] for c in cols[j]}
             for i,j in product(range(4),range(3))]
    @lru_cache(None)
    def cover(mask):
        votes = Counter()
        for i,cell in enumerate(sites):
            if mask >> i & 1:
                votes.update(cell)
        return sum(min(2,n) for n in votes.values())
    profiles = {(lit.bit_count(),phi.bit_count(),max((lit&r).bit_count() for r in ROWS),
                 max((lit&c).bit_count() for c in COLS),cover(phi)-cover(lit)) for lit,phi in reps}
    assert max(Q(inc,2*h-v) for h,v,_,_,inc in profiles) == 138
    return profiles


def group_budget(nv,nc,profiles,donors,guard=False):
    widths = (3,)*nv+(4,)*nc
    infinity = 10**6
    def canonical(x):
        return tuple(sorted(x[:nv]))+tuple(sorted(x[nv:]))
    @lru_cache(None)
    def debit(a,rank):
        total = sum(a)
        if sum(k > 0 for k in a) > donors:
            return infinity
        best = infinity
        for h,v,row,col,inc in profiles:
            if rank*v > 12*h or total > h:
                continue
            if rank == 12 and (v-total > 3 or 2*total < v):
                continue
            if max(a[:nv],default=0) > row or max(a[nv:],default=0) > col:
                continue
            delta = 2*h-v
            best = min(best,138*delta-inc if guard else delta)
        return best
    @lru_cache(None)
    def minimal(r):
        if not any(r):
            return 0
        rank = max(12*k//n for k,n in zip(r,widths))
        best = infinity
        for a in product(*(range(k+1) for k in r)):
            if not any(a):
                continue
            if any(j != nv and r[j] == r[j-1] and a[j] < a[j-1] for j in range(1,len(r))):
                continue
            pay = debit(canonical(a),rank)
            if pay < best:
                rest = canonical(tuple(k-x for k,x in zip(r,a)))
                best = min(best,pay+minimal(rest))
        return best
    best = 0
    for vv,cc in product(combinations_with_replacement(range(4),nv),
                         combinations_with_replacement(range(5),nc)):
        r = vv+cc
        s = sum(r)
        balance = s-(sum(widths)-s)%2
        best = max(best,(138 if guard else 1)*balance-minimal(r))
    return best


def crowded_coefficients(trees):
    """Árbol con j lectores distintos marcados en un rectángulo común."""
    def multiply(a,b):
        out = Counter()
        for (v,c),x in a.items():
            for (w,d),y in b.items():
                if v+c+w+d <= LIMIT:
                    out[v+w,c+d] += x*y
        return out
    center = Counter()
    for i in range(2):
        for key,x in trees[i].items():
            center[key] += (40,27)[i]*x
    output = {}
    power = {(0,0):Q(1)}
    for j in range(2,LIMIT+1):
        power = multiply(power,center)
        if j == 2:
            continue
        level = [{k:(30,36)[i]*v/factorial(j-1) for k,v in multiply(trees[i],power).items()}
                 for i in range(2)]
        total = [Counter(poly) for poly in level]
        for step in range(LIMIT-j):
            mix = [Counter(),Counter()]
            for i,t in product(range(2),repeat=2):
                for key,x in level[t].items():
                    mix[i][key] += N[i][t]*x
            level = [multiply(trees[i],mix[i]) for i in range(2)]
            for i in range(2):
                total[i].update(level[i])
        output[j] = total
    return output


def kernels(ball,Y,s):
    outer,rho = Y/s**2,Y*s
    assert 0 < outer < 1 and 0 < 36*rho < 1 and 0 < s < 1
    def kernel(v1,r1,v2,r2,ends):
        def size(l):
            return ball(v1,r1+2*l)*ball(v2,r2+2*l)
        c36 = Q(size(2)-7*size(1)+6*size(0),1050)
        c6 = Q(size(1)-size(0)-35*c36,5)
        c1 = size(0)-c6-c36
        assert c6 <= 0 and c6+c1 <= 0
        assert all(size(l) == c1+c6*6**l+c36*36**l for l in range(9))
        def weight(l):
            return outer**min(l,ends)*rho**max(0,l-ends)
        exact = size(0)+sum((size(l)-size(l-1))*weight(l) for l in range(1,ends+1))
        exact += outer**ends*(c6*5*6**ends*rho/(1-6*rho)+c36*35*36**ends*rho/(1-36*rho))
        cut = ends+14
        marked = max([size(l)*weight(l) for l in range(cut)]+[c36*36**cut*weight(cut)])
        return ceil(exact),ceil(marked)
    dd = ((kernel(1,4,1,6,2),kernel(1,5,1,5,2)),(kernel(0,5,0,5,2),kernel(0,6,0,4,2)))
    dg = (kernel(1,3,1,4,1),kernel(0,4,0,3,1))
    gd = (kernel(0,3,1,4,1),kernel(0,4,1,3,1))
    return dd,dg,gd,kernel(0,2,1,2,0)


def anchor(ball,i,k):
    if k < 1:
        return 0
    steps = 1
    extra = 0
    while (3*steps+1)//2 <= k:
        steps = (3*steps+1)//2
        extra += 1
    return ball(0,(3,4)[i]+2*extra)*ball(1,(4,3)[i]+2*extra)


def group_polynomials(Y,ball,tables,triples,series,pairs,geometry,trees):
    caps,credit = {},{}
    for n in range(1,LIMIT+1):
        for nv in range(n+1):
            nc = n-nv
            saved = max(0,246*nv+340*nc-group_budget(nv,nc,geometry,n,guard=True))
            # Techo racional monótono de exp(-.001·saved).
            k = Q(saved,1000)
            credit[nv,nc] = 1/sum(k**j/factorial(j) for j in range(6))
            if n >= 3:
                caps[nv,nc] = {j:group_budget(nv,nc,geometry,j) for j in range(2,n+1)}
                seq = list(caps[nv,nc].values())
                assert seq == sorted(seq)
    crowded = crowded_coefficients(trees)
    poly,marked = [Counter(),Counter()],[Counter(),Counter()]
    for i in range(2):
        cost = (lambda k:Y**(-k),lambda k:anchor(ball,i,k))
        target = (poly[i],marked[i])
        for mode,price in enumerate(cost):
            part = target[mode]
            part[1-i,i] = (Y if mode == 0 else 0)*credit[1-i,i]
            # Las cotas por caminos también valen con el precio físico de la marca.
            assert 2*price(0) <= price(1) and 2*price(1) <= price(2)
            for j in range(2):
                key = 2-i-j,i+j
                part[key] = dot(pairs[i][j],[price(k) for k in range(3)])*credit[key]
            low = Counter()
            for (nv,nc,total,k),mul in tables[i].items():
                if nv+nc in (3,4):
                    low[nv,nc] += mul*price(k)
            alt = Counter()
            for (nv,nc,k),mul in triples[i].items():
                alt[nv,nc] += mul*price(k)
            for key,val in alt.items():
                low[key] = min(low[key],val)
            sparse = Counter()
            for (nv,nc,nu),mul in series[i].items():
                n = nv+nc
                if n >= 5:
                    low[nv,nc] += mul*price(nu//6)
                if 3 <= n <= LIMIT:
                    sparse[nv,nc] += mul*price(min(nu//6,caps[nv,nc][2]))
            for key,val in low.items():
                if key in sparse:
                    bound = caps[key]
                    extra = sum(crowded[j][i].get(key,0)*(price(bound[j])-price(bound[j-1]))
                                for j in range(3,sum(key)+1))
                    val = min(val,sparse[key]+extra)
                part[key] = val*credit.get(key,1)
    return poly,marked


def certify(a,s,halos,series,poly,pure,kernel,exp_upper):
    p,q,b,alpha,eta = map(Q,("3.2e-7",".001",".001155",".001",".0001"))
    Y,outer = a**6,a**6/s**2
    m = tuple(map(Q,("6.32379443599e-7","7.48130074773e-7")))
    u = tuple(map(Q,("4.32327960619e-7","4.71195023121e-7")))
    h = tuple(map(Q,("5.92784029774e-8","2.88927830255e-7")))
    physical = tuple(map(Q,("2.71677227720e-9","1.70822883877e-8")))
    d = tuple(map(Q,("5.21362317064e-5",".000132561958282")))
    e = tuple(map(Q,("5.21362317064e-5",".000387140922708")))
    visible,rooted = Q("4.00577118323e-7"),Q("1.08244980135e-5")
    theta = Q(".001160695045")
    def upper(x):
        return Q(ceil(10**24*exp_upper(x)),10**24)
    dd,dg,gd,gg = kernel
    Ldd = tuple(tuple(v[0] for v in row) for row in dd)
    g = 40*physical[0]+27*physical[1]+rooted
    R = max(dot(row,u) for row in N)
    I = 280*u[0]+189*u[1]+588*g
    y = q*upper(Q(".015"))/(1-2*q)
    z = p*upper(Q(".015")+alpha)/(1-2*p)
    ell = 89*b+I
    eps = upper(ell)*(upper(ell)-1)
    assert eps < 1 and max(u) < eta and R < Q(".00058067") < alpha
    assert R*Q(13,16)**3 < Q(".00031147")
    ratios = [y*(1+b)*(1+eps)/b]
    def log_bound(t):
        return t+t*t/(2*(1-t))
    c = Q(89,2)*log_bound(theta**2)
    guard,ordinary = log_bound(theta)-c,b+2*c-log_bound(theta)
    assert guard >= Q(".001") and ordinary >= 0
    ratios.append(y*(1-theta)*upper(99*guard+90*ordinary+2*I)/theta)
    for i in range(2):
        interaction = 2*dot([r[i] for r in halos[1:]],(u[0],u[1],g))
        ratios.append(z*upper(eta+2*R+guard*(504,628)[i]+ordinary*halos[0][i]+interaction)/m[i])
    w = tuple(m[i]*upper(dot(Ldd[i],h)+dg[i][0]*visible) for i in range(2))
    depart = tuple(dot([v[1] for v in dd[i]],h)+dg[i][1]*visible for i in range(2))
    initial = [(210,252),(640,784),(441,513),(1344,1596)]
    trees = tuple(mi*Q("1.003") for mi in m)
    tails,removed = [],[]
    def evaluate(tab,values):
        return sum(k*values[0]**nv*values[1]**nc for (nv,nc),k in tab.items())
    for i in range(2):
        X,Z = dot(W[i],d),dot(B[i],e)
        low,high = w[i]/a**(7+i),w[i]/a**(7+3*i)
        full = upper(X)*(low+high*(upper(Z)-1))
        ratios.extend((full/d[i],high*upper(X+Z)/e[i]))
        base_ten = sum(n*w[0]**nv*w[1]**nc/a**nu for (nv,nc,nu),n in series[i].items())
        tail = full-base_ten
        assert tail >= Q(11**9,factorial(10))*W[i][i]**10*low**11
        tails.append(tail)
        branches = evaluate(poly[i],w)+tail
        removed.append(evaluate(poly[i],m))
        pointed = sum(k*w[0]**nv*w[1]**nc*(nv*depart[0]+nc*depart[1])
                      for (nv,nc),k in poly[i].items())
        connected = upper(dot(N[i],trees))
        assert m[i]*connected <= trees[i]
        fixed = min(branches-removed[i],pointed+tail)+m[i]*(connected-1)
        factors = b,2*u[0],2*u[1],2*g
        simple = z*upper(eta+2*R+dot(factors,[row[i] for row in initial]))*(1+y)**(3+i)
        ratios.extend((branches/h[i],fixed/physical[i],(simple+physical[i])/u[i]))
    J = dot([k[0] for k in gd],h)+gg[0]*visible
    one = dot([k[1] for k in gd],h)+gg[1]*visible
    without_exterior = dot([k[1] for k in gd],removed)
    direct = sum(evaluate(tab,m) for tab in pure)
    root = outer*(upper(J)*one-without_exterior)+direct
    ratios.extend((outer*(upper(J)-1-J)/visible,root/rooted))
    for j,ratio in enumerate(ratios):
        assert 0 < ratio < 1,(j,float(ratio))
    assert max(ratios) < Q(".99999981")
    assert tails[0] < Q("8.844e-11") and tails[1] < Q("5.403e-10")
    print(f"p=3.2e-7, q=.001. Fila uniforme {float(R):.12g} < .00058067.")
    print(f"Respuesta pasiva {float(ratios[0]):.12f}, medidas activas {float(ratios[1]):.12f}, "
          f"mayor cociente {float(max(ratios)):.12f}.")
    print(f"Restos desde once datos por raíz VV y CC {float(tails[0]):.10g}, {float(tails[1]):.10g}.")
    print(f"Raíz lejana {float(root):.12g}. Tres o más componentes en p=2.6e-7 < .00031147.")


def test_target(a,kernel,poly):
    """Condición necesaria en G3 para los pesos geométricos y descuentos de esta nota."""
    def down(x):
        return Q(floor(x*10**18),10**18)
    def exp_low(x):
        x = min(x,Q(20))
        return sum(x**j/factorial(j) for j in range(9))
    dd = [[val[0] for val in row] for row in kernel[0]]
    p = Q("1e-6")
    m = [down(p*exp_low(Q(k,1000))) for k in (504,628)]
    d,e,h = [Q(0)]*2,[Q(0)]*2,[Q(0)]*2
    for step in range(1,151):
        w = [down(m[i]*exp_low(dot(dd[i],h))) for i in range(2)]
        h = [down(sum(k*w[0]**nv*w[1]**nc for (nv,nc),k in tab.items() if nv+nc <= 2))
             for tab in poly]
        X,Z = [dot(row,d) for row in W],[dot(row,e) for row in B]
        d = [down(w[i]*exp_low(X[i])*(a**-(7+i)+a**-(7+3*i)*(exp_low(Z[i])-1))) for i in range(2)]
        e = [down(w[i]*a**-(7+3*i)*exp_low(X[i]+Z[i])) for i in range(2)]
        if e[1] >= Q(".1"):
            break
    assert e[1] >= Q(".1") and m[1]*a**-10*exp_low(Q("12.9")) > Q(".1")
    print(f"G3 con estos pesos geométricos: las filas de ramas son incompatibles "
          f"tras {step} iteraciones inferiores, incluso con el presupuesto de fila libre.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para verificar las desigualdades.")
    contacts,striped,base,history,old = map(previous,("contacto_de_checks","franjas","medidas_y_un_dato",
                                                        "medidas_por_historia","medidas_aisladas"))
    profiles,stable,reps = contacts["geometry"]()
    assert max(k.bit_count() for k,ok in enumerate(stable) if ok) == 3
    pairs = pair_counts(winners(reps),stable)
    _,prior_caps = striped["cost_tables"](striped["light_profiles"](),previous("fondo_pasivo")["step_costs"],
                                          previous("pares_adyacentes")["stripe_profiles"])
    hist = contacts["positions"](profiles,stable,prior_caps,striped["coefficients"](4))
    tables = history["data_coefficients"](hist)
    triples = old["count_triples"](profiles,stable,contacts["rooted_trees"])
    series = old["group_series"](10)
    a,s = Q(".5488178589126"),Q(".86")
    areas = kernels(history["ball"],a**6,s)
    poly,pure = group_polynomials(a**6,history["ball"],tables,triples,series,pairs,
                                  cover_profiles(reps),striped["coefficients"](LIMIT))
    certify(a,s,contacts["exposure"](reps,base["budgets"]),series,poly,pure,areas,base["exp_upper"])
    test_target(a,areas,poly)
