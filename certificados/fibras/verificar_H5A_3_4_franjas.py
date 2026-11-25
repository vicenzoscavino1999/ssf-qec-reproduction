"""HGP (3,4), rectángulo (1e-8,.001). Python 3.10+, biblioteca estándar."""

from fractions import Fraction as Q
from itertools import product
from math import ceil
from pathlib import Path
from runpy import run_path


A = ((640, 441), (784, 513))
FAR_COST = 11


def light_profiles():
    rows = tuple(7 << (3*i) for i in range(4))
    cols = tuple(sum(1 << (3*i+j) for i in range(4)) for j in range(3))
    stripes = rows+cols
    cuts = []
    for subset in range(1, 128):
        phi = 0
        for i, stripe in enumerate(stripes):
            if subset & (1 << i):
                phi ^= stripe
        size = phi.bit_count()
        if size and 2*size >= 3*subset.bit_count():
            assert size > 4 or subset.bit_count() == 1
            cuts.append((phi, size))
    assert len(cuts) == 101
    patterns = set()
    for state in range(4096):
        best_h, best_v, winners = -1, 1, []
        for phi, v in cuts:
            h = (state & phi).bit_count()
            cmp = h*best_v-best_h*v
            if cmp > 0:
                best_h, best_v, winners = h, v, []
            if cmp >= 0:
                winners.append((phi & state, phi))
        if 100*(2*best_h-best_v) >= best_v:
            patterns.update(winners)
    assert len(patterns) == 771

    def counts(lit, axis):
        return tuple(sorted(((lit & stripe).bit_count() for stripe in axis), reverse=True))

    return sorted({(phi.bit_count(), lit.bit_count(), counts(lit,rows), counts(lit,cols))
                   for lit, phi in patterns})


def partitions(size):
    """Particiones por igualdad de coordenada, en orden de primera aparición."""
    def extend(prefix):
        if len(prefix) == size:
            yield prefix
            return
        for label in range(max(prefix, default=-1)+2):
            yield from extend(prefix+(label,))
    return list(extend(()))


def loads(assignment, partition):
    totals = [0]*(max(partition, default=-1)+1)
    for count, group in zip(assignment, partition):
        totals[group] += count
    return tuple(sorted((v for v in totals if v), reverse=True))


def fits(demand, capacity):
    return len(demand) <= len(capacity) and all(x <= y for x,y in zip(demand,capacity))


def group_costs(nv, nc, profiles):
    widths = (3,)*nv+(4,)*nc
    states = list(product(*(range(n+1) for n in widths)))
    index = {r: i for i,r in enumerate(states)}
    totals = list(map(sum, states))
    # 12 es múltiplo de ambas anchuras. La comparación de rangos es entera.
    ranks = [max(12*k//n for k,n in zip(r,widths)) for r in states]
    subsets = [[index[a] for a in product(*(range(x+1) for x in r)) if any(a)]
               for r in states]
    result = {}
    for rows, cols in product(partitions(nv), partitions(nc)):
        # Tabla del menor descuento admisible para cada movimiento y rango.
        allowed = {rank: {} for rank in set(ranks)}
        for move,a in enumerate(states[1:], 1):
            row_demand = loads(a[:nv], rows)
            col_demand = loads(a[nv:], cols)
            for v,h,row_cap,col_cap in profiles:
                if totals[move] > h or not fits(row_demand,row_cap) or not fits(col_demand,col_cap):
                    continue
                for rank, costs in allowed.items():
                    if rank*v > 12*h or (rank == 12 and 2*totals[move] < v):
                        continue
                    cost = FAR_COST*(2*h-v)-1
                    if move not in costs or cost < costs[move]:
                        costs[move] = cost

        credit = [0]
        for i in range(1, len(states)):
            costs = allowed[ranks[i]]
            # El índice es de base mixta. Para a <= r, idx(r-a) = idx(r)-idx(a).
            credit.append(min(costs[k]+credit[i-k] for k in subsets[i] if k in costs))
        any_masks = max(FAR_COST*s-c for s,c in zip(totals,credit))
        disjoint = max(FAR_COST*(2*s-sum(widths))-c for s,c in zip(totals,credit))
        result[rows,cols] = (any_masks, disjoint)
    return result


def cost_tables(profiles, step_costs, stripe_profiles):
    bulk, old_pairs = step_costs(stripe_profiles(), FAR_COST)
    assert bulk == [14,Q(59,3)] and old_pairs == [(14,24),(14,24),(25,35)]
    vv = group_costs(2,0,profiles)
    cc = group_costs(0,2,profiles)
    v_different, v_same = vv[((0,1),())][1], vv[((0,0),())][1]
    c_different, c_same = cc[((),(0,1))][1], cc[((),(0,0))][1]
    c_overlap = max(cost[0] for cost in cc.values())
    assert (v_different,v_same,c_different,c_same,c_overlap) == (3,14,14,25,25)
    cap = {}
    for size in (3,4):
        for nc in range(size+1):
            costs = group_costs(size-nc,nc,profiles)
            cap[size-nc,nc] = max(any_masks for any_masks,_ in costs.values())
    assert [cap[3-c,c] for c in range(4)] == [37,37,47,48]
    assert [cap[4-c,c] for c in range(5)] == [49,49,60,61,71]
    # (número de posiciones, coste), primero VV y luego CC para raíz e hijo.
    pair = (
        (((640,v_different),(64,v_same),(9,old_pairs[0][1])),
         ((441,old_pairs[1][0]),(9,old_pairs[1][1]))),
        (((784,old_pairs[1][0]),(16,old_pairs[1][1])),
         ((513,c_different),(57,c_same),(8,c_overlap))),
    )
    return pair, cap


def coefficients(degree):
    """Árboles de datos enraizados T_i = z_i exp(sum_j A_ij T_j)."""
    trees = [{}, {}]
    forest = [{(0,0): Q(1)}, {(0,0): Q(1)}]
    for size in range(1,degree+1):
        for v in range(size+1):
            c = size-v
            trees[0][v,c] = forest[0].get((v-1,c),0)
            trees[1][v,c] = forest[1].get((v,c-1),0)
        for v in range(size+1):
            c = size-v
            for i in (0,1):
                # Identidad de Euler para la exponencial bivariada.
                terms = ((a+b)*(A[i][0]*trees[0][a,b]+A[i][1]*trees[1][a,b])
                         *forest[i][v-a,c-b]
                         for a in range(v+1) for b in range(c+1) if a+b)
                forest[i][v,c] = sum(terms,Q(0))/size
    assert trees[0][1,0] == trees[1][0,1] == 1
    assert trees[1][0,3] == Q(3,2)*513**2
    return trees


def certificate(exp_upper, endpoints, pairs, caps):
    p, q, alpha, eta, pin = map(Q, ("1e-8",".001",".001",".0001",".00117"))
    U,s = map(Q, ("0.00000285918840032462","0.00113259403193381"))
    C = tuple(map(Q, ("1.68954770737199e-8","1.95802514295680e-8")))
    d = tuple(map(Q, ("0.00000687198478039258","0.0000811526827105124")))
    u = tuple(map(Q, ("1.93893583047814e-8","3.14419535336581e-8")))
    un,uf,vn,vf,cn,cf = map(Q, ("2.90715985279773e-7","4.04948406454983e-7",
        ".00237288356039209",".0000470430220249001","1.20093188686862e-8","1.67281979138648e-8"))
    non_simple = tuple(map(Q, ("2.38661346248063e-9","1.17116131155342e-8")))
    root = Q(".9742")
    base = root**5
    t, y = base**3, base**33
    assert y == t**FAR_COST and 36*y < 1
    j,f1,f2,k1,k2,m1,m2,mg = endpoints(root**11)
    assert (j,f1,f2,k1,k2,m1,m2,mg) == (182,2578,3058,3437,2294,420,501,1149)

    # Redondeo racional por exceso, incluso al congelar V en ambos términos de la cola.
    def exp_bound(x):
        scale = 10**24
        return Q(ceil(scale*exp_upper(x)),scale)

    extra = un+uf
    measure = q*exp_bound(Q(".015"))/(1-2*q)
    data = p*exp_bound(Q(".015")+alpha)/(1-2*p)
    assert U < eta
    S = 89*pin+469*U
    eps = S*exp_bound(2*S)
    passive = measure*(1+pin)*(1+eps)/pin
    assert eps < 1 and passive < 1
    ratios = [measure*exp_bound(90*pin)/s]
    simple, V = [], []
    for i,(M,J,f,m,n,sigma,halo) in enumerate((
            (1344,30,f1,m1,3,Q(17,3),414), (1596,36,f2,m2,4,Q(103,15),536))):
        ratios.append(data*exp_bound(eta+sigma*1297*U+pin*halo)/C[i])
        simple.append(C[i]*(1+s)**n*exp_bound(M*extra))
        V.append((1+s)**n*exp_bound(M*extra+J*vn+t*(f*vf+m*s)))

    w = (C[0]*V[0]/t**14, C[1]*V[1]/base**59)
    trees = coefficients(5)
    tail, bulk_ratios = [], []
    for i in (0,1):
        full = w[i]*exp_bound(A[i][0]*d[0]+A[i][1]*d[1])
        bulk_ratios.append(full/d[i])
        small = sum(value*w[0]**a*w[1]**b for (a,b),value in trees[i].items() if a+b <= 4)
        tail.append(full-small)
        assert tail[i] >= sum(value*w[0]**a*w[1]**b for (a,b),value in trees[i].items() if a+b == 5)
        pair = C[i]*V[i]*sum(C[k]*V[k]*sum(count/t**cost for count,cost in pairs[i][k])
                              for k in (0,1))
        short = sum(trees[i][a,b]*(C[0]*V[0])**a*(C[1]*V[1])**b/t**cost
                    for (a,b),cost in caps.items())
        exceptional = C[i]*t**(FAR_COST-2)*V[i]
        ratios.append((exceptional+pair+short+tail[i])/non_simple[i])
        ratios.append((simple[i]+non_simple[i])/u[i])

    L = j*vf+y*mg*s
    G = (1+s)**12*exp_bound(1792*u[0]+1197*u[1]+3648*extra+L)
    Hn = 40*non_simple[0]+27*non_simple[1]
    Hf = j*(cn+cf)+t*(k1*non_simple[0]+k2*non_simple[1])
    ratios.extend((t*Hn*G/un, y*Hf*G/uf, t*(G-1-L-12*s)/vn,
                   y*(G-1-L-12*s)/vf, t*Hn*(G-1)/cn, y*Hf*(G-1)/cf,
                   max(u[0]+3*simple[0]+3*extra,u[1]+2*simple[1]+4*extra)/U))
    assert 0 < min(ratios) and max(ratios) < 1
    assert 0 < min(bulk_ratios) and max(bulk_ratios) < 1
    row = max(sum(a*x for a,x in zip(ar,u)) for ar in A)
    assert row < Q(".00003134") < alpha
    assert row/64 < Q("4.90e-7")

    # La clase pura de cinco CC por sí sola, dentro del mayorante actual.
    assert trees[1][0,5] == Q(125,24)*513**4
    lower = Q("1e-7")*(1+536*pin)/base**59
    obstruction = 513*trees[1][0,5]*lower**5
    assert obstruction > Q(".00117") > alpha
    print(f"p=1e-8, q=.001. Fila {float(row):.12g} < .00003134. "
          f"Tres o más en p=2.5e-9: {float(row/64):.12g}.")
    print(f"Fondo {float(passive):.10f}. Mayor cociente de las otras filas {float(max(ratios)):.12f}. "
          f"Ramas auxiliares {float(max(bulk_ratios)):.12f}.")
    print(f"Colas desde cinco, VV/CC: {float(tail[0]):.10g}, {float(tail[1]):.10g}. "
          f"Cota inferior para cinco CC en el mayorante fijo a p=1e-7: "
          f"{float(obstruction):.12g}, superior a .00117.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para comprobar las desigualdades.")
    folder = Path(__file__).resolve().parent
    def previous(name):
        return run_path(str(folder/f"verificar_H5A_3_4_{name}.py"))
    old = previous("fondo_pasivo")
    adjacent = previous("pares_adyacentes")
    pairs, caps = cost_tables(light_profiles(), old["step_costs"], adjacent["stripe_profiles"])
    certificate(previous("medidas_y_un_dato")["exp_upper"], previous("esquinas")["endpoint_counts"], pairs, caps)
