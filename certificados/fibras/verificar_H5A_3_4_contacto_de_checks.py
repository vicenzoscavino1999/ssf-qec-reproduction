"""HGP (3,4), p=4.5e-8, q=.001. Python 3.10+ y biblioteca estándar."""

from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import permutations, product
from math import ceil, factorial
from pathlib import Path
from runpy import run_path


WEAK = ((486,324),(576,384))
STRONG = ((154,117),(208,129))
# Créditos multiplicados por seis. Se indexan por tipo y por rango desde uno.
ORDINARY = ((20,33,72),(20,20,66,68))
PREVIOUS = ((20,30,64),(20,20,42,64))
GOOD_COST = (73,90)
BULK_COST = (84,118)


def geometry():
    rows = tuple(7 << (3*i) for i in range(4))
    cols = tuple(sum(1 << (3*i+j) for i in range(4)) for j in range(3))
    stripes = rows+cols
    footprints = []
    for bits in range(1,128):
        mask = 0
        for i,stripe in enumerate(stripes):
            if bits >> i & 1:
                mask ^= stripe
        if mask and 2*mask.bit_count() >= 3*bits.bit_count():
            footprints.append(mask)
    assert len(footprints) == 101
    cuts = [(p,p.bit_count()) for p in set(footprints)]
    stable, winners = [], set()
    for state in range(4096):
        hmax, vmax, tied = -1, 1, []
        for phi,v in cuts:
            lit = phi & state
            h = lit.bit_count()
            cmp = h*vmax-hmax*v
            if cmp > 0:
                hmax, vmax, tied = h,v,[]
            if cmp >= 0:
                tied.append((lit,phi ^ lit))
        eligible = 100*(2*hmax-vmax) >= vmax
        stable.append(not eligible)
        if eligible:
            winners.update(tied)
    assert len(winners) == 771 and sum(stable) == 91
    perm = [tuple(1 << (3*r[i//3]+c[i%3]) for i in range(12))
            for r in permutations(range(4)) for c in permutations(range(3))]
    def change(mask, mapping):
        return sum(bit for i,bit in enumerate(mapping) if mask >> i & 1)
    seen, profiles, representatives = set(), [], []
    for lit,dark in sorted(winners):
        if (lit,dark) in seen:
            continue
        orbit = {(change(lit,m),change(dark,m)) for m in perm}
        assert orbit <= winners
        seen.update(orbit)
        representatives.append((lit,lit|dark))
        h, v = lit.bit_count(), (lit|dark).bit_count()
        options = []
        for k,axis in enumerate((rows,cols)):
            choices = {0}
            for stripe in axis:
                available = stripe & lit
                mask = available
                while mask:
                    if v*mask.bit_count() <= (3+k)*h:
                        choices.add(mask)
                    mask = (mask-1) & available
            assigned = []
            for mask in sorted(choices):
                blocked = 0
                for stripe in stripes:
                    if stripe & mask:
                        blocked |= stripe
                assigned.append((mask,mask.bit_count(),blocked))
            options.append(assigned)
        profiles.append((lit,h,v,options))
    assert seen == winners and len(profiles) == 25
    return profiles,stable,representatives


def exposure(representatives, budgets):
    """Centros externos a distancia de lectura de una huella realmente conmutada."""
    def axis(root_degree):
        graph, level = {0: set()}, [0]
        for k in range(4):
            following = []
            degree = root_degree if k%2 == 0 else 7-root_degree
            for parent in level:
                for _ in range(degree-len(graph[parent])):
                    node = len(graph)
                    graph[node] = {parent}
                    graph[parent].add(node)
                    following.append(node)
            level = following
        def ball(start, distance):
            reached = {start}
            for _ in range(distance):
                reached = {v for x in reached for v in graph[x]}
            return reached
        splits = []
        for distance in (2,3):
            balls = [ball(neighbor,distance) for neighbor in graph[0]]
            common = set.intersection(*balls)
            private = [part-common for part in balls]
            assert len(set.union(*private)) == sum(map(len,private))
            sizes = {len(part) for part in private}
            assert len(sizes) == 1
            splits.append((len(common),sizes.pop()))
        return splits, [len(ball(0,k)) for k in range(5)]
    row_axis,cr = axis(4)
    col_axis,vr = axis(3)
    assert row_axis == [(4,6),(9,12)] and col_axis == [(3,6),(10,18)]
    assert (cr,vr) == ([1,4,9,28,57],[1,3,10,21,64])
    assert WEAK == (((vr[2]-1)*(vr[4]-vr[2]),(vr[3]-vr[1])**2),
                    ((cr[3]-cr[1])**2,(cr[2]-1)*(cr[4]-cr[2])))
    closed = ((vr[2]*vr[4],vr[3]**2),(cr[3]**2,cr[2]*cr[4]))
    assert all(WEAK[i][j]+STRONG[i][j] == closed[i][j] for i,j in product(range(2),repeat=2))
    rows = tuple(7 << (3*i) for i in range(4))
    cols = tuple(sum(1 << (3*i+j) for i in range(4)) for j in range(3))
    # Común y privado por eje. Los puntos de partida cubren la huella original.
    result = []
    for ar,br,ac,bc,initial in ((4,6,3,6,(210,252)),(4,6,10,18,(640,784)),
                                (9,12,3,6,(441,513)),(9,12,10,18,(1344,1596))):
        def cover(mask):
            if not mask:
                return 0
            return (ar*ac+br*ac*sum(bool(x&mask) for x in rows)
                    +ar*bc*sum(bool(x&mask) for x in cols)+br*bc*mask.bit_count())
        signatures = {(0,phi.bit_count(),lit.bit_count(),cover(phi)-cover(lit))
                      for lit,phi in representatives}
        ranks = budgets(signatures,lambda _,increment: increment)
        result.append(tuple(x+sum(r) for x,r in zip(initial,ranks)))
    assert result[0] == (414,536)
    return result


def universal_discounts(profiles, stable):
    for lit,h,v,options in profiles:
        credit = 6*(11*(2*h-v)-1)
        top = (3*h//v,4*h//v)
        old = max(max(row[:k]) for row,k in zip(PREVIOUS,top))
        assert h*old <= credit
        choices = [(k,mask,a,blocked) for k in range(2)
                   for mask,a,blocked in options[k] if mask]
        def assign(start, used, blocked, charge, pristine):
            if not pristine or stable[lit & blocked & ~used]:
                assert charge+old*(lit & ~blocked).bit_count() <= credit
            for j in range(start,len(choices)):
                k,mask,a,hit = choices[j]
                if mask & blocked:
                    continue
                for r in range(a,top[k]+1):
                    value = sum(ORDINARY[k][r-a:r])
                    assign(j+1,used|mask,blocked|hit,charge+value,pristine or r == 3+k)
        assign(0,0,0,0,False)
    for i,(good,bad) in enumerate(zip(ORDINARY,PREVIOUS)):
        n = 3+i
        assert max(66*(2*r-n)-sum(good[:r]) for r in range(n+1)) == GOOD_COST[i]
        assert max(66*r-sum(bad[:r]) for r in range(n+1)) == BULK_COST[i]


def group_bound(types, blocks, profiles, stable):
    widths = tuple(3+i for i in types)
    states = list(product(*(range(n+1) for n in widths)))
    index = {r:i for i,r in enumerate(states)}
    rank = [max(12*x//n for x,n in zip(r,widths)) for r in states]
    allowed = {k:{} for k in set(rank)}
    for lit,h,v,options in profiles:
        regular, perfect = set(),set()
        def assign(i, used, claimed, counts):
            if i == len(types):
                if used:
                    move = index[counts]
                    regular.add(move)
                    if h == v and stable[lit & ~used]:
                        perfect.add(move)
                return
            forbidden = used
            for j,mask in enumerate(claimed):
                if blocks[j] != blocks[i]:
                    forbidden |= mask
            for mask,num,blocked in options[types[i]]:
                if not mask & forbidden:
                    assign(i+1,used|mask,claimed+(blocked,),counts+(num,))
        assign(0,0,(),())
        for k,tbl in allowed.items():
            if k*v > 12*h:
                continue
            for move in perfect if k == 12 else regular:
                cost = 11*(2*h-v)-1
                tbl[move] = min(tbl.get(move,cost),cost)
    components = [[i for i,g in enumerate(blocks) if g == b] for b in sorted(set(blocks))]
    discount = [0]
    best = -10**9
    for i,r in enumerate(states):
        if i:
            tbl = allowed[rank[i]]
            moves = (index[a] for a in product(*(range(x+1) for x in r)))
            discount.append(min(tbl[m]+discount[i-m] for m in moves if m in tbl))
        balance = 0
        for ids in components:
            n = sum(widths[j] for j in ids)
            s = sum(r[j] for j in ids)
            balance += 2*s-n if len(ids) == 1 else s-(n-s)%2
        best = max(best,11*balance-discount[i])
    return best


def rooted_trees(n):
    if n == 1:
        yield ()
        return
    for seq in product(range(n),repeat=n-2):
        degree = [1+seq.count(i) for i in range(n)]
        edges = []
        for j in seq:
            i = degree.index(1)
            edges.append((i,j))
            degree[i] -= 1
            degree[j] -= 1
        edges.append(tuple(i for i,d in enumerate(degree) if d))
        todo = [(0,-1)]
        oriented = []
        while todo:
            i,parent = todo.pop()
            for a,b in edges:
                if a == i or b == i:
                    j = b if a == i else a
                    if j != parent:
                        todo.append((j,i))
                        oriented.append((i,j))
        yield tuple(oriented)


def positions(profiles, stable, old_caps, old_trees):
    @lru_cache(None)
    def cap(groups):
        types, blocks = [], []
        for b,(nv,nc) in enumerate(groups):
            types.extend((0,)*nv+(1,)*nc)
            blocks.extend((b,)*(nv+nc))
        return group_bound(types,blocks,profiles,stable)
    pair_weak = [cap(tuple(sorted(((1-i,i),(1-j,j))))) for i,j in ((0,0),(0,1),(1,1))]
    assert pair_weak == [2,13,3]
    counts = [Counter(),Counter()]
    for n in range(1,5):
        for edges in rooted_trees(n):
            for types in product(range(2),repeat=n):
                nv = types.count(0)
                nc = n-nv
                for flags in product((0,1),repeat=n-1):
                    labels = list(range(n))
                    multiplicity = 1
                    for (i,j),strong in zip(edges,flags):
                        multiplicity *= (STRONG if strong else WEAK)[types[i]][types[j]]
                        if strong:
                            source, target = labels[j], labels[i]
                            labels = [target if x == source else x for x in labels]
                    groups = tuple(sorted((sum(g == b and k == 0 for g,k in zip(labels,types)),
                                           sum(g == b and k == 1 for g,k in zip(labels,types)))
                                          for b in set(labels)))
                    gv,gc = groups.count((1,0)),groups.count((0,1))
                    total = gv*73+gc*90+(nv-gv)*84+(nc-gc)*118
                    short = min(total,6*old_caps[nv,nc],6*cap(groups)) if n >= 3 else 0
                    counts[types[0]][nv,nc,total,short] += Q(multiplicity,factorial(n-1))
    for i,table in enumerate(counts):
        for (nv,nc),value in old_trees[i].items():
            assert sum(q for (v,c,_,_),q in table.items() if (v,c) == (nv,nc)) == value
    assert cap(((0,1),(0,1),(0,1))) == 15 and cap(((0,1),)*4) == 16
    return counts


def check_history_merger():
    """Contraste de la comparación local con un rival que necesita tres historias."""
    def mask(*sites):
        return sum(1 << j for j in sites)

    pieces = []
    for i in range(3):
        a,b,l,z,d,e,h,f1,f2 = range(9*i,9*i+9)
        first,second = mask(a,b,l,z),mask(z,l,d,e,h,f1,f2)
        raw,noise = first^second,mask(l,f1,f2)
        pieces.append((first,second,raw,noise))
    passive = (0,1<<27,1<<28,mask(27,28))
    rival = mask(3,12,21,27,28)
    checked = 0
    for order in permutations(range(3)):
        cuts = [rival]+[phi for i in order for phi in pieces[i][:3]]
        def winner(state, choices):
            eligible = [j for j in choices if 100*(2*(state&cuts[j]).bit_count()
                        -cuts[j].bit_count()) >= cuts[j].bit_count()]
            return max(eligible, key=lambda j: (Q((state&cuts[j]).bit_count(),cuts[j].bit_count()),-j),
                       default=None)
        def decode(state):
            trace = []
            while (j := winner(state,range(len(cuts)))) is not None:
                trace.append(j)
                state ^= cuts[j]
            return tuple(trace)
        histories = []
        for first,second,raw,noise in pieces:
            initial = raw^noise
            actions = decode(initial)
            assert actions == (cuts.index(first),cuts.index(second))
            histories.append((first|second,initial,actions))
        def merge(chosen, background):
            state = background
            pending = []
            for _,initial,actions in chosen:
                state ^= initial
                pending.append(list(actions))
            while any(pending):
                j = winner(state,[seq[0] for seq in pending if seq])
                assert j is not None
                yield state,j
                next(seq for seq in pending if seq and seq[0] == j).pop(0)
                state ^= cuts[j]
            yield state,None
        for bits in range(8):
            chosen = [c for i,c in enumerate(histories) if bits >> i & 1]
            for bg in passive:
                full = list(merge(chosen,bg))
                proposed = tuple(j for _,j in full if j is not None)
                agrees = decode(full[0][0]) == proposed
                local = True
                for k,phi in enumerate(cuts):
                    readers = [c for c in chosen if c[0]&phi]
                    for state,j in merge(readers,bg&phi):
                        if winner(state,[k] if j is None else [k,j]) != j:
                            local = False
                assert agrees == local
                if not bg:
                    assert agrees == (len(chosen) < 3)
                checked += 1
    assert checked == 192
    print(f"Recomposición local en {checked} comparaciones abstractas, con un rival "
          "que cambia la historia al reunir tres núcleos.")


def certificate(exp_upper,endpoints,counts,halos):
    p,q,b,alpha,eta = map(Q,("4.5e-8",".001",".001155",".001",".0001"))
    r = Q(".98735")
    a = r**5
    t,Y = a**6,a**66
    assert t**11 == Y and 36*Y < 1
    j,fv,fc,kv,kc,mv,mc,mg = endpoints(r**22)
    assert (j,fv,fc,kv,kc,mv,mc,mg) == (206,2902,3443,3870,2582,474,565,1293)
    def E(x):
        scale = 10**24
        return Q(ceil(scale*exp_upper(x)),scale)
    s = Q(".0011507323105217158")
    C = tuple(map(Q,("8.046034046992379e-8","9.494937762602184e-8")))
    u = tuple(map(Q,("9.009995802141036e-8","3.193504736969793e-7")))
    h = tuple(map(Q,("2.895857605971868e-8","2.545436094070671e-7")))
    d = tuple(map(Q,("1.5580721322445828e-5","7.198586642051703e-5")))
    e = tuple(map(Q,("3.0044774988916885e-5",".0003439303978754843")))
    un,uf,vn,vf,cn,cf = map(Q,("5.735617881321943e-6","1.0552958258595964e-5",
        ".001691642631227813","3.7123402405713755e-5","2.5339406907999296e-7","4.6621952321562814e-7"))
    g = un+uf
    row = max(640*u[0]+441*u[1],784*u[0]+513*u[1])
    weights = b,2*u[0],2*u[1],2*g
    assert halos == [(414,536),(1264,1652),(Q(1719,2),1095),(Q(5247,2),3374)]
    initial = [(210,252),(640,784),(441,513),(1344,1596)]
    per_measure = 280*u[0]+189*u[1]+588*g
    assert max(u) < eta and row < Q(".00023447") < alpha
    assert row*(Q(2,9)**3) < Q(".000002574")
    L0 = 89*b+per_measure
    yq = q*E(Q(".015"))/(1-2*q)
    eps = E(L0)*(E(L0)-1)
    passive = yq*(1+b)*(1+eps)/b
    assert eps < 1 and passive < Q(".9991")
    z = p*E(Q(".015")+alpha)/(1-2*p)
    ratios = [yq*E(90*b+2*per_measure)/s]
    V, simple = [], []
    for i,(n,f,m,J) in enumerate(((3,fv,mv,30),(4,fc,mc,36))):
        def price(table):
            return z*E(eta+2*row+sum(w*line[i] for w,line in zip(weights,table)))
        ratios.append(price(halos)/C[i])
        simple.append(price(initial)*(1+yq)**n)
        V.append((1+s)**n*E(J*vn+t*(f*vf+m*s)))
    x = tuple(v*c for v,c in zip(V,C))
    pair = ((486/t**2+90/t**3+64/t**14+9/t**24,324/t**13+117/t**14+9/t**24),
            (576/t**13+208/t**14+16/t**24,384/t**3+72/t**14+65/t**25))

    def small(i, weights):
        return sum(qt*weights[0]**v*weights[1]**c/a**short
                   for (v,c,_,short),qt in counts[i].items() if v+c >= 3)

    bulk_ratios, tails = [], []
    for i in range(2):
        X = sum(k*y for k,y in zip(WEAK[i],d))
        Z = sum(k*y for k,y in zip(STRONG[i],e))
        good,bad = x[i]/a**GOOD_COST[i], x[i]/a**BULK_COST[i]
        full = E(X)*(good+bad*(E(Z)-1))
        bulk_ratios += [full/d[i],bad*E(X+Z)/e[i]]
        through_four = sum(qt*x[0]**v*x[1]**c/a**total for (v,c,total,_),qt in counts[i].items())
        tails.append(full-through_four)
        assert tails[i] > Q(125,24)*WEAK[i][i]**4*good**5
        two = x[i]*sum(k*y for k,y in zip(pair[i],x))
        ratios.append((x[i]*t**9+two+small(i,x)+tails[i])/h[i])
        ratios.append((simple[i]+h[i])/u[i])
    L = j*vf+Y*mg*s
    G = (1+s)**12*E(89*g+L)
    Hn,Hf = 40*h[0]+27*h[1],j*(cn+cf)+t*(kv*h[0]+kc*h[1])
    ratios.extend((t*Hn*G/un,Y*Hf*G/uf,t*(G-1-L-12*s)/vn,Y*(G-1-L-12*s)/vf,
                   t*Hn*(G-1)/cn,Y*Hf*(G-1)/cf))
    assert 0 < min(ratios) and max(ratios) < 1
    assert 0 < min(bulk_ratios) and max(bulk_ratios) < 1

    # Condición necesaria de estas mismas filas, manteniendo fijos t, b y alfa.
    lower_p = Q("1e-7")
    def exp_lower(w):
        return 1+w+w*w/2+w**3/6
    slo = q*exp_lower(90*b)
    xlo = tuple(lower_p*exp_lower(b*k)*(1+slo)**n*exp_lower(t*m*slo)
                for k,n,m in ((414,3,mv),(536,4,mc)))
    hlo = [xlo[i]*(t**9+sum(k*xj for k,xj in zip(pair[i],xlo)))+small(i,xlo) for i in range(2)]
    glo = (t*(40*hlo[0]+27*hlo[1])+Y*t*(kv*hlo[0]+kc*hlo[1]))*(1+slo)**12*exp_lower(Y*mg*slo)
    level = 89*b+280*hlo[0]+189*hlo[1]+588*glo
    eps_lo = level+3*level**2/2+7*level**3/6
    ratio_lower = q*(1+Q(".015"))*(1+b)*(1+eps_lo)/((1-2*q)*b)
    assert ratio_lower > Q("1.009")
    print(f"p=4.5e-8, q=.001. Fila general {float(row):.12g}. En p=1e-8, tres o más "
          "componentes < 2.574e-6.")
    print(f"Respuesta pasiva {float(passive):.12f}. Máximos de las filas de núcleo y "
          f"auxiliares {float(max(ratios)):.12f}, {float(max(bulk_ratios)):.12f}.")
    print(f"Colas desde cinco por raíz VV y CC: {float(tails[0]):.10g}, {float(tails[1]):.10g}. "
          "En p=1e-7 el cociente pasivo necesario con estos precios supera 1.009.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para comprobar las desigualdades.")
    folder = Path(__file__).resolve().parent
    def previous(name):
        return run_path(str(folder/f"verificar_H5A_3_4_{name}.py"))
    prev = previous("franjas")
    _,old_caps = prev["cost_tables"](prev["light_profiles"](), previous("fondo_pasivo")["step_costs"],
                                    previous("pares_adyacentes")["stripe_profiles"])
    profiles,stable,representatives = geometry()
    universal_discounts(profiles,stable)
    hist = positions(profiles,stable,old_caps,prev["coefficients"](4))
    check_history_merger()
    base = previous("medidas_y_un_dato")
    certificate(base["exp_upper"],previous("esquinas")["endpoint_counts"],hist,
                exposure(representatives,base["budgets"]))
