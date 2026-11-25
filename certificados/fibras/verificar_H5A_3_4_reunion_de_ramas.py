"""HGP (3,4). Reunión de descendientes y saldo de componentes de franjas solapadas."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from math import comb, factorial, prod
from pathlib import Path
from runpy import run_path


HERE = Path(__file__).resolve().parent
prior = run_path(str(HERE / "verificar_H5A_3_4_confluencia_en_generador.py"))
OVERLAP = ((9, 9), (16, 8))
DISJOINT = ((630, 432), (768, 504))
WEAK = ((486,324),(576,384))
STRONG_DISJOINT = ((144,108),(192,120))
TYPES = (16, 9)


def causal_rows(z):
    """Filas de un paso con uno o con dos predecesores de salidas prefijadas."""
    old = prior["previous"]("verificar_H5A_3_4_cadenas_un_predecesor.py")
    winners, rows, cols = old["catalogo"]()
    def kind(mask):
        if mask.bit_count() <= 1:
            return mask.bit_count()
        return 2 if any(mask & col == mask for col in cols) else 3
    def histogram(items, cost):
        hist = [Counter() for _ in range(4)]
        for h,foot in items:
            mult = 3 if foot in rows else 4 if foot in cols else 1
            hist[kind(foot^h)][cost(h)] += Q(1,mult)
        return [sum(n*z**e for e,n in price.items()) for price in hist]
    example = (1, 1|8, 1|16)
    unary = []
    for i,mask in enumerate(example):
        p = histogram(((h,f) for h,f in winners if h&mask), lambda h:(h&~mask).bit_count())
        if i == 0:
            first = p
            unary.append(tuple(12*x for x in p))
        else:
            unary.append(tuple(max(24*x,12*y) for x,y in zip(first,p)))
    # Fijar una celda distinta de cada brazo. En un rectángulo que lee ambas
    # quedan a lo sumo a+b-2 celdas suministradas por los predecesores.
    joint = {}
    for a,b in product((1,2),repeat=2):
        estimates = []
        for other in (1,3,4):
            pair = 1|(1<<other)
            estimates.append(histogram(((h,f) for h,f in winners if h&pair == pair),
                                       lambda h:max(0,h.bit_count()-a-b)))
        joint[a,b] = tuple(12*a*b*max(col) for col in zip(*estimates))
    return unary,joint


def fork_rows(z, tilt, arms):
    """Matrices positivas para composiciones en serie y en paralelo.

    arms[t] domina las redes no vacías que arrancan en una celda y
    terminan con salida del tipo t. Las ramas interiores pueden anidarse.
    """
    unary, joint = causal_rows(z)
    sizes = (1,2,2)
    matrix = []
    for i,row in enumerate(unary):
        costs = []
        for t in range(4):
            branches = 0
            if i:
                branches = 2*sum(x*joint[1,a][t] for x,a in zip(arms,sizes))
                branches += sum(x*y*joint[a,b][t] for x,a in zip(arms,sizes)
                                for y,b in zip(arms,sizes))
            costs.append(tilt*(row[t]+branches))
        matrix.append(costs)
    return matrix


def check_causal():
    arms = tuple(map(Q,(".01585",".00000170",".000776")))
    z = Q(".0012")
    weight = Q(8,5)
    matrix = fork_rows(z,weight,arms)
    bound = [tuple(map(Q,row)) for row in (
        (".00002527225375231",".0158393595024",".000001684080885158",".0007736214116403"),
        (".03125426400785",".4241184896092",".001774760104012",".5214826729192"),
        (".03125570192048",".4125518615121",".001720021611434",".5219492178773"))]
    for i,t in product(range(3),range(4)):
        assert matrix[i][t]+sum(matrix[i][j+1]*bound[j][t] for j in range(3)) < bound[i][t]
    assert all(b<a for b,a in zip(bound[0][1:],arms))
    assert sum(bound[0])<Q(".01664") and max(map(sum,bound))<Q(".979")
    rc = Q("1e-6")
    assert Q(".001")+5*2*rc/3+7*4*rc/(15-8*rc) < z
    print("Redes causales en serie y paralelo: fila ponderada < .979, factor por paso 5/8.",flush=True)
    return bound


def second_origin(bound):
    """Sumar el origen geométrico del segundo brazo cuando éste tiene un solo predecesor por paso."""
    old = prior["previous"]("verificar_H5A_3_4_cadenas_un_predecesor.py")
    winners,rows,cols = old["catalogo"]()
    z,tilt = Q(".0012"),Q(8,5)
    unary,_ = causal_rows(z)
    types = (tuple(1<<j for j in range(12)),
             tuple(x for x in range(1<<12) if x.bit_count()==2 and any(x&c==x for c in cols)),
             tuple(x for x in range(1<<12) if x.bit_count()==2 and
                   not any(x&s==x for s in rows+cols)))
    counts = tuple(map(len,types))
    through = (1,3,6)
    assert counts == (12,18,36)
    def multiplicity(f):
        return 3 if f in rows else 4 if f in cols else 1
    reverse = [[Q(0)]*3 for _ in types]
    for t,outputs in enumerate(types):
        for a,parents in enumerate(types):
            reverse[t][a] = Q(counts[a],counts[t])*unary[a][t+1]
            for out in outputs:
                upper = 0
                for h,f in winners:
                    if f^h != out:
                        continue
                    size = h.bit_count()
                    value = through[a]*size*z**(size-1)
                    if a:
                        pairs = sum(s&h==s for s in parents)
                        if pairs:
                            value += pairs*(1-2*z)*z**(size-2)
                    upper += 12*value/multiplicity(f)
                assert upper <= reverse[t][a]
    # g[t] cuenta los caminos inversos desde un estado terminal, incluido el vacío.
    g = tuple(map(Q,("1.191","1.001","1.156")))
    assert all(1+tilt*sum(v*x for v,x in zip(row,g)) <= cap for row,cap in zip(reverse,g))
    one_rect = []
    for a in (1,1|8,1|16):
        price = 0
        for h,f in winners:
            if not a&h:
                continue
            for j in range(12):
                y = 1<<j
                if not h&y or not (a&h&~y):
                    continue
                fresh = h & ~(a|y)
                e = fresh.bit_count()
                for k,choices in enumerate(types):
                    other = sum(bool(s&y and s&fresh) for s in choices) if k else 0
                    row = 12*through[k]*z**e
                    if other:
                        row += 12*other*(1-z)*z**(e-1)
                    price += g[k]*row/multiplicity(f)
        one_rect.append(price)
    arrival = [12*one_rect[0]]+[max(24*one_rect[0],12*x) for x in one_rect[1:]]
    costs = [tilt*sum(((i==j)+bound[i][j+1])*x for j,x in enumerate(arrival)) for i in range(3)]
    assert all(x < limit for x,limit in zip(costs,(11724,74563,77747)))
    print("Dos orígenes, sumando las posiciones del brazo lineal: "
          f"constantes por tipo {tuple(float(x) for x in costs)}.",flush=True)
    return costs


def shapes(v, c, minimum=(0, 1)):
    if v+c == 0:
        yield ()
        return
    for i in range(v+1):
        for j in range(c+1):
            if (i,j) >= minimum:
                for rest in shapes(v-i, c-j, (i,j)):
                    yield ((i,j),)+rest


def seed_tree(v, c, degrees):
    """Peso de árboles etiquetados cuando la arista i,j cuesta grados[i,j]/tipos[j]."""
    a, b = Q(degrees[0][0], 16), Q(degrees[0][1], 9)
    d = Q(degrees[1][1], 9)
    assert b == Q(degrees[1][0], 16)
    if not v or not c:
        n, weight = (v, a) if v else (c, d)
        return 1 if n == 1 else n**(n-2)*weight**(n-1)
    return b*(v*a+c*b)**(v-1)*(v*b+c*d)**(c-1)


def connecting_tree(blocks, degrees=DISJOINT):
    """Cofactor del laplaciano de los bloques con sus posibles extremos."""
    n = len(blocks)-1
    mat = [[Q(0) for _ in range(n)] for _ in range(n)]
    for i, (v,c) in enumerate(blocks):
        for j in range(i):
            w, d = blocks[j]
            edge = v*w*Q(degrees[0][0],16)+(v*d+c*w)*Q(degrees[0][1],9)+c*d*Q(degrees[1][1],9)
            mat[j][j] += edge
            if i < n:
                mat[i][i] += edge
                mat[i][j] -= edge
                mat[j][i] -= edge
    det = 1
    for i in range(n):
        pivot = mat[i][i]
        assert pivot > 0
        det *= pivot
        for j in range(i+1,n):
            ratio = mat[j][i]/pivot
            for k in range(i+1,n):
                mat[j][k] -= ratio*mat[i][k]
    return det


def block_count(shape):
    weight = connecting_tree(shape)
    for v,c in shape:
        weight *= Q(seed_tree(v,c,OVERLAP),factorial(v)*factorial(c))
    return weight/prod(factorial(m) for m in Counter(shape).values())


def composition(blocks):
    return tuple(map(sum,zip(*blocks)))


def nested_shapes(v,c):
    for outer in shapes(v,c):
        choices = [combinations_with_replacement(tuple(shapes(*size)),count)
                   for size,count in Counter(outer).items()]
        for selection in product(*choices):
            yield tuple(sorted(sum(selection,())))


@lru_cache(None)
def internal_count(blocks):
    weight = connecting_tree(blocks,STRONG_DISJOINT)
    for v,c in blocks:
        weight *= Q(seed_tree(v,c,OVERLAP),factorial(v)*factorial(c))
    return weight/prod(factorial(m) for m in Counter(blocks).values())


def nested_count(groups):
    sizes = tuple(map(composition,groups))
    top = connecting_tree(sizes,WEAK)
    return Q(top*prod(map(internal_count,groups)),prod(factorial(m) for m in Counter(groups).values()))


@lru_cache(None)
def vectors(n, width):
    return tuple(combinations_with_replacement(range(width+1), n))


def clean_tuple(values):
    return tuple(sorted(x for x in values if x))


@lru_cache(None)
def splits(ranks):
    """Un reparto por órbita de donantes que tienen el mismo rango."""
    options = []
    for r,m in sorted(Counter(ranks).items()):
        options.append([(clean_tuple(chosen), clean_tuple(r-x for x in chosen))
                        for chosen in vectors(m,r)])
    return tuple((tuple(sorted(sum((x for x,_ in choices),()))),
                  tuple(sorted(sum((x for _,x in choices),())))) for choices in product(*options))


def local_prices():
    """Catálogo conjunto por número de luces de cada donante VV/CC.

    Cada luz asignada tiene un donante. Las dos orientaciones se asignan
    en celdas disjuntas; en una misma línea se permiten varios donantes.
    """
    contact = prior["previous"]("verificar_H5A_3_4_contacto_de_checks.py")
    far = prior["previous"]("verificar_H5A_3_4_propagacion_lejana.py")
    _, stable, reps = contact["geometry"]()
    rows, cols = far["ROWS"], far["COLS"]

    @lru_cache(None)
    def sums_on_lines(counts):
        if not counts:
            return {()}
        left = counts[0]
        parts = {x for n in range(1,left+1) for x in combinations_with_replacement(range(1,left+1),n)
                 if sum(x) == left}
        return {tuple(sorted(p+r)) for p in parts for r in sums_on_lines(counts[1:])}

    def options(axis):
        return [tuple(sorted(x for stripe in axis if (x:=(mask&stripe).bit_count()))) for mask in range(4096)]

    by_row, by_col = options(rows), options(cols)
    prices = {r:{} for r in (3,4,6,8,9,12)}
    for lit, foot in reps:
        h, m = lit.bit_count(), foot.bit_count()
        delta = 2*h-m
        signatures = set()
        perfect = set()
        used = lit
        while used:
            accepted = h==m and stable[lit ^ used]
            rowmask = used
            while True:
                pair = (by_row[rowmask],by_col[used ^ rowmask])
                signatures.add(pair)
                if accepted:
                    perfect.add(pair)
                if not rowmask:
                    break
                rowmask = (rowmask-1)&used
            used = (used-1)&lit
        for level,out in prices.items():
            if level*m > 12*h:
                continue
            for r,cidx in perfect if level==12 else signatures:
                vv,cc = sums_on_lines(r),sums_on_lines(cidx)
                for v in vv:
                    if v and v[-1]*m > 3*h:
                        continue
                    for c in cc:
                        if c and c[-1]*m > 4*h:
                            continue
                        key = v,c
                        if delta < out.get(key, 1000):
                            out[key] = delta
    return prices


def make_debit(prices):
    @lru_cache(None)
    def debit(vv, cc):
        if not vv and not cc:
            return (0,0)
        rank = max(4*max(vv,default=0),3*max(cc,default=0))
        costs = prices[rank]
        limited, free = 1000, 1000
        for v,vr in splits(vv):
            for c,cr in splits(cc):
                delta = costs.get((v,c))
                if delta is None:
                    continue
                left, right = debit(vr,cr)
                free = min(free,delta+right)
                if len(v)+len(c)<=2:
                    limited = min(limited,delta+left)
        assert free <= limited < 1000
        return limited, free
    return debit


def separated_prices():
    """Dos conjuntos de donantes sin una fila ni una columna común entre sus asignaciones."""
    co = prior["previous"]("verificar_H5A_3_4_contacto_de_checks.py")
    far = prior["previous"]("verificar_H5A_3_4_propagacion_lejana.py")
    _,stable,reps = co["geometry"]()
    rows,cols = far["ROWS"],far["COLS"]
    blocked = [sum((1<<k) for k in range(12) if any((1<<k)&s and s&m for s in rows+cols))
               for m in range(4096)]
    @lru_cache(None)
    def parts(n):
        return [x for k in range(1,n+1) for x in combinations_with_replacement(range(1,n+1),k)
                if sum(x)==n] if n else [()]
    @lru_cache(None)
    def dissect(mask,axis):
        return {tuple(sorted(sum(s,()))) for s in product(*(parts((b&mask).bit_count()) for b in axis))}
    @lru_cache(None)
    def allocation(mask):
        out = set()
        a = mask
        while True:
            out.update(product(dissect(a,rows),dissect(mask^a,cols)))
            if not a:
                break
            a = (a-1)&mask
        return out
    result = {r:{} for r in (3,4,6,8,9,12)}
    for lit,foot in reps:
        h,w = lit.bit_count(),foot.bit_count()
        a = lit
        regular,perfect = set(),set()
        while a:
            available = lit & ~blocked[a]
            b = available
            while b:
                if a<b:
                    for left,right in product(allocation(a),allocation(b)):
                        max_v = max(left[0]+right[0],default=0)
                        max_c = max(left[1]+right[1],default=0)
                        if max_v*w<=3*h and max_c*w<=4*h:
                            key = tuple(sorted((left,right)))
                            regular.add(key)
                            if h==w and stable[lit & ~(a|b)]:
                                perfect.add(key)
                b = (b-1)&available
            a = (a-1)&lit
        for rank,costs in result.items():
            if rank*w<=12*h:
                for key in perfect if rank==12 else regular:
                    costs[key] = min(costs.get(key,1000),2*h-w)
    return result


def cut_debit(ordinary, both):
    @lru_cache(None)
    def moves(side):
        vv,cc=side
        return tuple(((v,c),(x,y)) for v,x in splits(vv) for c,y in splits(cc))
    empty = (),()
    @lru_cache(None)
    def cost(left,right):
        if left==right==empty:
            return 0
        rank = max(4*max(left[0]+right[0],default=0),3*max(left[1]+right[1],default=0))
        single,paired = ordinary[rank],both[rank]
        best = 1000
        for a, ar in moves(left):
            for b, br in moves(right):
                if a==empty or b==empty:
                    pay = single.get(b if a==empty else a)
                else:
                    pay = paired.get((a,b) if a<=b else (b,a))
                if pay is not None and pay<best:
                    remaining = cost(ar,br) if ar<=br else cost(br,ar)
                    best = min(best,pay+remaining)
        assert best<1000
        return best
    return cost


@lru_cache(None)
def block_states(v,c):
    n, width = v+c,3*v+4*c
    options = {}
    for vv in vectors(v,3):
        for cc in vectors(c,4):
            r = sum(vv)+sum(cc)
            if n == 1:
                if 2*r < width:
                    continue
                balance = 2*r-width
            else:
                if r > width-2*(n//2):
                    continue
                balance = r-(width-r)%2
            options[clean_tuple(vv),clean_tuple(cc)] = balance
    return options


@lru_cache(None)
def combined_states(shape, block=block_states):
    if not shape:
        return {((),()):0}
    v,c = shape[0]
    out = {}
    for (x,y),first in block(v,c).items():
        for (z,t),rest in combined_states(shape[1:],block).items():
            key = tuple(sorted(x+z)),tuple(sorted(y+t))
            out[key] = max(out.get(key,-100),first+rest)
    return out


def triple_reader(trees, nmax):
    """Extender una terna de lectores distintos hasta cualquier grupo conexo de grado ≤ nmax."""
    def multiply(left, right):
        dest = Counter()
        for (v,c),x in left.items():
            for (w,d),y in right.items():
                if v+c+w+d <= nmax:
                    dest[v+w,c+d] += x*y
        return dest
    # Serie de cada árbol enraizado. Elegir el trío en el rectángulo y
    # decorar sus tres vértices. Una raíz exterior llega por un camino.
    level = [Counter(),Counter()]
    for v in range(4):
        c = 3-v
        forest = Counter({(0,0):Q(1)})
        for kind in (0,)*v+(1,)*c:
            forest = multiply(forest,trees[kind])
        for i,m in enumerate((v,c)):
            if not m:
                continue
            count = (30,36)[i]*comb(40-(i==0),v-(i==0))*comb(27-(i==1),c-(i==1))
            for key,branches in forest.items():
                level[i][key] += count*branches
    output = [row.copy() for row in level]
    for depth in range(nmax-3):
        extended = []
        for i in (0,1):
            children = Counter()
            for j in (0,1):
                degree = OVERLAP[i][j]+DISJOINT[i][j]
                for key,n in level[j].items():
                    children[key] += degree*n
            extended.append(multiply(trees[i],children))
        level = extended
        for dst,add in zip(output,level):
            dst.update(add)
    return output


def block_caps(hist, triple_volume, nmax=6, circuit_hist=None, *, states=combined_states):
    debit = make_debit(local_prices())
    caps, circuit_caps = {}, {}
    total_degree = tuple(tuple(o+d for o,d in zip(ro,rd)) for ro,rd in zip(OVERLAP,DISJOINT))
    # Un subconjunto de paridad física uno puede tener varias componentes conexas.
    best_subset = {(0,0):0}
    for n in range(1,nmax+1):
        for v in range(n+1):
            c = n-v
            proper = max(best_subset.get((v-1,c),0),best_subset.get((v,c-1),0))
            separated = proper
            for x in range(v+1):
                for y in range(c+1):
                    if 0<x+y<n:
                        separated = max(separated,best_subset[x,y]+best_subset[v-x,c-y])
            idx = 0 if v else 1
            table = hist if circuit_hist is None else circuit_hist
            baseline = max(0,max(table[idx][v,c]))
            if v and c:
                baseline = min(baseline,max(table[1][v,c]))
            best_subset[v,c] = max(baseline,separated)
            if n<3:
                continue
            totals = Counter()
            by_level = [Counter(), Counter()]
            for shape in shapes(v,c):
                weight = block_count(shape)
                totals[0] += weight
                lim = free = 0
                for (vv,cc),balance in states(shape).items():
                    pay_l, pay_f = debit(vv,cc)
                    lim,free = max(lim,balance-pay_l),max(free,balance-pay_f)
                by_level[0][lim] += weight
                by_level[1][free] += weight
            assert totals[0] == seed_tree(v,c,total_degree)/Q(factorial(v)*factorial(c))
            for i,m in enumerate((v,c)):
                if not m:
                    continue
                key = i,v,c
                norm = Q(m*16**v*9**c,TYPES[i])
                for k in range(1,1+max(max(hist[i][v,c]),max(table[i][v,c]))):
                    low, high = [norm*sum(w for h,w in row.items() if h>=k) for row in by_level]
                    value = min(high,low+triple_volume[i].get((v,c),0))
                    caps[key,k] = value
                    if k>proper:
                        circuit_caps[key,k] = value
    return caps, circuit_caps


def separated_caps(hist, circuit_hist, nmax, *, states=combined_states):
    costs = local_prices()
    base = make_debit(costs)
    debit = cut_debit(costs,separated_prices())
    @lru_cache(None)
    def ceiling(left,right):
        peak = 0
        if not right:
            return max([0]+[b-base(*r)[1] for r,b in states(left).items()])
        for a,x in states(left).items():
            for b,y in states(right).items():
                if x+y>peak:
                    pay = debit(a,b) if a<=b else debit(b,a)
                    peak = max(peak,x+y-pay)
        return peak
    caps, circuit = {},{}
    best = {(0,0):0}
    full = tuple(tuple(o+d for o,d in zip(ro,rd)) for ro,rd in zip(OVERLAP,DISJOINT))
    for n in range(1,nmax+1):
        for v in range(n+1):
            c=n-v
            proper = max(best.get((v-1,c),0),best.get((v,c-1),0))
            mx = min(max(circuit_hist[i][v,c]) for i,m in enumerate((v,c)) if m)
            split = [best[x,y]+best[v-x,c-y] for x in range(v+1) for y in range(c+1) if 0<x+y<n]
            best[v,c] = max([proper,mx]+split)
            if n<3:
                continue
            by_cap = Counter()
            for nested in nested_shapes(v,c):
                pieces = tuple(sorted(sum(nested,())))
                bound = ceiling(pieces,())
                cuts = set()
                # Cada componente fuerte individual y, además, los cortes equilibrados.
                for j in range(1,len(nested)//2+1):
                    for ids in combinations(range(len(nested)),j):
                        size = sum(sum(composition(nested[k])) for k in ids)
                        if j>1 and abs(2*size-n)>2:
                            continue
                        left = tuple(sorted(sum((nested[k] for k in ids),())))
                        right = tuple(sorted(sum((g for k,g in enumerate(nested) if k not in ids),())))
                        cuts.add(tuple(sorted((left,right))))
                for cut in cuts:
                    bound = min(bound,ceiling(*cut))
                by_cap[bound] += nested_count(nested)
            assert sum(by_cap.values()) == seed_tree(v,c,full)/Q(factorial(v)*factorial(c))
            for i,m in enumerate((v,c)):
                if not m:
                    continue
                norm = Q(m*16**v*9**c,TYPES[i])
                high = max(max(hist[i][v,c]),max(circuit_hist[i][v,c]))
                for k in range(1,high+1):
                    value = norm*sum(w for t,w in by_cap.items() if t>=k)
                    caps[(i,v,c),k] = value
                    if k>proper:
                        circuit[(i,v,c),k] = value
    return caps,circuit


def tails(row):
    return [sum(x for j,x in row.items() if j>=k) for k in range(max(row)+1)]


def circuit_histograms(hc, hp, geometric):
    """Toda paridad y la subclase con paridad uno en cada dato marcado."""
    universal = prior["apply_caps"](hc,geometric)
    common = {}
    for i in (0,1):
        for key, row in hp[i].items():
            if sum(key)<2:
                continue
            t = tails(row)
            for k in range(1,max(universal[i][key])+1):
                common[(i,*key),k] = t[k] if k<len(t) else 0
    # En una pareja, cero o un real deja saldo a lo sumo cero en todo patrón.
    universal = prior["apply_caps"](universal,{key:cap for key,cap in common.items() if sum(key[0][1:])==2})
    all_real = prior["apply_caps"](universal,common)
    return universal,all_real


def surcharge(universal, real, law, discounts, Y, anchor, ball):
    ordinary, marked = [Counter(),Counter()], [Counter(),Counter()]
    for i,tab in enumerate(universal):
        for key,row in tab.items():
            v,c = key
            if v+c<2:
                continue
            top, lower = tails(row),tails(real[i][key])
            many = 1-(1-law[0][1])**v*(1-law[1][1])**c
            factor = discounts.get(key,1)*many
            for k,t in enumerate(top[1:],1):
                if k>=len(lower) or lower[k]<t:
                    ordinary[i][key] += factor*t*(Y**(-k)-Y**(1-k))
                    marked[i][key] += factor*t*(anchor(ball,i,k)-anchor(ball,i,k-1))
    return ordinary,marked


def circuit_law(kind,rate,first,theta):
    """Ley mayorante de la fase 5, con la masa múltiple calculada en la esquina elegida."""
    width = 3+kind
    wait,ctrl = 2*rate/(3-2*rate),8*rate/(15-8*rate)
    one = (12-width)*wait+width*ctrl
    total = (1+wait)**(12-width)*(1+ctrl)**width-1
    many = 1-one/total
    assert 0 < many < 1
    masses = [first*wait+ctrl/2]+[ctrl]*(width-1)+[(12-width-first)*wait+ctrl/2]
    assert sum(masses)==one
    ordinary = (2,3,6) if not kind else (2,2,4,6)
    isolated = (2,3,6) if not kind else (2,2,6,6)
    capacity,general,disjoint,weak = (Counter() for _ in range(4))
    for b in range(width+1):
        prob = comb(width,b)*theta**b*(1-theta)**(width-b)
        for h,mass in enumerate(masses):
            lit = width-h
            k = min(width,lit+b)
            val = prob*mass/total
            capacity[k] += val
            general[6*lit-sum(ordinary[:lit])+4*b] += val
            if 2*k>=width:
                disjoint[6*(2*k-width)-sum(ordinary[:k])] += val
                weak[6*(2*k-width)-sum(isolated[:k])] += val
        general[6*width-sum(ordinary)+4*b] += many*prob
        disjoint[6*width-sum(ordinary)] += many*prob
        weak[6*width-sum(isolated)] += many*prob
    assert sum(capacity.values())==1-many and sum(general.values())==1
    return total,many,capacity,general,disjoint,weak


def singleton_cost(kind,rate,first,theta,Y):
    """Precio para una fibra aislada en un núcleo no simple, sin el descuento de cobertura."""
    n = 3+kind
    wait,ctrl = 2*rate/(3-2*rate),8*rate/(15-8*rate)
    total = (1+wait)**(12-n)*(1+ctrl)**n-1
    multiple = total-(12-n)*wait-n*ctrl
    # Con h omisiones el síndrome observado comienza con n-h luces
    # intrínsecas y los bits de lectura se suman por XOR.
    def binomial(k,j):
        return comb(k,j)*theta**j*(1-theta)**(k-j)
    prices = [min(Y,Y**2+(7,12)[kind]*theta*Y)]
    for h in range(1,n+1):
        admit = sum(binomial(n-h,x)*binomial(h,y) for x in range(n-h+1) for y in range(h+1)
                    if 2*(n-h-x+y)>=n)
        prices.append(Y*admit)
    gate = ctrl*sum((x+y)/2 for x,y in zip(prices,prices[1:]))
    # Hay doce ubicaciones del cable. Antes del primer color y después del
    # último también hay una espera. Faltas primitivas sólo quitan monomios.
    candidates = []
    for gates in combinations(range(1,11),n):
        if any(t>first+j for j,t in enumerate(gates)):
            continue
        waits = sum(prices[sum(g<=t for g in gates)] for t in range(12) if t not in gates)
        candidates.append((wait*waits+gate+multiple*Y)/total)
    return max(candidates)


PHEN = dict(name="Reunión y rangos por componentes",
    rates=("5.05e-7",".001145",".004",".00089","1.01"),
    m=("9.49861845035664e-7","1.10689693964587e-6"),
    u=("6.67063054983046e-7","7.13512899270768e-7"),
    h=("3.67738122532257e-8","1.19301583924547e-7"),
    physical=("4.27336433978205e-9","1.6434786637199e-8"),
    outer=(".00115","4.12173989557303e-7","1.44535127107231e-6"),
    d=(".00006588173634",".00016234177090"),e=(".00006588173634",".00044142241972"))

CIRCUIT = (
    ("VV primero","9.1e-8",4,".561",".73",
     (("1.30066186496466e-6","1.4900078306743e-6"),
      ("4.19474589177925e-7","7.85459887202955e-7"),
      ("4.38753822596231e-8","2.16386790209275e-7"),
      ("5.71314934016118e-9","3.22332640564612e-8"),
      "6.35296551279425e-7","1.47196462480623e-6")),
    ("Cualquier orden de columnas","8.3e-8",7,".558",".75",
     (("1.18607220300996e-6","1.35858588969823e-6"),
      ("6.02388864288153e-7","7.19479830733358e-7"),
      ("6.89987543842338e-8","2.18445739174788e-7"),
      ("7.39749426522782e-9","3.24940051949566e-8"),
      "6.12887399318424e-7","1.24979275715443e-6")))


def verify_data(v15,ball,halos,series,disc,hp):
    a,s = Q(".567"),Q(".685")
    w = PHEN
    reg,marked = prior["physical_prices"](v15,a,Q(w["outer"][0]),ball,disc,hp)
    marked[0][0,0] = v15["single_tail"](tuple(map(Q,w["m"])))
    v11 = v15["v11"]
    ker = v15["corrected_kernel"](v11["prop"]["kernels"](ball,a**6,s),w,reg,a**6/s**2)
    d,e = [tuple(map(Q,w[key])) for key in ("d","e")]
    eff = v11["effective"]
    eff["certify"](a**6,s,halos,reg,marked,ker,lambda x,i:eff["full_tail"](series,a,x,i,d,e),w)


def verify_circuit(fase,far,ball,anchor,halos,trees,crowded,disc,geom,trip,hc,hr):
    theta = fase["THETA"]
    for name,rate,first,a,s,w in CIRCUIT:
        rate,a,s = map(Q,(rate,a,s))
        law = []
        for i,g in ((0,first),(1,7)):
            # Misma ley, incluidos los controles de calendario de la fase 5.
            assert circuit_law(i,Q("7.5e-8"),g,theta) == fase["laws"](i,Q("7.5e-8"),g)
            law.append(circuit_law(i,rate,g,theta))
        reg,marked = fase["polys"](hr,trees,crowded,disc,geom,trip,law,far["N"],anchor,ball,a**6)
        regular_add,marked_add = surcharge(hc,hr,law,disc,a**6,anchor,ball)
        for i,g in ((0,first),(1,7)):
            reg[i].update(regular_add[i])
            marked[i].update(marked_add[i])
            key = (1-i,i)
            new = disc[key]*singleton_cost(i,rate,g,theta,a**6)
            assert 0 < new < reg[i][key]
            reg[i][key] = new
        row,highest,rest,passive,active = fase["certificate"](law,reg,marked,halos,far,ball,rate,a,s,w)
        print(f"{name}, rc={float(rate):g}, rm=.001. Fila {float(row):.12g} sobre .002.")
        print(f"  Máximo {float(highest):.12f}, pasivo {float(passive):.12f}, activo {float(active):.12f}.")
        print(f"  Colas VV, CC {float(rest[0]):.9g}, {float(rest[1]):.9g}. "
              f"Masa múltiple máxima {float(max(x[1] for x in law)):.9g}.",flush=True)


def main():
    second_origin(check_causal())
    v15 = prior["previous"]("verificar_H5A_3_4_contacto_oculto.py")
    fase = prior["previous"]("verificar_FASE5_circuito_sindrome.py")
    ball,halo,series,geometry,hist,_ = v15["v11"]["data_tables"]()
    discounts = v15["v11"]["geometry_discounts"](geometry)
    base,trees,crowded,disc,geom,trip,hc_halo = prior["circuit_data"](fase)
    limits = prior["triple_limits"]()
    hp = prior["refine"](hist,limits)
    hc = prior["refine"](base,limits,circuit=True)
    old_p,old_c = {},{}
    for i in (0,1):
        cp,cc = prior["aligned_limits"](hp,i)
        old_p.update(cp)
        old_c.update(cc)
    hp,hc = prior["apply_caps"](hp,old_p),prior["apply_caps"](hc,old_c)
    volume = triple_reader(trees,8)
    cp,cc = block_caps(hp,volume,8,hc)
    more_p,more_c = separated_caps(hp,hc,8)
    for dst,extra in ((cp,more_p),(cc,more_c)):
        for key,value in extra.items():
            dst[key] = min(value,dst.get(key,value))
    address = (1,2,3),7
    before = hp[1][2,3][7]
    after = min(before,cp[address])
    assert after < before*Q(".047")
    print(f"Dos VV y tres CC, raíz CC, saldo siete: {before} -> {after}.",flush=True)
    hp = prior["apply_caps"](hp,cp)
    hc,hr = circuit_histograms(hc,hp,cc)
    verify_data(v15,ball,halo,series,discounts,hp)
    verify_circuit(fase,v15["v11"]["prop"],ball,v15["v10"]["anchor"],hc_halo,trees,crowded,
                   disc,geom,trip,hc,hr)


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para comprobar las desigualdades.")
    main()
