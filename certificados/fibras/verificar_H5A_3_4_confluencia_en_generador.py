"""Certificado HGP (3,4) con un lector común para las primeras extinciones."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import permutations, product
from math import comb, factorial, prod
from pathlib import Path
from runpy import run_path


DIR = Path(__file__).resolve().parent
PAIRS = ((0, 1), (0, 2), (1, 2))


def previous(name):
    return run_path(str(DIR / name))


def limited_budget(types, edges, profiles, stable, lines):
    """Techo de saldo cuando cada primer apagado tiene a lo sumo dos donantes."""
    widths = tuple(3+i for i in types)
    ranks = list(product(*(range(n+1) for n in widths)))
    level = {r: max(12*x//n for x, n in zip(r, widths)) for r in ranks}
    moves = {k: {} for k in set(level.values())}
    for lit, h, size, options in profiles:
        generic, perfect = set(), set()
        for (i, j), relation in zip(PAIRS, edges):
            for left, _, left_block in options[types[i]]:
                for right, _, _ in options[types[j]]:
                    used = left | right
                    if not used or left & right:
                        continue
                    if left and right:
                        if relation == 0 or (relation == 1 and left_block & right):
                            continue
                        if types[i] == types[j]:
                            same = any(left & a and right & a for a in lines[types[i]])
                            if same != (relation >= 3):
                                continue
                    entry = [0, 0, 0]
                    entry[i], entry[j] = left.bit_count(), right.bit_count()
                    step = tuple(entry)
                    generic.add(step)
                    if h == size and stable[lit & ~used]:
                        perfect.add(step)
        cost = 2*h-size
        for k, choices in moves.items():
            if k*size <= 12*h:
                for step in perfect if k == 12 else generic:
                    choices[step] = min(cost, choices.get(step, cost))
    blocks = [{i} for i in range(3)]
    for (i, j), relation in zip(PAIRS, edges):
        if relation == 4:
            touched = [g for g in blocks if i in g or j in g]
            blocks = [g for g in blocks if g not in touched] + [set.union(*touched)]
    debit = {(0, 0, 0): 0}
    answer = 0
    for r in ranks[1:]:
        debit[r] = min(cost+debit[tuple(x-y for x, y in zip(r, step))]
                       for step, cost in moves[level[r]].items() if all(y<=x for x,y in zip(r,step)))
        balance = 0
        for g in blocks:
            n, s = sum(widths[i] for i in g), sum(r[i] for i in g)
            if len(g) == 1:
                balance += 2*s-n
            elif s > n-2:
                break
            else:
                balance += s - (n-s)%2
        else:
            answer = max(answer, balance-debit[r])
    return answer


def canonical(types, edges):
    rel = dict(zip(PAIRS, edges))
    return min((tuple(types[i] for i in p),
                tuple(rel[tuple(sorted((p[i],p[j])))] for i,j in PAIRS))
               for p in permutations(range(3)))


def triple_limits():
    iso = previous("verificar_H5A_3_4_medidas_aisladas.py")
    contact = previous("verificar_H5A_3_4_contacto_de_checks.py")
    profiles, stable, _ = contact["geometry"]()
    lines = (iso["ROWS"],iso["COLS"])

    @lru_cache(None)
    def caps(types, edges):
        return (limited_budget(types, edges, profiles, stable, lines),
                iso["triple_budget"](types, edges, profiles, stable))

    ordinary, joint = Counter(), Counter()
    joint_shapes = {}
    trees = tuple(contact["rooted_trees"](3))
    for types in product((0,1), repeat=3):
        options = [iso["positions"](types[i],types[j]) for i,j in PAIRS]
        for edges in product(*[(0,*opt.keys()) for opt in options]):
            if sum(bool(x) for x in edges) < 2:
                continue
            # Dos igualdades de coordenada fija fuerzan la tercera si también es vecina.
            if len(set(types))==1 and sum(x>=3 for x in edges)==2 and 0 not in edges:
                continue
            rel = dict(zip(PAIRS,edges))
            costs = []
            for tree in trees:
                count = 1
                for i,j in tree:
                    edge = rel[tuple(sorted((i,j)))]
                    if not edge:
                        break
                    count *= iso["positions"](types[i], types[j])[edge]
                else:
                    costs.append(count)
            low, high = caps(*canonical(types,edges))
            assert low <= high
            key = (types[0],types.count(0),types.count(1))
            for k in range(1,high+1):
                weight = Q(min(costs),2)
                if k <= low:
                    ordinary[key,k] += weight
                elif 0 not in edges:
                    joint[key,k] += weight
                    joint_shapes.setdefault((key,k), set()).add(tuple(sorted(edges)))
    result={}
    for key,k in ordinary.keys() | joint.keys():
        i,v,c=key
        max_joint = (30,36)[i]*comb(40-(i==0),v-(i==0))*comb(27-(i==1),c-(i==1))
        shapes=joint_shapes.get((key,k),set())
        if key==(1,0,3) and shapes and all(max(e)<=2 for e in shapes):
            # Tres columnas físicas distintas e incluidas en la misma vecindad de tamaño tres.
            max_joint=min(max_joint, 4*9**3)
        if key==(0,3,0) and shapes=={(3,3,3)}:
            # Una fila fija y tres vecindades de columna disjuntas leídas por el mismo centro.
            max_joint=min(max_joint, 9*3*3)
        result[key,k]=ordinary[key,k]+min(joint[key,k],max_joint)
    assert joint_shapes[(1,0,3),4]=={(1,1,2),(1,2,2),(2,2,2)}
    for shape in joint_shapes[(1,0,3),4]:
        assert caps((1,1,1),shape) == (1,4)
    assert ordinary[(1,0,3),4]==0 and result[(1,0,3),4]==2916
    assert ordinary[(0,3,0),3]==0 and result[(0,3,0),3]==81
    return result


def refine(hist, limits, circuit=False):
    """Opera sobre techos acumulados, conservando las colas y todos los órdenes altos."""
    out = [dict(table) for table in hist]
    for i,table in enumerate(hist):
        for (v,c),old in table.items():
            if i==1 and (v,c)==(0,2):
                tails=[sum(n for j,n in old.items() if j>=k) for k in range(3)]
                tails[1]=min(tails[1],104)
                tails[2]=min(tails[1],tails[2],12)
                out[i][v,c]={j:x-y for j,(x,y) in enumerate(zip(tails,tails[1:]+[0])) if x>y}
                continue
            if v+c != 3:
                continue
            high=max(old)
            tails=[sum(n for j,n in old.items() if j>=k) for k in range(high+1)]
            for k in range(3 if circuit else 1, high+1):
                tails[k]=min(tails[k-1],tails[k],limits.get(((i,v,c),k),0))
            out[i][v,c]={j:x-y for j,(x,y) in enumerate(zip(tails,tails[1:]+[0])) if x>y}
    return out


def physical_prices(v15, a, theta, ball, discounts, hist):
    y = a**6
    ordinary, marked = v15["v12"]["prices"](y,Q(1),theta,ball,discounts,hist)
    clean = lambda k: y**2 + Q(".528")/y**k
    pair = 384*clean(0)+72*clean(1)+24+20/y+12/y**2
    ordinary[1][0,2]=min(ordinary[1][0,2], discounts[0,2]*pair)
    return ordinary, marked


def check_pair_bounds():
    # Un centro t de solape u aporta el propio t y los extremos disjuntos
    # que alcanzan el producto v*u exigido. Hay como máximo 2(4-u) incidencias libres.
    for u in range(1, 5):
        at_least_one = 1 + 2*(4-u)//(2 if u == 1 else 1)
        at_least_two = 1 + (2 if u == 2 else 0)
        assert at_least_one <= 4*u and 2*at_least_two <= 3*u
    assert 72+4*8 == 104 and 3*8//2 == 12


def circuit_data(circuit):
    far,histories,contact,stripes,base,iso,prior=map(circuit["previous"],(
        "propagacion_lejana","medidas_por_historia","contacto_de_checks","franjas",
        "medidas_y_un_dato","medidas_aisladas","exceso_sin_compensar"))
    profiles,stable,reps=contact["geometry"]()
    wins=far["winners"](reps)
    _,cap=stripes["cost_tables"](stripes["light_profiles"](),circuit["previous"]("fondo_pasivo")["step_costs"],
                                 circuit["previous"]("pares_adyacentes")["stripe_profiles"])
    small=histories["data_coefficients"](contact["positions"](profiles,stable,cap,stripes["coefficients"](4)))
    triples=iso["count_triples"](profiles,stable,contact["rooted_trees"])
    trees,six=stripes["coefficients"](10),stripes["coefficients"](6)
    hist=prior["cumulative"](trees,small,triples,iso["group_series"](10),far["pair_counts"](wins,stable))
    prior["refine_tails"](hist,six,reps)
    geometry=tuple(sorted(far["cover_profiles"](reps)))
    v11=previous("verificar_H5A_3_4_ramificaciones_dos_ejes.py")
    discounts=v11["geometry_discounts"](geometry)
    triple=circuit["triple_counts"](iso,profiles,stable,contact["rooted_trees"])
    halos=contact["exposure"](reps,base["budgets"])
    return hist,trees,far["crowded_coefficients"](six),discounts,geometry,triple,halos


def partitions(n):
    """Particiones de n etiquetas, con cada bloque en orden de aparición."""
    out=[()]
    for i in range(n):
        new=[]
        for p in out:
            new.append(p+((i,),))
            for j in range(len(p)):
                new.append(p[:j]+(p[j]+(i,),)+p[j+1:])
        out=new
    return out


def tree_links(sizes, degree):
    if len(sizes)==1:
        return 1
    return sum(sizes)**(len(sizes)-2)*prod(sizes)*degree**(len(sizes)-1)


def shape_counts(n,kind):
    # Aristas de solape, de igual coordenada sin solape, de coordenadas diferentes.
    overlap,same,other=((9,54,576),(8,48,456))[kind]
    shapes=Counter()
    for blocks in partitions(n):
        for groups in partitions(len(blocks)):
            shape=tuple(sorted(tuple(sorted(len(blocks[i]) for i in g)) for g in groups))
            shapes[shape]+=1
    for shape in shapes:
        weight=tree_links([sum(g) for g in shape],other)
        for g in shape:
            weight*=tree_links(g,same)
            for k in g:
                weight*=tree_links([1]*k,overlap)
        shapes[shape]=Q(weight*shapes[shape],factorial(n-1))
    assert sum(shapes.values())==Q(n**(n-2)*(overlap+same+other)**(n-1),factorial(n-1))
    return shapes


def local_costs(kind, stable, representatives, axis):
    """Menor descenso por patrón de primeras extinciones en las líneas."""
    width=3+kind
    price={}
    for lit,foot in representatives:
        h,v=lit.bit_count(),foot.bit_count()
        subset=lit
        while subset:
            counts=tuple(sorted(k for stripe in axis if (k:=(stripe&subset).bit_count())))
            for rank in range(1,width+1):
                if h*width<rank*v or (rank==width and (h!=v or not stable[lit ^ subset])):
                    continue
                for largest in range(1,1+width*h//v):
                    key=rank,largest,counts
                    price[key]=min(price.get(key,1000),2*h-v)
            subset=(subset-1)&lit
    return price


@lru_cache(None)
def rank_grid(n,width):
    """Cada vector de rangos tiene un índice en base width+1. Restar un paso resta su índice."""
    from array import array
    states=list(product(range(width+1),repeat=n))
    ones=[(width+1)**j for j in reversed(range(n))]
    offsets=[]
    for r in states:
        sub=[0]
        for v,p in zip(r,ones):
            sub=[i+j*p for i in sub for j in range(v+1)]
        offsets.append(array("H",sub[1:]))
    return states,offsets


def rank_debits(n,kind,sizes,prices):
    width=3+kind
    states,offsets=rank_grid(n,width)
    blocks = []
    begin = 0
    for size in sizes:
        blocks.append(slice(begin, begin+size))
        begin += size
    ident={r:i for i,r in enumerate(states)}
    pr=[[1000]*len(states) for _ in range(width+1)]
    lines=[0]*len(states)
    for i,step in enumerate(states[1:],1):
        counts=tuple(sorted(x for block in blocks if (x:=sum(step[block]))))
        lines[i]=len(counts)
        for rank in range(1,width+1):
            pr[rank][i]=prices.get((rank,max(step),counts),1000)
    ans=[]
    for limit in (2,4):
        dp=[0]*len(states)
        for i,r in enumerate(states[1:],1):
            # Permutar los donantes de una línea conserva sus sumas y los máximos
            # de rango y de paso. El estado ordenado ya se resolvió.
            sorted_r=tuple(k for block in blocks for k in sorted(r[block]))
            if r!=sorted_r:
                dp[i]=dp[ident[sorted_r]]
                continue
            row=pr[max(r)]
            dp[i]=min(row[j]+dp[i-j] for j in offsets[i] if lines[j]<=limit)
            assert dp[i] < 1000
        ans.append(dp)
    return ans


def aligned_limits(hist, kind, nmax=5):
    co=previous("verificar_H5A_3_4_contacto_de_checks.py")
    far=previous("verificar_H5A_3_4_propagacion_lejana.py")
    _,stable,reps=co["geometry"]()
    axis=(far["ROWS"],far["COLS"])[kind]
    width = 3+kind
    prices = local_costs(kind, stable, reps, axis)
    ordinary, circuit = {}, {}
    for n in range(4,nmax+1):
        grid,_=rank_grid(n,width)
        shapes = shape_counts(n, kind)
        memo = {}
        counts=[Counter(),Counter()]
        for shape,weight in shapes.items():
            # Ordenar por tamaño los bloques que comparten coordenada.
            shape=tuple(sorted(shape,key=lambda g:(sum(g),g)))
            sizes=tuple(map(sum,shape))
            if sizes not in memo:
                memo[sizes]=rank_debits(n,kind,sizes,prices)
            dep=memo[sizes]
            pieces = []
            start = 0
            for g in shape:
                for sz in g:
                    pieces.append((sz, slice(start,start+sz)))
                    start += sz
            top=[0,0]
            for i,r in enumerate(grid):
                balance=0
                for sz,sl in pieces:
                    total=sum(r[sl])
                    if sz==1:
                        balance+=2*total-width
                    elif total>width*sz-2:
                        break
                    else:
                        balance+=total-(width*sz-total)%2
                else:
                    for j in (0,1):
                        top[j]=max(top[j],balance-dep[j][i])
            for j,t in enumerate(top):
                counts[j][t]+=weight
        key=(n,0) if kind==0 else (0,n)
        # Tres donantes de coordenadas fijas distintas en el mismo rectángulo.
        degree=(639,512)[kind]
        motif=(9000,2916)[kind]*Q((degree*n)**(n-3),factorial(n-3))
        for k in range(1,max(hist[kind][key])+1):
            sums=[sum(w for h,w in row.items() if h>=k) for row in counts]
            ordinary[(kind,*key),k]=min(sums[1],sums[0]+motif)
        # Una paridad de dato nula deja a lo sumo n-1 originales reales.
        proper=[0]*n
        for size in range(1,n):
            vals=[0]
            for j in range(2,size+1):
                small=(j,0) if kind==0 else (0,j)
                vals.append(max(hist[kind][small])+proper[size-j])
            proper[size]=max(vals)
        for k in range(proper[-1]+1,max(hist[kind][key])+1):
            circuit[(kind,*key),k]=ordinary[(kind,*key),k]
    return ordinary,circuit


def apply_caps(hist,limits):
    new=[dict(tab) for tab in hist]
    addressed = {key for key, _ in limits}
    for i,tab in enumerate(hist):
        for (v,c),old in tab.items():
            if (i,v,c) not in addressed:
                continue
            tails=[sum(w for j,w in old.items() if j>=k) for k in range(max(old)+1)]
            for k in range(1,len(tails)):
                tails[k]=min(tails[k-1],tails[k],limits.get(((i,v,c),k),tails[k]))
            assert all(a >= b >= 0 for a,b in zip(tails, tails[1:]))
            new[i][v,c]={k:a-b for k,(a,b) in enumerate(zip(tails,tails[1:]+[0])) if a>b}
    return new


PHEN = dict(name="Lector común", rates=("4.4e-7", ".001143", ".004", ".000777", "1.01"),
    m=("8.28593485203236e-7", "9.66004196291571e-7"),
    u=("5.82167417582085e-7", "6.23147766829749e-7"),
    h=("3.54793050333231e-8", "1.11269902992774e-7"),
    physical=("4.32075158412186e-9", "1.5313510711388e-8"),
    outer=(".001148", "4.87115611667827e-7", "2.02017754058527e-6"),
    d=(".00005484831685", ".00013117678238"),
    e=(".00005484831685", ".00036140734522"))

CIRCUITS = (
    ("VV primero", "8.3e-8", 4, ".56", ".74", ".001141", ".001148",
     (("1.18452883015478e-6", "1.35646782758077e-6"),
      ("3.83256569339477e-7", "7.1967748289793e-7"),
      ("5.57242251072595e-8", "2.37366511744496e-7"),
      ("6.13340780858609e-9", "3.32299420645426e-8"),
      "7.67818840335149e-7", "1.56586830888333e-6")),
    ("Columnas en cualquier orden", "7.5e-8", 7, ".557", ".76", ".001145", ".00115",
     (("1.0736060046244e-6", "1.2304204741898e-6"),
      ("5.4677513665222e-7", "6.58261024227923e-7"),
      ("8.51662065816268e-8", "2.57921177398537e-7"),
      ("8.69509924776902e-9", "3.68745410131238e-8"),
      "8.26420250603591e-7", "1.48446176671018e-6")))


def verify_phen(v15, ball, halos, series, hist, discounts):
    cert = PHEN
    a, s = Q(".57"), Q(".662")
    y = a**6
    reg, marked = physical_prices(v15, a, Q(cert["outer"][0]), ball, discounts, hist)
    marked[0][0, 0] = v15["single_tail"](tuple(map(Q, cert["m"])))
    v11 = v15["v11"]
    kernel = v15["corrected_kernel"](v11["prop"]["kernels"](ball, y, s), cert, reg, y/s**2)
    d, e = [tuple(map(Q, cert[field])) for field in ("d", "e")]
    effective = v11["effective"]
    effective["certify"](y, s, halos, reg, marked, kernel,
                         lambda w, i: effective["full_tail"](series, a, w, i, d, e), cert)


def verify_circuits(ball, halos, hist, trees, crowded, discounts, geometry, triple, far, anchor):
    for name, rate, first, a, s, b, theta, witness in CIRCUITS:
        c = previous("verificar_FASE5_circuito_sindrome.py")
        rate, a, s = map(Q, (rate, a, s))
        # Cada certificado carga su propia instancia del verificador previo.
        # b y theta son parámetros auxiliares libres del recuento por historias.
        params = c["certificate"].__globals__
        params.update(B_PASS=Q(b), THETA=Q(theta))
        assert c["laws"].__globals__ is params
        law = [c["laws"](0, rate, first), c["laws"](1, rate, 7)]
        reg, marked = c["polys"](hist, trees, crowded, discounts, geometry, triple, law,
                                 far["N"], anchor, ball, a**6)
        row, largest, tails, passive, active = c["certificate"](
            law, reg, marked, halos, far, ball, rate, a, s, witness)
        print(f"{name}, rc={float(rate):g}, rm=.001. Fila {float(row):.12g} sobre .002.")
        print(f"  Cociente máximo {float(largest):.12f}. Pasivo {float(passive):.12f}, "
              f"activo {float(active):.12f}.")
        print(f"  Colas VV, CC {float(tails[0]):.8g}, {float(tails[1]):.8g}.", flush=True)


def main():
    check_pair_bounds()
    limits = triple_limits()
    chain = 2916*6*64*6**14*Q(1,10**6)**3
    assert chain == Q(334731302496,3814697265625) < 1
    print(f"Tríos CC a saldo 4: 2916. Tríos VV a saldo 3: 81. Factor de (4): {float(chain):.12g}.",
          flush=True)
    v15 = previous("verificar_H5A_3_4_contacto_oculto.py")
    ball, halos, series, geometry, hist, _ = v15["v11"]["data_tables"]()
    discounts = v15["v11"]["geometry_discounts"](geometry)
    physical, circuit = {}, {}
    hist_p = refine(hist, limits)
    for kind in (0, 1):
        cp, cc = aligned_limits(hist_p, kind)
        physical.update(cp)
        circuit.update(cc)
    assert physical[((1,0,4),6)] == 58368 and physical[((1,0,5),8)] == 0
    verify_phen(v15, ball, halos, series, apply_caps(hist_p, physical), discounts)
    c = previous("verificar_FASE5_circuito_sindrome.py")
    # La fase 5 usa su propio histograma de pares mixtos y pesos por fibra.
    base_c, trees, crowded, discounts_c, geometry_c, triple, halos_c = circuit_data(c)
    far = v15["v11"]["prop"]
    for i, tab in enumerate(base_c):
        assert all(row == hist[i][key] for key,row in tab.items() if 0 in key)
    hist_c = apply_caps(refine(base_c, limits, circuit=True), circuit)
    verify_circuits(ball, halos_c, hist_c, trees, crowded, discounts_c, geometry_c,
                    triple, far, v15["v10"]["anchor"])


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para comprobar las desigualdades.")
    main()
