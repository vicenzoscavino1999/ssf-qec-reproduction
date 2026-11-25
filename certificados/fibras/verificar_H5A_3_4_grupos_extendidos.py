"""Grupos extendidos en G3 para HGP (3,4). Python 3.10+, biblioteca estándar."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import factorial
from pathlib import Path
from runpy import run_path


def previous(stem):
    return run_path(str(Path(__file__).with_name(f"verificar_H5A_3_4_{stem}.py")))


effective = previous("sindrome_efectivo")
prop = effective["propagation"]
N, dot, upper = prop["N"], prop["dot"], effective["upper"]
ROWS, COLS = prop["ROWS"], prop["COLS"]
GENERAL = ((2,3,6), (2,2,4,6))
SEPARATED = ((2,4,8), (2,2,7,9))
SINGLE_CELL = ((2,4,12), (2,2,8,12))


def check_credits(profiles, stable):
    """Una fila/columna física por fuente y ningún doble cobro de un primer apagado."""
    assert all(a<=b<=c for i in range(2) for a,b,c in zip(GENERAL[i],SEPARATED[i],SINGLE_CELL[i]))
    tested = [set(),set()]
    for lit,h,v,options in profiles:
        ranks = (3*h//v,4*h//v)
        assert h*max(max(row[:r]) for row,r in zip(GENERAL,ranks)) <= 6*(2*h-v)
        choices = [(i,mask,n,blocked) for i,opts in enumerate(options)
                   for mask,n,blocked in opts if mask]
        def place(start,used,blocked,counts,pristine):
            if counts and (not pristine or stable[lit & ~used]):
                for mode,fees in enumerate((SEPARATED,SINGLE_CELL)):
                    if mode and len(counts)>1 and any(n>1 for _,_,n in counts):
                        continue
                    total = sum(sum(fees[i][r-n:r]) for i,r,n in counts)
                    assert total <= 6*(2*h-v)
                    tested[mode].add((tuple(sorted(counts)),2*h-v))
            for j in range(start,len(choices)):
                i,mask,n,hit = choices[j]
                if mask & blocked:
                    continue
                for rank in range(n,ranks[i]+1):
                    place(j+1,used|mask,blocked|hit,counts+[(i,rank,n)],pristine or rank==3+i)
        place(0,0,0,[],False)
    for fees,expected in ((GENERAL,(7,10)),(SEPARATED,(4,4)),(SINGLE_CELL,(0,0))):
        for typ,row in enumerate(fees):
            signed = [6*(r if fees is GENERAL else 2*r-len(row))-sum(row[:r])
                      for r in range(len(row)+1)]
            assert max(signed) == expected[typ]
    assert all(tested)


def homogeneous_credits(profiles,stable):
    # Créditos en sextos. Se permiten lectores alineados de dos en dos.
    cases=(((3,3,9),(2,3,9)),((2,2,6,9),(3,3,6,9)))
    for typ,alternatives in enumerate(cases):
        for mode,fees in enumerate(alternatives):
            own_limit,transverse_limit=((2,1),(1,2))[mode]
            for lit,h,v,options in profiles:
                highest=(3+typ)*h//v
                for mask,n,_ in options[typ]:
                    if n:
                        for r in range(n,highest+1):
                            if r<3+typ or stable[lit & ~mask]:
                                assert sum(fees[r-n:r])<=6*(2*h-v)
                firsts=lit
                while firsts:
                    own=(ROWS,COLS)[typ]
                    cross=(COLS,ROWS)[typ]
                    if (all((firsts & line).bit_count()<=own_limit for line in own) and
                        all((firsts & line).bit_count()<=transverse_limit for line in cross)):
                        for ranks in combinations_with_replacement(range(1,highest+1),firsts.bit_count()):
                            if max(ranks)<3+typ or stable[lit & ~firsts]:
                                assert sum(fees[r-1] for r in ranks)<=6*(2*h-v)
                    firsts=(firsts-1) & lit
            assert max(6*(2*r-len(fees))-sum(fees[:r]) for r in range(len(fees)+1))==((3,4),(5,3))[typ][mode]


def previous_patch(mask,winners):
    """Descuento g(B) de v8 para auditar el ejemplo que cruza varios rectángulos."""
    cuts = set()
    for bits in range(1,128):
        phi = 0
        for j,line in enumerate(ROWS+COLS):
            if bits >> j & 1:
                phi ^= line
        if phi and 2*phi.bit_count() >= 3*bits.bit_count():
            cuts.add(phi)
    @lru_cache(None)
    def excess(pending):
        if not pending:
            return 0
        rate = max(Q((pending & c).bit_count(),c.bit_count()) for c in cuts)
        firsts, best = pending, 10000
        while firsts:
            cost = min((3*(2*h.bit_count()-phi.bit_count())-h.bit_count() for h,phi in winners
                        if h & firsts == firsts and h.bit_count() >= rate*phi.bit_count()),default=10000)
            best = min(best,cost+excess(pending ^ firsts))
            firsts = (firsts-1) & pending
        return best
    pending,best = mask,5*mask.bit_count()
    while pending:
        best = min(best,5*(mask ^ pending).bit_count()+excess(pending))
        pending = (pending-1) & mask
    return best


def check_example(winners):
    points = set(product(range(4),repeat=2))
    lines = {(kind,t):{(x,y) for x,y in points if (y,x,(y-x)%4)[kind]==t}
             for kind,t in product(range(3),range(4))}
    incident = {p:{l for l,pts in lines.items() if p in pts} for p in points}
    assert all(len(L)==4 for L in lines.values()) and all(len(I)==3 for I in incident.values())
    sources = [(lines[0,i],(i,0)) for i in range(3)]
    stripes = [{(p,c) for p in line} for line,c in sources]
    sigma = set.union(*stripes)
    assert len(sigma)==12
    common = False
    for line,cols in product(lines.values(),incident.values()):
        rect = set(product(line,cols))
        assert not sigma <= rect
        reads = [stripe & rect for stripe in stripes if stripe & rect]
        if len(reads)>1:
            assert all(len(block)==1 for block in reads)
            hits = set.union(*reads)
            assert len({p for p,c in hits}) == len({c for p,c in hits}) == len(reads)
            common |= len(reads)==3
    assert common
    for j,value in enumerate((0,0,0,3,8)):
        assert previous_patch(sum(1<<(3*i) for i in range(j)),winners)==value
    for j in range(4):
        assert previous_patch(sum(1<<(4*i) for i in range(j)),winners)==0
    for mask in (1+8,1+8+128,1+16+64):
        assert previous_patch(mask,winners)==0
    assert (24-(3*3+(8-3)))//3 == 3 > 2
    # Dos columnas distintas, una de ellas compartida por dos de los tres CC.
    repeated=[{(p,c) for p in lines[0,i]} for i,c in enumerate(((0,0),(0,0),(1,0)))]
    assert len(set.union(*repeated))==12
    common=False
    for line,cols in product(lines.values(),incident.values()):
        rect=set(product(line,cols))
        reads=[stripe & rect for stripe in repeated if stripe & rect]
        if len(reads)>1:
            assert all(len(block)==1 for block in reads)
            hits=set.union(*reads)
            assert len({p for p,c in hits})==len(hits)
            assert max(Counter(c for p,c in hits).values())<=2
            common |= len(reads)==3
    assert common and 3*5//6==2


def line_checks(profiles,stable,reps,old):
    """Descenso exterior y créditos en una sola fila o columna, incluidos los empates."""
    for lit,phi in reps:
        for line in ROWS+COLS:
            assert 2*(lit & ~line).bit_count() >= (phi & ~line).bit_count()
    for typ,fees in enumerate(((3,3,6),(3,3,4,6))):
        for lit,h,v,options in profiles:
            rank_max=(3+typ)*h//v
            for line in (ROWS,COLS)[typ]:
                choices=[(mask,n) for mask,n,_ in options[typ] if mask and (mask & line)==mask]
                def place(start,used,credit,pristine):
                    if not pristine or stable[lit & ~used]:
                        assert credit <= 6*(2*h-v)
                    for j in range(start,len(choices)):
                        mask,n=choices[j]
                        if used & mask:
                            continue
                        for r in range(n,rank_max+1):
                            place(j+1,used|mask,credit+sum(fees[r-n:r]),pristine or r==3+typ)
                place(0,0,0,False)
        assert max(6*r-sum(fees[:r]) for r in range(len(fees)+1)) == (6,8)[typ]
        # Clases del trío en la línea: sin arista, disjuntas, con solapamiento.
        for edges in combinations_with_replacement((0,3,4),3):
            if edges.count(0)<=1:
                assert old["triple_budget"]((typ,)*3,edges,profiles,stable)<=3
        for edges in combinations_with_replacement((0,2),3):
            if edges.count(0)<=1:
                assert old["triple_budget"]((typ,)*3,edges,profiles,stable)<=(3,4)[typ]


def stripe_anchor(ball,i,k):
    if k<=0:
        return 0
    depth,spent=0,1
    while (3*spent+1)//2 <= k:
        spent=(3*spent+1)//2
        depth+=1
    return 4*ball(1,4+2*depth) if i==0 else 3*ball(0,4+2*depth)


FIBER_DEGREES=((63,9),(56,8))
FIBER_SCALES=((50,96),(185,679))


def compensated_polynomials(Y,ball,trees,pairs,tables,triples):
    branches,physical = effective["sector_polynomials"](Y,ball,trees,pairs)
    for i in range(2):
        anchor = lambda k: prop["anchor"](ball,i,k)
        # Coeficientes con raíz de dato fija y con saldo suficiente para la capa física.
        by_history = Counter()
        by_triple = Counter()
        for (nv,nc,_,cap),count in tables[i].items():
            for layer in range(1,cap+1):
                by_history[nv,nc,layer] += count
        for (nv,nc,cap),count in triples[i].items():
            for layer in range(1,cap+1):
                by_triple[nv,nc,layer] += count
        for key,multiplicity in trees[i].items():
            n=sum(key)
            if n<3:
                continue
            threshold=2*n//3
            branches[i][key]=multiplicity/Y**threshold
            physical[i][key] = 0
            nv,nc=key
            for layer in range(1,(7*nv+10*nc)//6+1):
                count=multiplicity
                if n in (3,4):
                    count=min(count,by_history[nv,nc,layer])
                if n==3:
                    count=min(count,by_triple[nv,nc,layer])
                shell=anchor(layer)-anchor(layer-1)
                weight=min(1/Y,Y**(layer-threshold))
                physical[i][key]+=count*weight*shell
            if key[1-i]:
                continue
            for direction,degree in enumerate(FIBER_DEGREES[i]):
                if direction==0:
                    cap=3 if n==3 else (6,8)[i]*n//6
                    shell_extra=stripe_anchor(ball,i,cap)-stripe_anchor(ball,i,threshold)
                else:
                    cap=(3,4)[i] if n==3 else (7,10)[i]*n//6
                    shell_extra=anchor(cap)-anchor(threshold)
                rooted=Q(n**(n-2)*degree**(n-1),factorial(n-1))
                branches[i][key]+=rooted*max(0,Y**(-cap)-Y**(-threshold))
                physical[i][key]+=rooted*max(0,shell_extra)
    return branches,physical


def bounded_tail(trees,Y,w,i,d):
    scale = Q("13.58")
    assert scale**3 >= Y**-2
    total = scale*w[i]*upper(dot(N[i],d))
    assert total < d[i]
    tail = total-effective["evaluate"](trees[i],tuple(scale*t for t in w))
    assert tail > 0
    for direction,(deg,price) in enumerate(zip(FIBER_DEGREES[i],FIBER_SCALES[i])):
        sixths=((6,7),(8,10))[i][direction]
        assert price**6 >= Y**(-sixths)
        root=price*w[i]
        ratio=deg*root*upper(Q(1))
        assert ratio<1
        tail+=Q(11**9,factorial(10))*root*(deg*root)**10/(1-ratio)
    return tail


COMPENSATED = dict(m=("2.28542702745e-6","2.81202337324e-6"), theta=".00119700171114",
                   u=("1.47878926788e-6","1.62198699624e-6"), h=("4.03462971383e-7","9.19699124856e-7"),
                   physical=("3.17519371738e-8","7.23500889204e-8"),
                   visible="4.52171252808e-7", rooted="3.15301244770e-5",
                   d=("3.54857049042e-5","4.47805789899e-5"))
DENSITY = dict(m=("7.5231347255e-6","8.7732310928e-6"), theta=".00114755089158",
               u=("2.23635214431e-11","2.92329097529e-11"), h=("1.51970716754e-11","5.26916427477e-11"),
               physical=("2.35654737323e-12","8.18460676878e-12"), visible="2.12565207415e-9",
               rooted="1.10141944786e-6", d=(".000186999217572",".000406866605892"),
               e=(".000186999217572",".000786664754115"))


def certify(witness,Y,s,ball,halos,poly,pure,make_tail,*,scale,group,single_group_factor,tree_factor,b):
    """Fijación y marcas físicas de v7/v8, con pesos por dato y por grupo explícitos."""
    p,q,alpha,eta=map(Q,("1e-6",".001",".004",".0001"))
    m,u,h,physical=(tuple(map(Q,witness[k])) for k in ("m","u","h","physical"))
    theta,visible,rooted=(Q(witness[k]) for k in ("theta","visible","rooted"))
    dd,dg,gd,gg=prop["kernels"](ball,Y,s)
    g=40*physical[0]+27*physical[1]+rooted
    row_mass=max(dot(row,u) for row in N)
    I=280*u[0]+189*u[1]+588*g
    ell=89*b+I
    eps=upper(ell)*(upper(ell)-1)
    y=q*upper(Q(".015"))/(1-2*q)
    z=scale*p*upper(Q(".019"))/(1-2*p)
    assert eps<1 and max(u)<eta and row_mass<alpha
    ratios=[y*(1+b)*(1+eps)/b]
    def log_bound(t):
        return t+t*t/(2*(1-t))
    c=Q(89,2)*log_bound(theta**2)
    guard,ordinary=log_bound(theta)-c,b+2*c-log_bound(theta)
    assert guard>=q and ordinary>=0
    ratios.append(y*(1-theta)*upper(99*guard+90*ordinary+2*I)/theta)
    for i in range(2):
        interaction=2*dot([row[i] for row in halos[1:]],(u[0],u[1],g))
        ratios.append(z*upper(eta+2*row_mass+guard*(504,628)[i]+ordinary*halos[0][i]+interaction)/m[i])
    w=tuple(m[i]*upper(dot([v[0] for v in dd[i]],h)+dg[i][0]*visible) for i in range(2))
    depart=tuple(dot([v[1] for v in dd[i]],h)+dg[i][1]*visible for i in range(2))
    trees=tuple(mi*tree_factor for mi in m)
    removed=[]
    for i in range(2):
        tail=make_tail(w,i)
        branches=group*(effective["evaluate"](poly[i],w)+tail)
        removed.append(group*effective["evaluate"](poly[i],m))
        pointed=group*sum(v*w[0]**nv*w[1]**nc*(nv*depart[0]+nc*depart[1]) for (nv,nc),v in poly[i].items())
        connected=upper(dot(N[i],trees))
        assert m[i]*connected<trees[i]
        two_or_more=m[i]*(connected-1)
        three_or_more=two_or_more-m[i]*dot(N[i],m)
        assert three_or_more>0
        fixed=min(branches-removed[i],pointed+group*tail)
        fixed+=group*(two_or_more+(single_group_factor-1)*three_or_more)
        initial=(210,252),(640,784),(441,513),(1344,1596)
        simple=group*z*upper(eta+2*row_mass+dot((b,2*u[0],2*u[1],2*g),
                                                [row[i] for row in initial]))*(1+y)**(3+i)
        ratios.extend((branches/h[i],fixed/physical[i],(simple+physical[i])/u[i]))
    J=dot([v[0] for v in gd],h)+gg[0]*visible
    one=dot([v[1] for v in gd],h)+gg[1]*visible
    without=dot([v[1] for v in gd],removed)
    direct=group*sum(effective["evaluate"](tab,m) for tab in pure)
    outer=Y/s**2
    ratios.extend((outer*(upper(J)-1-J)/visible,(outer*(upper(J)*one-without)+direct)/rooted))
    assert all(0<t<1 for t in ratios),[(j,float(t)) for j,t in enumerate(ratios) if not 0<t<1]
    # Se descarta el monomio individual de un núcleo simple, antes de sustituir variables físicas.
    rate=Q(".019")
    z_lower=group*scale*p*(1+rate+rate**2/2)/(1-2*p)
    residual=[v-z_lower*(1+q)**(3+i) for i,v in enumerate(u)]
    assert all(v>0 for v in residual)
    projected=max(dot(row,residual) for row in N)
    assert 0 < projected < row_mass
    print(f"Fila auxiliar {float(row_mass):.12g}, proyección {float(projected):.12g}, "
          f"mayor cociente exterior {float(max(ratios)):.12f}.")
    return row_mass,projected


def main():
    contacts,striped,base,history,old=map(previous,("contacto_de_checks","franjas","medidas_y_un_dato",
                                                     "medidas_por_historia","medidas_aisladas"))
    profiles,stable,reps=contacts["geometry"]()
    winners=prop["winners"](reps)
    check_credits(profiles,stable)
    homogeneous_credits(profiles,stable)
    line_checks(profiles,stable,reps,old)
    check_example(winners)
    caps=effective["effective_budgets"](winners,stable)
    pairs=effective["pair_classes"](winners,stable,caps)
    ball=history["ball"]
    halos=contacts["exposure"](reps,base["budgets"])
    trees=striped["coefficients"](10)
    _,prior=striped["cost_tables"](striped["light_profiles"](),previous("fondo_pasivo")["step_costs"],
                                   previous("pares_adyacentes")["stripe_profiles"])
    hist=contacts["positions"](profiles,stable,prior,striped["coefficients"](4))
    tables=history["data_coefficients"](hist)
    triples=old["count_triples"](profiles,stable,contacts["rooted_trees"])
    Y=Q(1,50)
    poly,pure=compensated_polynomials(Y,ball,trees,pairs,tables,triples)
    d=tuple(map(Q,COMPENSATED["d"]))
    complete_first,first=certify(COMPENSATED,Y,Q(".9999"),ball,halos,poly,pure,lambda w,i: bounded_tail(trees,Y,w,i,d),
                  scale=1,group=Q(1),single_group_factor=50,tree_factor=Q("1.006"),b=Q(".001174"))
    assert first<Q(".000666")
    series=old["group_series"](10)
    a=Q(2,3)
    poly,pure=prop["group_polynomials"](a**6,ball,tables,triples,series,pairs,
                                       prop["cover_profiles"](reps),striped["coefficients"](prop["LIMIT"]))
    d,e=(tuple(map(Q,DENSITY[k])) for k in ("d","e"))
    _,second=certify(DENSITY,a**6,Q(".3"),ball,halos,poly,pure,lambda w,i: effective["full_tail"](series,a,w,i,d,e),
                   scale=4,group=Q(1,4**9),single_group_factor=1,tree_factor=Q("1.06"),b=Q(".00117"))
    dense=4**8*second
    print(f"Parte cubierta del complemento <= {float(first):.12g} + {float(dense):.12g} "
          f"= {float(first+dense):.12g}.")
    assert dense < Q(".000806") and first+dense < Q(".001471")
    assert Q(".001471")+Q(".000537")==Q(".002008")
    assert complete_first+dense+Q(".001203") < Q(".004")
    print(f"Mayorante de la parte cubierta de la fila completa: {float(complete_first+dense):.12g}. "
          "Para alcanzar G3 basta que el mismo resto sea <= 0.001203.")
    # Si todos los CC pagan el precio bruto y tienen 513 continuaciones, ni su fila escalar cierra.
    assert 670**3 < 50**5 and 513*Q("1e-6")*(1+Q(".628"))*670*Q(5,2)>1
    print("G3 sigue abierto. Queda la clase con exceso positivo y n <= 9m-9.")


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para verificar las desigualdades.")
    main()
