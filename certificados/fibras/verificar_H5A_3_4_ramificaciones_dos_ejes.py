"""HGP binario (3,4). Ramificaciones y saldo de los dos ejes, con aritmética racional."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import ceil, factorial
from pathlib import Path
from runpy import run_path

prior = run_path(str(Path(__file__).with_name("verificar_H5A_3_4_exceso_sin_compensar.py")))
old, prop, effective = prior["old"], prior["prop"], prior["effective"]
previous, dot, upper, N = prior["previous"], prior["dot"], prior["upper"], prior["N"]
Q_RATE = Q(".001") * upper(Q(".015")) / (1-Q(".002"))
Z_RATE = Q("1e-6") * upper(Q(".019")) / (1-Q(".000002"))
ETA = Q(".0001")

SPARSE = dict(a=".52", xi="1", refund=0, b=".001155", scale=".13", group=Q(215,13),
              m=("2.53336973855e-7","2.97759463946e-7"), u=("2.88660334492e-6","3.07735469227e-6"),
              h=("1.62463759649e-7","6.59358213099e-7"), physical=("6.40679611474e-9","3.93795762100e-8"),
              outer=(".00115463477719","1.44600459180e-7","3.70129362246e-6"),
              d=("2.77949717767e-5","7.03213333142e-5"),e=("2.77949717767e-5",".000237198777884"))
AXES = dict(a=".6375", xi=".15938", refund=1, b=".001155", scale="1", group=Q(13,10),
            m=("2.04777241373e-6","2.44246582035e-6"), u=("1.79721361528e-6","1.92456037426e-6"),
            h=("2.49810504392e-7","4.87747793335e-7"), physical=("1.49244362814e-8","3.54897931213e-8"),
            outer=(".00116812989263","3.95581779891e-7","1.40478785947e-5"),
            d=("5.56123618250e-5",".000113807397393"),e=("5.56123618250e-5",".000262954272589"))


def kernel(ball, first, r1, second, r2, Y, xi, refund):
    """L suma todas las capas. M reserva el saldo sobrante del enlace distinguido."""
    A, B = ball(first,r1), ball(second,r2)
    da, db = ball(first,r1+2)-A, ball(second,r2+2)-B
    beta = Y / xi**refund
    assert 0 < beta < 1 and 6*Y < 1 and 36*Y*xi < 1
    total = A*B + (da*B+db*A)*Y/(1-6*Y)
    total += da*db*xi*Y*(1+12*Y/(1-6*Y))/(1-36*xi*Y)
    partial = best = moment = A*B
    inner = Q(0)
    for length in range(1,200):
        shell = 6**(length-1)*(da*B+db*A+2*da*db*inner) + da*db*36**(length-1)*xi**length
        if length <= 4:
            row = [A]+[ball(first,r1+2*j)-ball(first,r1+2*j-2) for j in range(1,length+1)]
            col = [B]+[ball(second,r2+2*j)-ball(second,r2+2*j-2) for j in range(1,length+1)]
            explicit = sum(x*y*xi**min(i,j) for i,x in enumerate(row)
                           for j,y in enumerate(col) if max(i,j)==length)
            assert shell == explicit
        weight = shell*Y**length
        partial += weight
        moment = beta*moment+weight
        best = max(best,moment)
        assert 0 <= total-partial
        if beta*moment+total-partial <= best:
            return ceil(total),ceil(best)
        inner += xi*(6*xi)**(length-1)
    raise ArithmeticError("No se certificó la cola del enlace.")


def kernels(ball,Y,xi,refund):
    def k(v1,r1,v2,r2):
        return kernel(ball,v1,r1,v2,r2,Y,xi,refund)
    return (((k(1,4,1,6),k(1,5,1,5)),(k(0,5,0,5),k(0,6,0,4))),
            (k(1,3,1,4),k(0,4,0,3)),(k(0,3,1,4),k(0,4,1,3)),k(0,2,1,2))


def data_tables():
    contact,stripes,base,history,earlier = map(previous,("contacto_de_checks","franjas","medidas_y_un_dato",
                                                         "medidas_por_historia","medidas_aisladas"))
    profiles,stable,reps = contact["geometry"]()
    winners = prop["winners"](reps)
    assert min(phi.bit_count() for _,phi in winners) == 3
    assert all(phi.bit_count()<=5*(2*h.bit_count()-phi.bit_count()) for h,phi in winners)
    stripes_local = prop["ROWS"]+prop["COLS"]
    footprints = {0}
    for stripe in stripes_local:
        footprints |= {mask^stripe for mask in footprints}
    assert len(footprints)==64 and min(mask.bit_count() for mask in footprints if mask)==3
    # Dos contactos en distinta fila no activan una acción si no hay otra luz.
    assert all(stable[(1<<i)|(1<<j)] for i in range(12) for j in range(i) if i//3!=j//3)
    pairs=effective["pair_classes"](winners,stable,effective["effective_budgets"](winners,stable))
    _,cost=stripes["cost_tables"](stripes["light_profiles"](),previous("fondo_pasivo")["step_costs"],
                                 previous("pares_adyacentes")["stripe_profiles"])
    positions=contact["positions"](profiles,stable,cost,stripes["coefficients"](4))
    small=history["data_coefficients"](positions)
    triples=earlier["count_triples"](profiles,stable,contact["rooted_trees"])
    series=earlier["group_series"](10)
    trees=stripes["coefficients"](prop["LIMIT"])
    hist=prior["cumulative"](stripes["coefficients"](10),small,triples,series,pairs)
    prior["refine_tails"](hist,trees,reps)
    return (history["ball"],contact["exposure"](reps,base["budgets"]),series,
            prop["cover_profiles"](reps),hist,stripes["coefficients"])


def geometry_discounts(geometry):
    result={}
    for n in range(1,prop["LIMIT"]+1):
        for v in range(n+1):
            c=n-v
            saving=max(0,246*v+340*c-prop["group_budget"](v,c,geometry,n,guard=True))
            result[v,c]=effective["discount"](saving)
    return result


def polynomials(Y,ball,discounts,hist):
    regular,marked=[{},{}],[{},{}]
    for i in (0,1):
        for key,counts in hist[i].items():
            discount=discounts.get(key,1)
            regular[i][key]=discount*sum(mu*Y**(-k) for k,mu in counts.items())
            marked[i][key]=discount*sum(mu*prior["anchor"](ball,i,k) for k,mu in counts.items() if k>=2)
    return regular,marked


def derivative_bound(d,e):
    weak,strong=prop["W"],prop["B"]
    dp,ep=[4*x for x in d],[4*x for x in e]
    def rounded(x):
        return Q(ceil(x*10**28),10**28)
    for _ in range(40):
        x=[dot(row,dp) for row in weak]
        y=[dot(row,ep) for row in strong]
        dn=[rounded(d[i]*(1+x[i])+e[i]*y[i]) for i in (0,1)]
        en=[rounded(e[i]*(1+x[i]+y[i])) for i in (0,1)]
        assert all(a<=b for a,b in zip(dn+en,dp+ep))
        dp,ep=dn,en
    return dp


def fit_row(witness,tables,discounts):
    ball,halos,series,_,hist,_=tables
    a,xi,b,scale=(Q(witness[k]) for k in ("a","xi","b","scale"))
    Y,tau=a**6,witness["group"]
    regular,marked=polynomials(Y,ball,discounts,hist)
    dd,dg,gd,gg=kernels(ball,Y,xi,witness["refund"])
    m,u,h,H,d,e=(tuple(map(Q,witness[k])) for k in ("m","u","h","physical","d","e"))
    theta,visible,rooted=map(Q,witness["outer"])
    g=40*H[0]+27*H[1]+rooted
    R=max(dot(row,u) for row in N)
    I=280*u[0]+189*u[1]+588*g
    ell=89*b+I
    eps=upper(ell)*(upper(ell)-1)
    assert max(u)<ETA and eps<1 and R<Q(".004")
    z=scale*Z_RATE
    def log_bound(x):
        return x+x*x/(2*(1-x))
    share=Q(89,2)*log_bound(theta**2)
    guard=log_bound(theta)-share
    ordinary=b+2*share-log_bound(theta)
    assert guard>=Q(".001") and ordinary>=0
    ratios=[Q_RATE*(1+b)*(1+eps)/b,Q_RATE*(1-theta)*upper(99*guard+90*ordinary+2*I)/theta]
    for i in (0,1):
        costs=2*dot([row[i] for row in halos[1:]],(u[0],u[1],g))
        ratios.append(z*upper(ETA+2*R+guard*(504,628)[i]+ordinary*halos[0][i]+costs)/m[i])
    w=tuple(m[i]*upper(dot([x[0] for x in dd[i]],h)+dg[i][0]*visible) for i in (0,1))
    depart=tuple(dot([x[1] for x in dd[i]],h)+dg[i][1]*visible for i in (0,1))
    t=tuple(Q("1.01")*mi for mi in m)
    removed=[]
    multi=[]
    derivative=derivative_bound(d,e)
    for i in (0,1):
        tail=effective["full_tail"](series,a,w,i,d,e)
        branches=tau*(effective["evaluate"](regular[i],w)+tail)
        removed.append(tau*effective["evaluate"](regular[i],m))
        pointed=tau*sum(mu*w[0]**v*w[1]**c*(v*depart[0]+c*depart[1]) for (v,c),mu in regular[i].items())
        ten=sum(mu*(v+c)*w[0]**v*w[1]**c/a**nu for (v,c,nu),mu in series[i].items())
        assert ten<derivative[i]
        multi.append(pointed+tau*max(depart)*(derivative[i]-ten))
        tree=upper(dot(N[i],t))
        assert m[i]*tree<t[i]
        nonsimple=min(branches-removed[i],pointed+tau*tail)+tau*m[i]*(tree-1)
        initial=((640,784)[i],(441,513)[i],(1344,1596)[i])
        simple=tau*z*upper(ETA+2*R+b*(210,252)[i]+2*dot(initial,(u[0],u[1],g)))*(1+Q_RATE)**(3+i)
        ratios.extend((branches/h[i],nonsimple/H[i],(simple+H[i])/u[i]))
    J=dot([x[0] for x in gd],h)+gg[0]*visible
    one=dot([x[1] for x in gd],h)+gg[1]*visible
    base=dot([x[1] for x in gd],removed)
    direct=tau*sum(effective["evaluate"](tab,m) for tab in marked)
    # Las especies de un grupo tienen marca diagonal uno a cualquier orden.
    largest=max(m)
    rate=upper(Q(1))*1297*largest*Q(12,11)**9
    assert rate<1
    single_tail=2*tau*2000*11**9*Q(11**9,factorial(10))*largest*(1297*largest)**10/(1-rate)
    ratios.extend((Y*(upper(J)-1-J)/visible,(Y*(upper(J)*one-base)+direct+single_tail)/rooted))
    assert all(0<x<1 for x in ratios),[(i,float(x)) for i,x in enumerate(ratios) if not 0<x<1]
    factor,passive=prior["marked_sensitivity"](b,Q_RATE,u,H,g)
    speeds=[dot(row,t) for row in N]
    assert max(speeds)<1
    fixed=[]
    for i in (0,1):
        moment=m[i]*upper(speeds[i])*(1+speeds[i]/(1-max(speeds)))-m[i]
        single=tau*moment*(2*(R+24*I)+2160*passive)
        fixed.append((u[i]-H[i])*factor[i]+min(H[i],multi[i]+single))
    projected=max(dot(row,fixed) for row in N)
    assert projected<R
    print(f"Fila completa con marca {float(R):.12g}, proyección m>=2 {float(projected):.12g}, "
          f"cociente máximo {float(max(ratios)):.12f}.",flush=True)
    return projected


def unweighted_obstruction(ball,hist,discounts):
    """Las filas de ramas de este método, sin marcas, fallan en G3 para todo Y."""
    precision=10**14
    def below(x):
        return Q(x*precision//1,precision)
    def exp_lower(x):
        term=total=Q(1)
        for j in range(1,15):
            term*=min(x,Q(14))/j
            total+=term
        return total
    base=(Q("1.68e-6"),Q("1.90e-6"))
    for mi,guard in zip(base,(504,628)):
        assert mi < Q("1e-6")*exp_lower(Q(".019")+Q(guard,1000))
    # Una rama CC con e_C>=.1 ya no admite supersolución finita.
    assert base[1]*exp_lower(Q("12.9"))>Q(".1")
    weak,strong=prop["W"],prop["B"]
    steps=[]
    for cell in range(32):
        lo,hi=Q(cell,36*32),Q(cell+1,36*32)
        left,right=100000,400000
        while right-left>1:
            mid=(left+right)//2
            if hi*mid**6<=100000**6:
                left=mid
            else:
                right=mid
        inv=Q(left,100000)
        assert hi*inv**6<=1
        def link(spec):
            v1,r1,v2,r2=spec
            A,B=ball(v1,r1),ball(v2,r2)
            da,db=ball(v1,r1+2)-A,ball(v2,r2+2)-B
            return A*B+(da*B+db*A)*lo/(1-6*lo)+da*db*lo*(1+12*lo/(1-6*lo))/(1-36*lo)
        links=((link((1,4,1,6)),link((1,5,1,5))),(link((0,5,0,5)),link((0,6,0,4))))
        one=[lo*discounts[1-i,i] for i in (0,1)]
        two=[[sum(n*hi**(-k) for k,n in hist[i][2-i-j,i+j].items())*
              discounts[2-i-j,i+j] for j in (0,1)] for i in (0,1)]
        h=d=e=(Q(0),Q(0))
        for step in range(1,101):
            w=[below(base[i]*exp_lower(dot(links[i],h))) for i in (0,1)]
            X,Z=[dot(row,d) for row in weak],[dot(row,e) for row in strong]
            h=[below(w[i]*(one[i]+dot(two[i],w))) for i in (0,1)]
            d=[below(w[i]*exp_lower(X[i])*(inv**(7+i)+inv**(7+3*i)*(exp_lower(Z[i])-1))) for i in (0,1)]
            e=[below(w[i]*inv**(7+3*i)*exp_lower(X[i]+Z[i])) for i in (0,1)]
            if e[1]>=Q(".1"):
                break
        assert e[1]>=Q(".1"),cell
        steps.append(step)
    print(f"Sin marca diagonal ni inclinación de talla: 32 intervalos incompatibles, máximo {max(steps)} iteraciones.")


def three_arm_bound(ball, hist, degree15, discounts):
    """Un núcleo, tres grupos, bifurcación en un rectángulo lejano."""
    limit = 17  # n <= 9m-10 en el resto de v10, aquí m=3.
    cap = 5*limit//3
    precision = 10**75
    mass = (Q("1.846e-6"), Q("2.141e-6"))
    prior["single_group_response"](ball,hist,mass)
    data_to_gen = ((1,3,1,4),(0,4,0,3))
    gen_to_data = ((0,3,1,4),(0,4,1,3))
    tallies = [dict(row) for row in hist]
    for i in (0,1):
        for (v,c),count in degree15[i].items():
            if 11 <= v+c <= limit-2 and count:
                tallies[i][v,c] = {(7*v+10*c)//6: count}

    costs = {}
    for lengths in sorted(product(range(cap+1), repeat=3), key=sum, reverse=True):
        if 1+sum(lengths) > cap:
            continue
        center = min(lengths)+1
        levels = Counter({center:1})
        for length in lengths:
            for j in range(1,length+1):
                levels[min(j,length+1-j+center)] += 1
        need = 0
        for level in range(max(levels),0,-1):
            need += max(levels[level],(need+1)//2)
        for axis in range(3):
            longer = list(lengths)
            longer[axis] += 1
            need = min(need,costs.get(tuple(longer),cap+1))
        costs[lengths] = need

    @lru_cache(None)
    def shell(spec,length):
        v,r,w,t = spec
        # En la capa cero se excluyen los centros que leen el dato.
        return (ball(v,r+2*length)*ball(w,t+2*length)
                - ball(v,r+2*length-2)*ball(w,t+2*length-2))

    @lru_cache(None)
    def volume(specs):
        values = [0]*(cap+1)
        for (a,b,c),need in costs.items():
            if need <= cap:
                values[need] += shell(specs[0],a)*shell(specs[1],b)*shell(specs[2],c)
        for k in range(1,cap+1):
            values[k] += values[k-1]
        return values

    @lru_cache(None)
    def group(root,endpoint=None):
        totals = Counter()
        for (v,c),counts in tallies[root].items():
            n = v+c
            anchors = 1 if endpoint is None else (v,c)[endpoint]
            weight = mass[0]**v*mass[1]**c*anchors*discounts.get((v,c),1)
            for k,count in counts.items():
                delta = min(4,max(-1,k-n+1)) if n>=3 else 0
                if weight:
                    totals[n,k,delta] += weight*count
        return {key:ceil(precision*x) for key,x in totals.items()}

    @lru_cache(None)
    def children(j,k):
        pairs = Counter()
        for (n,a,dx),x in group(j).items():
            for (r,b,dy),y in group(k).items():
                if n+r<limit:
                    pairs[n+r,a+b,min(3,dx+dy)] += x*y
        return pairs

    by_root = []
    for root in (0,1):
        total = 0
        for i,j,k in product((0,1),repeat=3):
            vol = volume((data_to_gen[i],gen_to_data[j],gen_to_data[k]))
            for (n,a,dx),x in group(root,i).items():
                for (r,b,dy),y in children(j,k).items():
                    if n+r<=limit and a+b>=0 and dx+dy>=2:
                        total += x*y*vol[a+b]
        by_root.append(Q(total,2*precision**3))
    return max(dot(row,by_root) for row in N)


def main():
    tables=data_tables()
    discounts=geometry_discounts(tables[3])
    unweighted_obstruction(tables[0],tables[4],discounts)
    sparse=fit_row(SPARSE,tables,discounts)
    axes=fit_row(AXES,tables,discounts)
    prod=SPARSE["group"]*Q(SPARSE["scale"])
    assert prod**8*Q(SPARSE["scale"])**3>1
    tau=AXES["group"]
    assert tau**7*Q(AXES["xi"])>1
    covered=sparse+tau**4*axes
    assert sparse<Q(".000120006") and axes<Q(".000165643") and tau**4*axes<Q(".000473092")
    assert covered<Q(".000593097")
    assert Q(".002903264")+covered+Q(".0005036") < Q(".004")
    print(f"Parte cubierta de la fila pendiente {float(covered):.12g}.",flush=True)
    fork=three_arm_bound(tables[0],tables[4],tables[5](15),discounts)
    assert fork<Q(".000024016")
    print(f"Tres grupos en un núcleo con bifurcación lejana {float(fork):.12g}.",flush=True)


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para comprobar las desigualdades.")
    main()
