"""Fila conjunta HGP (3,4). Paridad de las huellas y dos componentes originales."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
from pathlib import Path
from runpy import run_path

old = run_path(str(Path(__file__).with_name("verificar_H5A_3_4_grupos_extendidos.py")))
previous, prop = old["previous"], old["prop"]
N, dot, upper = old["N"], old["dot"], old["upper"]
DEGREES = old["FIBER_DEGREES"]
SCALES = ((Q("46.51163"),Q("88.20404")),(Q("167.26898"),Q("601.54655")))
BONUS = 6
effective = old["effective"]


def depth_from_cover(k):
    need, depth = 2, 0
    while need <= k:
        depth += 1
        need = (3*need+1)//2
    return depth


def anchor(ball,i,k,fiber=False):
    depth = depth_from_cover(k)
    if not depth:
        return 0
    extra=2*(depth-1)
    first=ball(0,(3,4)[i]+extra)
    second=ball(1,(4,3)[i]+extra)
    if fiber:
        return 4*second if i==0 else 3*first
    return first*second


def cumulative(trees,tables,triples,series,pairs):
    counts=[{},{}]
    for i in (0,1):
        for (v,c),total in trees[i].items():
            n=v+c
            if not total:
                continue
            if n==1:
                counts[i][v,c]={-1:total}
                continue
            tails=[total]
            for k in range(1,(7*v+10*c)//6+1):
                estimates=[total]
                if n==2:
                    pair=pairs[i][c-i]
                    estimates.append((pair[1]+2*pair[2],pair[2],0)[k-1])
                else:
                    estimates.append(sum(mu for (vv,cc,nu),mu in series[i].items()
                                         if (vv,cc)==(v,c) and nu//6>=k))
                    if n in (3,4):
                        estimates.append(sum(mu for (vv,cc,_,cap),mu in tables[i].items()
                                             if (vv,cc)==(v,c) and cap>=k))
                    if n==3:
                        estimates.append(sum(mu for (vv,cc,cap),mu in triples[i].items()
                                             if (vv,cc)==(v,c) and cap>=k))
                tails.append(min(estimates))
            assert tails==sorted(tails,reverse=True)
            counts[i][v,c]={k:x-y for k,(x,y) in enumerate(zip(tails,tails[1:]+[0])) if x>y}
    return counts


def refine_tails(tables,trees,reps):
    geometry=prop["cover_profiles"](reps)
    crowded=prop["crowded_coefficients"](trees)
    budgets={}
    for i,tab in enumerate(tables):
        for key,original in tab.items():
            n=sum(key)
            if not 3<=n<=prop["LIMIT"]:
                continue
            if key not in budgets:
                budgets[key]={j:prop["group_budget"](*key,geometry,j) for j in range(2,n+1)}
            caps=budgets[key]
            assert list(caps.values())==sorted(caps.values())
            tails=[sum(original.values())]
            for k in range(1,max(original)+1):
                count=sum(v for j,v in original.items() if j>=k)
                if k>caps[n]:
                    count=0
                else:
                    for j in range(3,n+1):
                        if k>caps[j-1]:
                            count=min(count,crowded[j][i].get(key,0))
                tails.append(count)
            assert tails==sorted(tails,reverse=True)
            tab[key]={j:x-y for j,(x,y) in enumerate(zip(tails,tails[1:]+[0])) if x>y}


def sector(Y,ball,hist,pairs):
    poly,pure=[Counter(),Counter()],[Counter(),Counter()]
    evaluate=effective["discount"]
    for i in (0,1):
        poly[i][1-i,i]=Y*evaluate((108,64)[i])
        for j in (0,1):
            key=(2-i-j,i+j)
            saving=evaluate(((12,118),(118,74))[i][j])
            for target,price in ((poly[i],lambda k:Y**(-k)),(pure[i],lambda k:anchor(ball,i,k))):
                assert 2*price(0)<=price(1) and 2*price(1)<=price(2)
                target[key]=sum(n*price(k) for k,n in enumerate(pairs[i][j]))*saving
        for key,by_cap in hist[i].items():
            n=sum(key)
            if n<3:
                continue
            threshold=n-1
            radius=lambda k:Y**max(-1,k-threshold)*(BONUS if k>threshold else 1)
            poly[i][key]=sum(mu*radius(k)*Y**(-k) for k,mu in by_cap.items())
            # Suma de capas de posiciones, cada una con su cota acumulada.
            pure[i][key]=sum(sum(mu for k,mu in by_cap.items() if k>=l)*radius(l)*
                             (anchor(ball,i,l)-anchor(ball,i,l-1)) for l in range(1,max(by_cap)+1))
            if key[1-i]:
                continue
            for axis,degree in enumerate(DEGREES[i]):
                cap=((3 if n==3 else (6,8)[i]*n//6) if axis==0 else
                     ((3,4)[i] if n==3 else (7,10)[i]*n//6))
                if cap>threshold:
                    count=Q(n**(n-2)*degree**(n-1),factorial(n-1))
                    poly[i][key]+=count*(Y**(-cap)-BONUS*Y**(-threshold))
                    pure[i][key]+=count*(anchor(ball,i,cap,axis==0)-anchor(ball,i,threshold,axis==0))
    return poly,pure


WITNESS=dict(m=("1.98480016535e-6","2.34735845908e-6"),theta=".00116039746203",
             u=("1.37391285409e-6","1.47781230869e-6"),h=("3.39821561052e-7","7.24712905821e-7"),
             physical=("2.47381416744e-8","5.18330475831e-8"),visible="4.98220098651e-7",
             rooted="8.40219933573e-6",d=(".000116216399417",".000143659726530"))


def tail(trees,Y,w,i,d):
    bulk=(w[i]/Y)*upper(dot(N[i],d))
    assert bulk<d[i]
    first_ten=effective["evaluate"](trees[i],tuple(x/Y for x in w))
    assert bulk>first_ten
    result=BONUS*Y*(bulk-first_ten)
    for j,degree in enumerate(DEGREES[i]):
        price=SCALES[i][j]
        assert price**6>=Y**(-((6,7),(8,10))[i][j])
        root=w[i]*price
        rate=upper(Q(1))*degree*root
        assert rate<1
        result+=Q(11**9,factorial(10))*root*(degree*root)**10/(1-rate)
    return result


@lru_cache(None)
def distance_budget(length):
    # En una cadena más corta entre dos grupos hay 1 ó 2 ocurrencias por capa.
    total=0
    for depth in range((length+1)//2,0,-1):
        on_path=1 if 2*depth-1==length else 2
        total+=max(on_path,(total+1)//2)
    return total


def link_volume(ball,i,j,budget):
    if budget<0:
        return 0
    length=0
    while distance_budget(length+1)<=budget:
        length+=1
    r=((4,5),(5,6))[i][j]
    return ball(1-i,r+2*length)*ball(1-i,10-r+2*length)


def marked_sensitivity(b,y,u,H,g):
    """Diferencia respecto del sistema sin núcleos libres, en la norma multivariada."""
    R=max(dot(row,u) for row in N)
    RH=max(dot(row,H) for row in N)
    I=280*u[0]+189*u[1]+588*g
    ell=89*b+I
    eps=upper(ell)*(upper(ell)-1)
    assert eps<1
    response=y*upper(2*ell)*(1+eps)
    contraction=89*response*(1+b)**2
    assert contraction<1
    D=response*I/(1-contraction)
    per=(2+b)*(1+b)*D
    rows=[]
    for i in (0,1):
        initial=((640,784)[i],(441,513)[i],(1344,1596)[i])
        full=R+dot(initial,(u[0],u[1],g))
        non=RH+dot(initial,(H[0],H[1],g))
        qbox=(210,252)[i]
        assert upper(Q(".0001")+2*full+b*qbox)<2
        rows.append(full+non+qbox*per)
    return rows,per


def single_group_response(ball,hist,m):
    y=Q(".001")*upper(Q(".015"))/(1-Q(".002"))
    z=Q("1e-6")*upper(Q(".019"))/(1-Q(".000002"))
    b,theta=Q(".001139"),Q(".001145")
    eta=Q(".0001")
    ell=89*b
    weak=upper(ell)*(upper(ell)-1)
    guard=lambda t:t+t*t/(2*(1-t))
    c=Q(89,2)*guard(theta**2)
    A,B=guard(theta)-c,b+2*c-guard(theta)
    assert A>Q(".001") and B>0 and weak<1 and y*(1+b)*(1+weak)<b
    assert y*(1-theta)*upper(99*A+90*B)<theta
    simple=[]
    for i in (0,1):
        assert z*upper(eta+A*(504,628)[i]+B*(414,536)[i])<m[i]
        simple.append(z*upper(eta+b*(210,252)[i])*(1+y)**(3+i))
    tree=[Q("1.01")*x for x in m]
    T=[dot(row,tree) for row in N]
    assert max(T)<1 and all(m[i]*upper(T[i])<tree[i] for i in (0,1))
    h=[m[i]*(upper(T[i])-1) for i in (0,1)]
    moments=[m[i]*upper(T[i])*(1+T[i]/(1-max(T)))-m[i] for i in (0,1)]
    u=[simple[i]+h[i] for i in (0,1)]
    R=max(dot(row,u) for row in N)
    marked=0
    for i in (0,1):
        for (v,c),by_cap in hist[i].items():
            marked+=m[0]**v*m[1]**c*sum(mu*anchor(ball,i,k) for k,mu in by_cap.items())
    # Para una componente desde 11 datos, f_i(k)<=2000 n^9 y árbol de grado 1297.
    assert Q(3,2)**9>36 and max(Q(4*9*6**4,25),Q(8*3*6**4,25))<=2000
    mx=max(m)
    ratio=upper(Q(1))*1297*mx*Q(12,11)**9
    assert ratio<1
    marked+=2*2000*11**9*Q(11**9,factorial(10))*mx*(1297*mx)**10/(1-ratio)
    g=40*h[0]+27*h[1]+marked
    I=280*u[0]+189*u[1]+588*g
    simple_factor,per_check=marked_sensitivity(b,y,u,h,g)
    response=[]
    for i in (0,1):
        response.append(simple[i]*simple_factor[i]
                        +moments[i]*(2*(R+24*I)+2160*per_check))
    return max(dot(row,response) for row in N),R


def direct_two(ball,hist,m):
    # Ambos grupos con raíz de tipo específico. En los grados altos basta el árbol unicolor mayorante.
    stop=26
    maxm=max(m)
    tabs=[[],[]]
    for i in (0,1):
        for (v,c),by_cap in hist[i].items():
            for k,mu in by_cap.items():
                tabs[i].append((v+c,k,(v,c),mu*m[0]**v*m[1]**c))
        for n in range(11,stop):
            weight=Q(n**(n-2)*1297**(n-1),factorial(n-1))*maxm**n
            # n posibles extremos de cada tipo dan una cota para ambos.
            tabs[i].append((n,5*n//3,(n,n),weight))
    volume=lru_cache(None)(lambda i,j,k:link_volume(ball,i,j,k))
    rows=[Q(0),Q(0)]
    hard=[Q(0),Q(0)]
    for root in (0,1):
        for size,k,ends,wt in tabs[root]:
            for j in (0,1):
                for size2,l,_,wt2 in tabs[j]:
                    if size+size2>=stop or max(size,size2)<3:
                        continue
                    count=sum(ends[i]*volume(i,j,k+l) for i in (0,1))
                    value=wt*wt2*count
                    rows[root]+=value
                    excess=(max(0,k-size+1) if size>=3 else 0)+(max(0,l-size2+1) if size2>=3 else 0)
                    if size+size2<=8 and excess>=2:
                        hard[root]+=value
    # Desde 26: longitud <= 2+2 log_(3/2)(5n/3); volumen <= 50000*36^2*(5n/3)^18.
    assert max(Q(d1*d2*6**e,25) for d1,d2,e in ((9,9,5),(3,3,6),(4,4,6),(8,8,5)))<50000
    rate=upper(Q(1))*1297*maxm
    ratio=rate*Q(stop+1,stop)**20
    assert ratio<1
    rest=2*stop**2*maxm**2*rate**(stop-2)*50000*36**2*Q(5*stop,3)**18/(1-ratio)
    return max(dot(row,[x+rest for x in rows]) for row in N),max(dot(row,hard) for row in N)


def density_projection(ball,regular,series):
    """Serie auxiliar de densidad de v9, restringida a dos o más componentes."""
    from math import ceil
    witness=old["DENSITY"]
    m,u,h,d,e=(tuple(map(Q,witness[k])) for k in ("m","u","h","d","e"))
    H=tuple(map(Q,witness["physical"]))
    visible,rooted=map(Q,(witness["visible"],witness["rooted"]))
    a=Q(2,3); b=Q(".00117"); tau=Q(1,4**9)
    dd,dg,_,_=prop["kernels"](ball,a**6,Q(".3"))
    w=tuple(m[i]*upper(dot([v[0] for v in dd[i]],h)+dg[i][0]*visible) for i in (0,1))
    depart=tuple(dot([v[1] for v in dd[i]],h)+dg[i][1]*visible for i in (0,1))
    W,B=prop["W"],prop["B"]
    def next_derivative(dp,ep):
        return ([d[i]*(1+dot(W[i],dp))+e[i]*dot(B[i],ep) for i in (0,1)],
                [e[i]*(1+dot(W[i],dp)+dot(B[i],ep)) for i in (0,1)])
    dp,ep=[4*x for x in d],[4*x for x in e]
    for _ in range(40):
        f,g=next_derivative(dp,ep)
        round_up=lambda v:Q(ceil(v*10**28),10**28)
        nxt,pnt=[round_up(v) for v in f],[round_up(v) for v in g]
        assert all(t<=s for t,s in zip(nxt+pnt,dp+ep))
        dp,ep=nxt,pnt
    g=40*H[0]+27*H[1]+rooted
    I=280*u[0]+189*u[1]+588*g
    R=max(dot(row,u) for row in N)
    y=Q(".001")*upper(Q(".015"))/(1-Q(".002"))
    simple_factor,per=marked_sensitivity(b,y,u,H,g)
    roots=[]
    tree=[Q("1.06")*mi for mi in m]
    speeds=[dot(row,tree) for row in N]
    assert max(speeds)<1 and all(m[i]*upper(speeds[i])<tree[i] for i in (0,1))
    for i in (0,1):
        p1=sum(mu*(v+c)*w[0]**v*w[1]**c/a**nu for (v,c,nu),mu in series[i].items())
        assert dp[i]>p1
        distinguished=sum(mu*w[0]**v*w[1]**c*(v*depart[0]+c*depart[1])
                          for (v,c),mu in regular[i].items())
        multi=tau*(distinguished+max(depart)*(dp[i]-p1))
        simple=(u[i]-H[i])*simple_factor[i]
        moment=m[i]*upper(speeds[i])*(1+speeds[i]/(1-max(speeds)))-m[i]
        other=tau*moment*(2*(R+24*I)+2160*per)
        roots.append(simple+min(H[i],multi+other))
    return max(dot(row,roots) for row in N)


def full_rectangle(ball,reg,prior_pure,hist,pairs,series,halos,reps):
    # Cotas por capas con una marca física necesaria para cubrir el síndrome.
    a,s=Q(".5488"),Q(".86")
    pure=[Counter(),Counter()]
    geometry=prop["cover_profiles"](reps)
    for i in (0,1):
        for (v,c),by_cap in hist[i].items():
            n=v+c
            if n==1:
                continue
            if n==2:
                j=c-i
                discount=effective["discount"](((12,118),(118,74))[i][j])
                pure[i][v,c]=discount*sum(mu*anchor(ball,i,k) for k,mu in enumerate(pairs[i][j]))
            else:
                price=sum(mu*anchor(ball,i,k) for k,mu in by_cap.items())
                ordinary=sum(mu*(a**6)**(-k) for k,mu in by_cap.items())
                if n<=prop["LIMIT"]:
                    saving=max(0,246*v+340*c-prop["group_budget"](v,c,geometry,n,guard=True))
                    discount=effective["discount"](saving)
                    price*=discount
                    ordinary*=discount
                reg[i][v,c]=min(reg[i][v,c],ordinary)
                pure[i][v,c]=min(price,prior_pure[i][v,c])
    w=dict(name="Rectángulo completo por paridad",rates=("3.28e-7",".0011475",".004",".000588","1.01"),
           m=("6.35112693404e-7","7.46425692065e-7"),u=("4.38872267915e-7","4.74900538897e-7"),
           h=("5.17875971060e-8","2.64899907414e-7"),physical=("2.28172494083e-9","1.44215071619e-8"),
           outer=(".00115472728643","3.24634619537e-7","7.18307491411e-6"))
    d=tuple(map(Q,("5.18215971492e-5",".000130454722605")))
    e=tuple(map(Q,("5.18215971492e-5",".000381581171565")))
    effective["certify"](a**6,s,halos,reg,pure,prop["kernels"](ball,a**6,s),
                         lambda x,i:effective["full_tail"](series,a,x,i,d,e),w)


def main():
    contact,stripe,base,history,earlier=map(previous,("contacto_de_checks","franjas","medidas_y_un_dato",
                                                     "medidas_por_historia","medidas_aisladas"))
    profiles,stable,reps=contact["geometry"]()
    winners=prop["winners"](reps)
    assert all(phi.bit_count()<=5*(2*h.bit_count()-phi.bit_count()) for h,phi in winners)
    pairs=effective["pair_classes"](winners,stable,effective["effective_budgets"](winners,stable))
    _,prior=stripe["cost_tables"](stripe["light_profiles"](),previous("fondo_pasivo")["step_costs"],
                                 previous("pares_adyacentes")["stripe_profiles"])
    coarse=contact["positions"](profiles,stable,prior,stripe["coefficients"](4))
    trees=stripe["coefficients"](10)
    tables=history["data_coefficients"](coarse)
    triples=earlier["count_triples"](profiles,stable,contact["rooted_trees"])
    series=earlier["group_series"](10)
    hist=cumulative(trees,tables,triples,series,pairs)
    refine_tails(hist,stripe["coefficients"](prop["LIMIT"]),reps)
    ball=history["ball"]
    assert [depth_from_cover(j) for j in range(6)]==[0,0,1,2,2,3]
    assert [distance_budget(l) for l in range(9)]==[0,1,2,3,4,5,6,8,9]
    m=Q("1.846e-6"),Q("2.141e-6")
    direct,hard=direct_two(ball,hist,m)
    response,first_group=single_group_response(ball,hist,m)
    pair=direct+response
    print(f"Dos componentes con alguna de tamaño >=3: directa {float(direct):.12g}, "
          f"interacción entre especies {float(response):.12g}, total {float(pair):.12g}.",flush=True)
    assert direct < Q('.000108826') and response < Q('.000012047') and pair < Q('.000121')
    Y=Q(43,2000)
    assert 1<BONUS<1/Y
    regular,marked=sector(Y,ball,hist,pairs)
    d=tuple(map(Q,WITNESS["d"]))
    old["certify"](WITNESS,Y,Q(".9999"),ball,contact["exposure"](reps,base["budgets"]),
                        regular,marked,lambda w,i:tail(trees,Y,w,i,d),
                        scale=1,group=Q(1),single_group_factor=1/Y,tree_factor=Q("1.006"),b=Q(".00115"))
    u,H=(tuple(map(Q,WITNESS[k])) for k in ("u","physical"))
    g=40*H[0]+27*H[1]+Q(WITNESS["rooted"])
    y=Q(".001")*upper(Q(".015"))/(1-Q(".002"))
    simple_factor,_=marked_sensitivity(Q(".00115"),y,u,H,g)
    weighted=[H[i]+(u[i]-H[i])*simple_factor[i] for i in (0,1)]
    broad=max(dot(row,weighted) for row in N)/(BONUS*Y)
    reg,_=prop["group_polynomials"](Q(2,3)**6,ball,tables,triples,series,pairs,
                                    prop["cover_profiles"](reps),stripe["coefficients"](prop["LIMIT"]))
    measure=density_projection(ball,reg,series)
    density=4**9*measure
    subtotal=broad+hard+response+density
    joint=first_group+subtotal
    assert broad < Q('.000911752') and hard+response < Q('.000111820') and density < Q('.000156022')
    assert subtotal < Q('.001180') and subtotal+Q('.0000234') < Q('.001203')
    assert first_group < Q('.001723672') and joint < Q('.002903264') and joint+Q('.0010967') < Q('.004')
    print(f"Dos o más con marca >=6Y: {float(broad):.12g}; dos fuera de las otras clases: "
          f"{float(hard+response):.12g}; capa de densidad: {float(density):.12g}.")
    print(f"Subtotal del resto v9 {float(subtotal):.12g}, todas las de un grupo {float(first_group):.12g}, "
          f"parte de la fila conjunta {float(joint):.12g}.")
    reg,pure0=prop["group_polynomials"](Q(".5488")**6,ball,tables,triples,series,pairs,
                                       prop["cover_profiles"](reps),stripe["coefficients"](prop["LIMIT"]))
    full_rectangle(ball,reg,pure0,hist,pairs,series,contact["exposure"](reps,base["budgets"]),reps)


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para verificar las desigualdades.")
    main()

