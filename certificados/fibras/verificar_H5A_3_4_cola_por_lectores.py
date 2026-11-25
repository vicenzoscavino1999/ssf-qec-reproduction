"""HGP (3,4). Árboles de grupos arbitrariamente grandes con lectores compartidos."""
from fractions import Fraction as Q
from functools import lru_cache
from math import ceil, comb
from pathlib import Path
from runpy import run_path


HERE = Path(__file__).resolve().parent
v21 = run_path(str(HERE/"verificar_H5A_3_4_salidas_del_mismo_grupo.py"))
v19 = v21["old"]
v18 = v21["prior"]
V15 = v18["previous"]("verificar_H5A_3_4_contacto_oculto.py")
F5 = v18["previous"]("verificar_FASE5_circuito_sindrome.py")
FAR = V15["v11"]["prop"]
BALL = V15["v11"]["previous"]("medidas_por_historia")["ball"]
ANCHOR = V15["v10"]["anchor"]
SEPARATED = ((2,4,8),(2,2,7,9))
SINGLE = ((2,4,12),(2,2,8,12))
SATELLITE = ((2,3,6),(2,2,6,6))
ORDINARY = ((2,3,6),(2,2,4,6))
FLOOR_FIELDS = tuple(map(Q,(".001","0","0","0","0","1e-8","3e-8","1e-7")))


def edge_kind(i,j,x,y):
    if i!=j:
        return 2 if x==y==1 else int(x==1 or y==1)
    if i:
        x,y=y,x
    return 2 if x==0 and y<=2 else int(x==0 or y<=2)


def sphere(i, r):
    return BALL(1-i,r)-(BALL(1-i,r-2) if r>=2 else 0)


@lru_cache(None)
def matrices(a,s,floor_fields=FLOOR_FIELDS):
    matrices = [[[Q(0),Q(0)] for i in (0,1)] for kind in range(3)]
    counts = [[[0,0] for i in (0,1)] for kind in range(3)]
    for i in (0,1):
        for j in (0,1):
            left,right = sorted((i,j))
            for x,y in v21["displacements"](left,right):
                cost = sum(n*w for n,w in zip(v21["initial_overlap"](left,right,x,y),floor_fields[:5]))
                cost += sum(v21["link_overlap"](left,right,x,y,t,a**6,s)*h
                            for t,h in enumerate(floor_fields[5:]))
                count=sphere(i,x)*sphere(i,y)
                kind=edge_kind(i,j,x,y)
                counts[kind][i][j]+=count
                matrices[kind][i][j]+=count*v21["upper_factor"](cost)
    for m,name in zip(counts,("WEAK","STRONG_DISJOINT","OVERLAP")):
        assert tuple(map(tuple,m))==v19[name]
    return tuple(tuple(map(tuple,m)) for m in matrices)


def exact_disjoint(i, law, rate, first, theta, a, fees):
    n=3+i
    wait,ctrl=2*rate/(3-2*rate),8*rate/(15-8*rate)
    masses=[first*wait+ctrl/2]+[ctrl]*(n-1)+[(12-n-first)*wait+ctrl/2]
    def weight(r):
        return a**(sum(fees[:r])-6*(2*r-n)) if 2*r>=n else Q(0)
    result=law[1]*max(Q(1),weight(n))
    for omitted, mass in enumerate(masses):
        for x in range(n-omitted+1):
            for y in range(omitted+1):
                prob=comb(n-omitted,x)*comb(omitted,y)*theta**(x+y)*(1-theta)**(n-x-y)
                result+=mass/law[0]*prob*weight(n-omitted-x+y)
    return result


def source_prices(a, rate=None, first=None):
    result=[]
    for i in (0,1):
        if rate is None:
            single=Q(1)
            separate=a**-4
            satellite=a**(-(7+i))
            plain=overlap=a**(-(7+3*i))
        else:
            gates=first if i==0 else 7
            theta=F5["THETA"]
            law=v19["circuit_law"](i,rate,gates,theta)
            single,separate,satellite,plain=(exact_disjoint(i,law,rate,gates,theta,a,fees[i])
                                      for fees in (SINGLE,SEPARATED,SATELLITE,ORDINARY))
            overlap=sum(p*a**(-k) for k,p in law[3].items())
        assert 0<single<=separate<=satellite<=plain<=overlap
        result.append((single,separate,(satellite-single)/2,plain,overlap))
    return tuple(result)


def pair_prices(i,j,energies,Y,discount,old_regular,old_mark,circuit):
    """Precio por desplazamiento de dos lectores, con las clases de saldo de v12/v18."""
    def count(xy):
        return sphere(i,xy[0])*sphere(i,xy[1])
    if circuit:
        allowed=set(energies)
        aligned=set()
        if i==j==0:
            allowed={p for p in energies if p[0]==0}
            aligned={(0,2)}
        elif i==j==1:
            allowed=allowed-{(4,2)}
            aligned={(2,0),(4,0)}
        else:
            allowed.discard((1,1))
        all_cost=max(energies.values())
        one_cost=max(energies[p] for p in allowed)
        # La fila de parejas comparte los escalones de k con la marca de generador.
        second=old_mark*(Y**-2-Y**-1)/ANCHOR(BALL,i,2) if aligned else 0
        assert 0<=second<=old_regular
        zero=min(old_regular-second,discount*sum(count(p) for p in energies))
        return (old_regular-second)*one_cost+zero*(all_cost-one_cost)+second*max(
            (energies[p] for p in aligned),default=0)
    if i==j==0:
        clean=Y**2+Q(".152")
        return discount*(sum(count(p)*energies[p]*clean for p in energies if p[0]!=0)
                         +54*energies[0,4]/Y+9*energies[0,2]/Y**2)
    if i==j==1:
        clean0=Y**2+Q(".528")
        clean1=Y**2+Q(".528")/Y
        def top(k):
            total=0
            for pt in sorted(((2,0),(4,0)),key=energies.get,reverse=True):
                charge=min(count(pt),k)
                total+=charge*energies[pt]
                k-=charge
            return total
        return discount*(384*clean0*energies[4,2]+sum(count(p)*clean1*energies[p]
                         for p in ((0,2),(2,2)))+top(56)+(Y**-1-1)*top(32)
                         +(Y**-2-Y**-1)*top(12))
    clean0,clean1=Y**2+Q(".167"),Y**2+Q(".167")/Y
    weak=(3,3)
    return discount*(sum(count(p)*energies[p]/Y for p in energies if p!=weak)
                     +24*energies[weak]/Y+count((1,1))*energies[1,1]
                     +count(weak)*energies[weak]*(clean0+clean1/2))


def with_multiple_limits(matrix, limits):
    """La cuarta matriz cuenta destinos débiles promovidos a lectura múltiple."""
    if limits is None:
        return matrix
    weak = v19["WEAK"]
    assert all(0<=2*limits[i][j]<=weak[i][j] for i in (0,1) for j in (0,1))
    multiple = tuple(tuple(matrix[0][i][j]*Q(limits[i][j],weak[i][j]) for j in (0,1)) for i in (0,1))
    return (*matrix,multiple)


def weak_matrices(matrix):
    weak = matrix[0]
    multiple = matrix[3] if len(matrix)==4 else tuple(tuple(x/2 for x in row) for row in weak)
    return weak,multiple


def transfer(t, coefficients, matrix, state, exponential, exp_minus_one):
    """U,D,O,C,T y el incremento por un enlace débil de lectura múltiple."""
    weak,multiple=weak_matrices(matrix)
    strong,overlap=matrix[1:3]
    free, inc_strong, inc_overlap, core, tagged, premium=state
    result=[[] for _ in state]
    for i,(single,base,extra,plain,high) in enumerate(coefficients):
        def dot(row,z):
            return row[0]*z[0]+row[1]*z[1]
        potential=dot(multiple[i],free)
        good=dot(weak[i],free)-potential
        bad=potential+dot(multiple[i],premium)
        ordinary=exponential(good)
        forest=exponential(good+bad)
        attached=exponential(dot(weak[i],tagged))
        ds,os=dot(strong[i],inc_strong), dot(overlap[i],inc_overlap)
        dc=t[i]*forest*attached*(plain*exp_minus_one(ds)+high*exponential(ds)*exp_minus_one(os))
        dt=t[i]*extra*forest
        result[0].append(t[i]*ordinary*(single+base*exp_minus_one(bad))+dc+dt*dot(weak[i],core))
        result[1].append(t[i]*plain*forest*attached+dc)
        result[2].append(t[i]*high*forest*attached*exponential(ds+os))
        result[3].append(dc)
        result[4].append(dt)
        result[5].append(t[i]*(base-single)*ordinary)
    return result


def total(t, prices, matrices):
    state=[[0,0] for _ in range(6)]
    exponential=F5["exp_upper"]
    for _ in range(500):
        new=transfer(t,prices,matrices,state,exponential,lambda x:exponential(x)-1)
        if max(max(row) for row in new)>Q(".02"):
            raise ValueError("La iteración rebasa el dominio de esta supersolución.")
        new=[[Q(ceil(x*10**24),10**24) for x in row] for row in new]
        assert all(x<=y for old,after in zip(state,new) for x,y in zip(old,after))
        if new==state:
            return state[0]
        state=new
    raise ValueError("La iteración no se estabilizó.")


def polynomials(prices,matrix,nmax):
    """Coeficientes bivariados de las mismas series, sin restas entre mayorantes."""
    free,inc_strong,inc_overlap,core,tagged,premium=({},{},{},{},{},{})
    logs=[[{} for _ in range(5)] for i in (0,1)]
    forests=[[{(0,0):1} for _ in range(5)] for i in (0,1)]
    weak,multiple=weak_matrices(matrix)
    strong,overlap=matrix[1:3]
    # Una sola escala por tipo de raíz. Las claves completas son (i,v,c).
    def value(poly,i,v,c):
        return poly.get((i,v,c),0)
    def linear(edge,poly,i,v,c):
        return sum(edge[i][j]*value(poly,j,v,c) for j in (0,1))
    for n in range(1,nmax+1):
        for v in range(n+1):
            c=n-v
            for i,(single,base,extra,plain,high) in enumerate(prices):
                u,b=v-(1-i),c-i
                if u<0 or b<0:
                    continue
                fg,f0,f1,f2,f3=[f[u,b] for f in forests[i]]
                key=(i,v,c)
                co=plain*(f2-f1)+high*(f3-f2)
                cross=sum(forests[i][1][u-x,b-y]*linear(weak,core,i,x,y)
                          for x in range(u+1) for y in range(b+1) if x+y)
                core[key]=co
                tagged[key]=extra*f0
                free[key]=single*fg+base*(f0-fg)+co+extra*cross
                inc_strong[key]=plain*f1+co
                inc_overlap[key]=high*f3
                premium[key]=(base-single)*fg
        for i in (0,1):
            for v in range(n+1):
                c=n-v
                all_weak=linear(weak,free,i,v,c)
                cg=all_weak-linear(multiple,free,i,v,c)
                co=all_weak+linear(multiple,premium,i,v,c)
                cc=co+linear(weak,tagged,i,v,c)
                cd=cc+linear(strong,inc_strong,i,v,c)
                ce=cd+linear(overlap,inc_overlap,i,v,c)
                for z,f,coef in zip(logs[i],forests[i],(cg,co,cc,cd,ce)):
                    z[v,c]=coef
                    agg=sum((x+y)*z[x,y]*f[v-x,c-y] for x in range(v+1) for y in range(c+1) if x+y)
                    f[v,c]=agg*Q(1,n)
    return tuple({(v,c):val for (j,v,c),val in free.items() if j==i} for i in (0,1))


CASES = (
    ("Bernoulli", "5.22e-7", None, "0.567", "0.685",
     (("9.843567200397559e-7", "0.0000011479151954169456"),
      ("6.907998652116276e-7", "7.409401423513036e-7"),
      ("3.938629566625392e-8", "1.3781647732343652e-7"),
      ("4.776774520600804e-9", "1.9253236511750095e-8"),
      "5.356289080880986e-7", "0.000001804633701797979")),
    ("VV primero", "9.75e-8", 4, "0.558", "0.75",
     (("0.0000013944453202328718", "0.000001597683814451405"),
      ("4.500377813455979e-7", "8.43247822776818e-7"),
      ("5.220621315603707e-8", "2.6528080034774725e-7"),
      ("6.555654762326946e-9", "3.587581818878111e-8"),
      "7.057503276460359e-7", "0.000001406004314147712")),
    ("Cualquier orden de columnas", "8.9e-8", 7, "0.555", "0.78",
     (("0.0000012760673279608317", "0.0000014630123547883024"),
      ("6.488191875129268e-7", "7.820761227784187e-7"),
      ("8.827848726141487e-8", "3.075013258652207e-7"),
      ("9.681416799569123e-9", "4.3890538432470264e-8"),
      "9.115077146388823e-7", "0.0000014182924222685875")))


def rounded_prices(prices):
    def up(x):
        return Q(ceil(x*10**13),10**13)
    # Se redondean hacia arriba las actividades primitivas: también el recargo.
    return tuple((up(single), up(single)+up(base-single), up(extra), up(plain), up(high))
                 for single,base,extra,plain,high in prices)


def ordinary_heads(data,circuit,rate,first,a):
    theta=F5["THETA"]
    if first is None:
        halo,_,discount,hist=data
        reg,mark=v18["physical_prices"](V15,a,theta,BALL,discount,hist)
        return halo,discount,reg,mark,None
    halo,trees,crowded,discount,geometry,trip,hc,hr=circuit
    law=[v19["circuit_law"](i,rate,first if i==0 else 7,theta) for i in (0,1)]
    reg,mark=F5["polys"](hr,trees,crowded,discount,geometry,trip,law,FAR["N"],ANCHOR,BALL,a**6)
    extra,extra_mark=v19["surcharge"](hc,hr,law,discount,a**6,ANCHOR,BALL)
    for i in (0,1):
        reg[i].update(extra[i])
        mark[i].update(extra_mark[i])
        key=(1-i,i)
        cost=(v19["singleton_cost"](0,rate,first,theta,a**6) if i==0
              else v21["singleton_cc"](rate,7,theta,a**6))
        reg[i][key]=min(reg[i][key],discount[key]*cost)
        assert max(hc[i][2-2*i,2*i])<=2 and max(hr[i][2-2*i,2*i])<=2
    return halo,discount,reg,mark,law


def field_coefficients():
    theta,b=F5["THETA"],F5["B_PASS"]
    def logbound(x):
        return x+x*x/(2*(1-x))
    collision=Q(89,2)*logbound(theta*theta)
    guard=logbound(theta)-collision
    return guard,b+collision-guard


def partial(poly,weight,small,large):
    return sum(coef*weight[0]**v*weight[1]**c for (v,c),coef in poly.items() if small<=v+c<=large)


def verify_model(case,data,circuit,*,floor_fields=FLOOR_FIELDS,disjoint_pair_factors=None,
                 multiple_limits=None,adjust_heads=None):
    """Las cotas opcionales por par se justifican para dos franjas físicamente disjuntas."""
    name,rate,first,a,s,witness=case
    rate,a,s=map(Q,(rate,a,s))
    m,u,h,H=[tuple(map(Q,row)) for row in witness[:4]]
    visible,root=map(Q,witness[4:])
    b,theta,eta,tilt,measurement=map(F5.get,("B_PASS","THETA","ETA","TILT","RM"))
    E,dot,evaluate=F5["exp_upper"],F5["dot"],F5["evaluate"]
    y0=a**6
    halo,discount,reg,mark,law=ordinary_heads(data,circuit,rate,first,a)
    if adjust_heads is not None:
        adjust_heads(reg,mark,discount,data,circuit,rate,first,a,s,floor_fields)
    raw_reg=[x.copy() for x in reg]
    raw_mark=[x.copy() for x in mark]
    alpha=Q(".004") if first is None else F5["ALPHA"]
    if law is None:
        y=measurement*E(tilt)/(1-2*measurement)
        z=[rate*E(tilt+alpha)/(1-2*rate)]*2
        adm=(1,1)
    else:
        target=4*rate/(15-8*rate)
        y=((1-2*rate/3)**-5/(1-measurement)/(1-target)**7-1)*E(tilt)
        z=[v[0]*E(tilt+alpha) for v in law]
        adm=[sum(v[4].values()) for v in law]
    assert y/(1+y)<=theta
    guard,other=field_coefficients()
    g=40*H[0]+27*H[1]+root
    children=(*h,visible)
    fields=(guard,other,2*u[0],2*u[1],2*g,*children)
    assert len(floor_fields)==len(fields) and min(floor_fields)>=0
    assert all(floor<=value for floor,value in zip(floor_fields,fields))
    dd,dg,gd,gg=FAR["kernels"](BALL,y0,s)
    v21["pair_mark"](mark,ANCHOR,BALL,discount,law is not None)
    v21["joint_rows"](reg,mark,b,theta,y0,s,u,h,H,visible,root)
    prices=rounded_prices(source_prices(a,rate if law else None,first))
    matrix=with_multiple_limits(matrices(a,s,floor_fields),multiple_limits)
    poly=polynomials(prices,matrix,13)
    for i in (0,1):
        forbidden=dot(v21["DIRECT"][i],children)
        assert all(value>=0 for value in poly[i].values())
        for key in reg[i]:
            reg[i][key]=min(reg[i][key],poly[i][key]*discount.get(key,1))
        for j in (0,1):
            left,right=sorted((i,j))
            factors={}
            for x,ypos in v21["displacements"](left,right):
                loss=dot(v21["initial_overlap"](left,right,x,ypos),fields[:5])+forbidden
                loss+=sum(u0*v21["link_overlap"](left,right,x,ypos,t,y0,s) for t,u0 in enumerate(children))
                factors[x,ypos]=v21["upper_factor"](loss)
                if disjoint_pair_factors is not None and edge_kind(left,right,x,ypos)!=2:
                    factor=disjoint_pair_factors[left,right]
                    assert 0<factor<=1
                    factors[x,ypos]*=factor
            key=(2-i-j,i+j)
            best=pair_prices(i,j,factors,y0,discount[key],raw_reg[i][key],
                             raw_mark[i].get(key,0),law is not None)
            reg[i][key]=min(reg[i][key],best)
    R=max(dot(row,u) for row in FAR["N"])
    I=280*u[0]+189*u[1]+588*g
    phi=E(89*b+I)
    eps=phi*(phi-1)
    assert eps<1 and R<alpha and max(u)<eta
    ratios=[y*(1+b)*(1+eps)/b,y*(1-theta)*E(99*guard+90*other+2*I)/theta]
    w,depart=[],[]
    for i in (0,1):
        exposure=2*dot([ar[i] for ar in halo[1:]],(u[0],u[1],g))
        rhs=z[i]*E(eta+2*R+guard*(504,628)[i]+other*halo[0][i]+exposure)
        ratios.append(rhs/m[i])
        w.append(F5["rounded"](m[i]*E(dot([v[0] for v in dd[i]],h)+dg[i][0]*visible)))
        depart.append(dot([v[1] for v in dd[i]],h)+dg[i][1]*visible)
    gross=total(w,prices,matrix)
    tail=[gross[i]-partial(poly[i],w,1,10) for i in (0,1)]
    naked_tail=[partial(poly[i],m,11,13) for i in (0,1)]
    assert all(t>=n>=0 for t,n in zip(tail,naked_tail))
    tree_factor=Q("1.01") if law is None else Q("1.014")
    trees=[mi*tree_factor for mi in m]
    init=((210,252),(640,784),(441,513),(1344,1596))
    removed=[]
    for i in (0,1):
        group,naked=evaluate(reg[i],w),evaluate(reg[i],m)
        removed.append(naked+naked_tail[i])
        pointed=sum(val*w[0]**v*w[1]**c*(v*depart[0]+c*depart[1]) for (v,c),val in reg[i].items())
        connected=E(dot(FAR["N"][i],trees))
        assert m[i]*connected<=trees[i]
        fixed=min(group-naked,pointed)+tail[i]-naked_tail[i]+m[i]*(connected-1)
        simple=z[i]*E(eta+2*R+dot((b,2*u[0],2*u[1],2*g),[r[i] for r in init]))*(1+y)**(3+i)*adm[i]
        ratios.extend(((group+tail[i])/h[i],fixed/H[i],(simple+H[i])/u[i]))
    J=dot([k[0] for k in gd],h)+gg[0]*visible
    one=dot([k[1] for k in gd],h)+gg[1]*visible
    raw=E(J)*one-dot([k[1] for k in gd],removed)
    assert raw>0
    outer=y0/s**2
    distant=outer*(outer*raw+(E(J)-1)*one/2)+V15["single_tail"](m)+sum(evaluate(tab,m) for tab in mark)
    ratios.extend((outer*(E(J)-1-J)/visible,distant/root))
    assert all(0<r<1 for r in ratios), (name,[(j,float(r)) for j,r in enumerate(ratios) if not 0<r<1])
    print(f"{name}, tasa={float(rate):g}, medición=.001. Fila {float(R):.12g}, máximo {float(max(ratios)):.12f}.")
    print(f"  Pasivo {float(ratios[0]):.12f}, activo {float(ratios[1]):.12f}, "
          f"colas {float(tail[0]):.9g}, {float(tail[1]):.9g}.",flush=True)


def check_local():
    contacts=V15["v11"]["previous"]("contacto_de_checks")
    profiles,stable,reps=contacts["geometry"]()
    assert len(FAR["winners"](reps))==771
    ext=v18["previous"]("verificar_H5A_3_4_grupos_extendidos.py")
    assert (ext["GENERAL"],ext["SEPARATED"],ext["SINGLE_CELL"])==(ORDINARY,SEPARATED,SINGLE)
    ext["check_credits"](profiles,stable)
    V15["v11"]["previous"]("medidas_por_historia")["far_credits"](profiles,stable)


def lower_exp(x):
    assert x>=0
    term=result=Q(1)
    for k in range(1,11):
        term*=x/k
        result+=term
    return result


def obstruction(data,a,s,*,floor_fields=FLOOR_FIELDS,multiple_limits=None,
                compositions=((1,3),(2,2)),adjust_heads=None):
    """Minorante de la fila ordinaria CC en G3, para las composiciones indicadas."""
    assert compositions and len(set(compositions))==len(compositions)
    assert all(c>=1 and v>=0 and 3<=v+c<=10 for v,c in compositions)
    halo,discount,reg,mark,_=ordinary_heads(data,None,Q("1e-6"),None,a)
    if adjust_heads is not None:
        adjust_heads(reg,mark,discount,data,None,Q("1e-6"),None,a,s,floor_fields)
    dd,dg,_,_=FAR["kernels"](BALL,a**6,s)
    matrix=with_multiple_limits(matrices(a,s,floor_fields),multiple_limits)
    poly=polynomials(rounded_prices(source_prices(a)),matrix,max(v+c for v,c in compositions))
    guard,other=field_coefficients()
    p=Q("1e-6")
    subtotal=Q(0)
    for v,c in compositions:
        # Centro CC; los v VV y los c-1 CC se añaden en una estrella. Se evalúa
        # una opción de la minimización usada por la fila, no una realización física.
        edges=((v,0,1,3,3),(c-1,1,1,4,2))
        cut=[sum(n*v21["initial_overlap"](i,j,x,y)[k] for n,i,j,x,y in edges) for k in range(5)]
        save=[v21["DIRECT"][1][t]+sum(n*v21["link_overlap"](i,j,x,y,t,a**6,s)
               for n,i,j,x,y in edges) for t in range(3)]
        assert all(cap<=v*r[0]+c*r[1] for cap,r in zip(cut[2:],halo[1:]))
        assert save[2]<=v*dg[0][0]+c*dg[1][0] and save[0]<=v*dd[0][0][0]+c*dd[1][0][0]
        scale=(v+c)*Q(".0191")+guard*(v*504+c*628-cut[0])+other*(v*halo[0][0]+c*halo[0][1]-cut[1])
        coef=p**(v+c)*lower_exp(scale)*min(reg[1][v,c],poly[1][v,c]*discount.get((v,c),1))
        slope=v*dd[0][1][0]+c*dd[1][1][0]-ceil(save[1])
        assert slope>0 and coef>0
        subtotal+=coef*slope
        print(f"  {v} VV + {c} CC: c={float(coef):.10g}, L={slope}, c L={float(coef*slope):.9f}.")
    # exp(Lh) >= e Lh para L,h>0. Se necesita e Σ cL <= 1.
    ratio=lower_exp(Q(1))*subtotal
    print(f"G3, a={float(a):g}, s={float(s):g}: minorante de e Σ cL = {float(ratio):.10f}.",flush=True)
    return ratio


def check_g3(data):
    a,s=Q(".569"),Q(".673")
    mass=(Q("1.84e-6"),Q("2.14e-6"))
    old=v21["bulk_degrees"](mass,a,FAR["W"],FAR["B"],11)
    series=polynomials(rounded_prices(source_prices(a)),matrices(a,s),11)
    for i in (0,1):
        after=partial(series[i],mass,11,11)
        before=old[i][11]
        assert 0<after<before
        print(f"G3, raíz {('VV','CC')[i]}, grado 11: {float(before):.10g} -> {float(after):.10g},"
              f" factor {float(before/after):.7g}.",flush=True)
    assert obstruction(data,a,s)>Q("1.0239")
    assert obstruction(data,Q(".567"),Q(".685"))>Q("1.0066")


def main():
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para comprobar las desigualdades.")
    check_local()
    data,circuit=v21["previous_tables"]()
    for case in CASES:
        verify_model(case,data,circuit)
    check_g3(data)


if __name__=="__main__":
    main()

