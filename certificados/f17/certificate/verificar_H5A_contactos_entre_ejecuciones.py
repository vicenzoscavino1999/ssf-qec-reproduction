"""Certificado de H5A_contactos_entre_ejecuciones.txt. Python 3.10 o posterior."""

from dataclasses import dataclass
from fractions import Fraction as F
from functools import cache
from math import ceil, comb
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verificar_H5A_resto_global as old
import verificar_H5A_propagacion_entre_fuentes as earlier

exp = old.prior.exp_arriba
dot = old.prior.escalar
odds = old.prior.actividad
Q = F('2e-11')
TARGET = F('1e-30')
SCALE = 10**6
WIDTHS = (17,18,1,17,18,306)
READERS = (18,17,1,18,17,306)
TOP = F('.04')
DIRECT = F('.05')
XI = F('2.055e-11')
Z_CAP = tuple(map(F, ('.075','.075','.001')))
LEVELS = sorted({F(2*k-n,n) for n in (17,18) for k in range(n//2+1,n+1)}, reverse=True)


def fractions(values):
    return tuple(map(F,values.split()))


@dataclass(frozen=True)
class Certificate:
    rate: F
    direction: tuple
    cutoff: int
    dual: tuple
    rebates: tuple
    sources: tuple
    branches: tuple


GENERAL = Certificate(F('1.2e-34'),fractions('.0018381 .0010339 .000000721'),1,
    (997066,0,30848,445312,65690,475000,109003,510045,165007,548077,239483,
     595052,337121,649148,471023,714062,628993,131944),
    ((142187,)*8+(534722,)+(562500,)*7+(983455,), (142187,)*9+(562500,)*8+(983455,)),
    fractions('3.27933854339e-24 1.5719183199e-23 7.90422623235e-28'),
    fractions('7.3938427233e-19 7.39384610412e-19 1.0398200053e-17 '
              '1.09060747771e-11 1.08785965592e-11 1.63006058771e-10 '
              '5.74755036827e-24 6.08563277183e-24 1.45480509221e-21 '
              '2.81200687712e-22 2.97742040766e-22 5.06228122804e-21 '
              '3.94080317744e-7 4.17260894629e-7 7.09370815268e-6'))

PARTIAL = Certificate(F('2.1e-30'),fractions('.0024 .0005 .0000006'),3,
    (1009199,0,23208,349609,50651,372917,84524,400670,128348,430288,186899,
     467448,263376,509943,368608,560938,343750,0),
    ((110937,)*8+(494791,)+(665625,)*7+(995588,), (110937,)*9+(665625,)*8+(995588,)),
    fractions('1.00777146e-24 3.95772622e-24 2.99269608e-28'),
    fractions('1.7667773533e-19 1.73213506619e-19 3.18041472431e-18 '
              '4.57386135764e-12 4.47561044683e-12 8.52226491963e-11 '
              '7.04166865898e-25 7.30968664037e-25 2.28163371184e-22 '
              '7.36248531587e-23 7.79557451469e-23 1.32533713488e-21 '
              '2.02034828976e-7 2.13919037907e-7 3.63668715222e-6'))


@cache
def bits(p):
    """Coste superior en bits, redondeado a 1/64, con prueba entera."""
    if p >= 1:
        return F(0)
    assert p > 0
    num,den = p.numerator**64,p.denominator**64
    k = max(0,den.bit_length()-num.bit_length())
    if num << k < den:
        k += 1
    assert num << k >= den and (not k or num << (k-1) < den)
    return F(k,64)


def up(x):
    return F(ceil(x*10**75),10**75)


@cache
def two_fraction(k):
    base = F('1.00271128')
    assert base**256 > 2 and 0 <= k < 256
    return up(base**k)


def two_upper(microbits):
    whole,k = divmod(ceil(F(256*microbits,SCALE)),256)
    assert whole >= 0
    return 2**whole * two_fraction(k)


@cache
def shapes():
    """Mínimos de luces alineadas por celda emitida. Las cotas admiten ciclos."""
    sig = {}
    for x in range(19):
        for y in range(18):
            a = 17*x+18*y-2*x*y
            if not x+y or 2*a < 17*(x+y):
                continue
            for h in range(a//2+1,a+1):
                d = 2*h-a
                if 100*d < a:
                    continue
                hf,hv = min(17-y,17*h//a),min(18-x,18*h//a)
                if x*hf+y*hv < h:
                    continue
                if h == a:
                    sig.setdefault((a,d),set())
                    continue
                arms = ((True,x*(17-y),17-y,x),(False,(18-x)*y,y,18-x))
                for row_arm,count,rf,rv in arms:
                    if not count:
                        continue
                    f = max(0,ceil(F((2*rf-17)*h,a)))
                    v = max(0,ceil(F((2*rv-18)*h,a)))
                    if row_arm:
                        f = max(f,h-(x-1)*hf-y*hv)
                    else:
                        v = max(v,h-x*hf-(y-1)*hv)
                    if f < rf and v < rv and f+v <= h:
                        sig.setdefault((a,d),set()).add((f,v))
    return sig


def check_local(cert):
    lam = cert.dual
    assert len(lam)==1+len(LEVELS) and min(lam)>=0
    for n,r in zip((17,18),cert.rebates,strict=True):
        assert len(r)==n and tuple(sorted(r))==r and 0 <= r[0] <= r[-1] <= lam[0]
    kr,kc,kd = cert.direction
    assert 0 < kd <= min(kr,kc)

    @cache
    def ordinary(h,f,v):
        counts = (f,v,h-f-v)
        limits = [TOP]
        for subset in range(8):
            num = sum(n for i,n in enumerate(counts) if subset>>i&1)
            if num < 4:
                available = sum(n*w for i,(n,w) in enumerate(zip(counts,cert.direction))
                                if not (subset>>i&1))
                limits.append(4*available/(4-num))
        return min(limits)

    charges = {1:F(0),2:F(0)}
    for (a,d),opts in shapes().items():
        h,c = (a+d)//2,(a-d)//2
        normal = min((ordinary(h,f,v) for f,v in opts),default=TOP)
        assert normal > 0 and (c==0 or opts)
        mu = lam[0]*d+sum(w*(d if d*u.denominator < a*u.numerator else -c)
                              for w,u in zip(lam[1:],LEVELS,strict=True))
        rebate = max(r[n*h//a-1] for n,r in zip((17,18),cert.rebates,strict=True))
        assert mu >= 0
        if d >= cert.cutoff or a <= 18:
            assert mu >= SCALE*bits(normal), (a,d,'sin luz inicial')
        else:
            assert mu >= SCALE*bits(normal+DIRECT)+rebate, (a,d,'una luz inicial')
            charges[d] = max(charges[d],bits(normal)-F(mu,SCALE))
        assert mu >= SCALE*bits(normal+h*DIRECT)+h*rebate, (a,d,'extremo superior')
    return charges


def group_profiles(cert):
    lam,ranks = cert.dual,cert.rebates
    fields = [(ceil(17*(1+u)/2)-1,ceil(18*(1+u)/2)-1,ceil(17*(1+u)/2),w)
              for u,w in zip(LEVELS,lam[1:],strict=True) if w]
    linear = tuple(n*lam[0]+sum(row[i]*row[3] for row in fields)-sum(ranks[i])
                   for i,n in enumerate((17,18)))

    def profile(v,c):
        @cache
        def tail(index,k,balance):
            if index == len(fields):
                return 0
            bv,bc,req,weight = fields[index]
            best = 0
            for j in range(min(k,bc if c else bv)+1):
                used = max(req,(17*k+17)//18) if j < k else 0
                remain = min(v*min(j,bv)+c*min(j,bc),max(0,balance-used))
                best = max(best,weight*remain+tail(index+1,j,remain))
            return best
        values = []
        for k in range(19 if c else 18):
            initial = v*min(17,k)+c*k
            rebate = v*sum(ranks[0][:k])+c*sum(ranks[1][:k])
            values.append(lam[0]*initial-rebate+tail(0,k,initial))
        assert max(values) <= linear[0]*v+linear[1]*c
        return values
    return linear,{e:[profile(v,e-v) for v in range(e+1)] for e in range(1,17)}


@cache
def overlap_rows():
    # Coeficientes por número de datos de un grupo con k inicial <=16.
    x,y = old.X,old.XI
    active,measured,site,tree = [[F(0)]*17 for _ in range(4)]
    measured[0] = y

    def exponential(series,degree):
        values = [up(exp(series[0]))]
        for j in range(1,degree+1):
            values.append(up(sum(k*series[k]*values[j-k] for k in range(1,j+1))/j))
        return values

    for e in range(1,17):
        neighbours = [46870021*v+1420146*d for v,d in zip(active,measured,strict=True)]
        site[e] = up(x*exponential(neighbours,e-1)[e-1])
        internal = exponential([612*t for t in tree],e-1)
        extra = sum(site[e-j]*internal[j] for j in range(1,e))
        tree[e] = up(site[e]+extra)
        active[e] = up(18*y*site[e]+extra)
        measured[e] = up(y*exponential([2766309*v for v in active],e)[e])
    return active,active,[u+35*y*b for u,b in zip(measured,site,strict=True)]


def check_groups(cert,n,i_price):
    linear,profile = group_profiles(cert)
    x = [i_price*exp(z) for z in Z_CAP[:2]]
    y = XI*exp(Z_CAP[2])
    ref = (old.X*F('.16452'),old.X)
    assert y < old.XI
    ratio = max(u*two_upper(b)/r for u,b,r in zip(x,linear,ref,strict=True))
    assert ratio < 1
    data_roots = []
    measure_root = F(0)
    for j,width in enumerate((17,18)):
        by_k = profile[1][1-j]
        factor = x[j]*exp(n[j][2]*y)
        data_roots.append(factor*sum(two_upper(by_k[width-m])*comb(width,m)*y**m
                                     for m in range(width+1)))
        measure_root += 35*factor*sum(two_upper(by_k[width-m])*comb(width-1,m-1)*y**m
                                      for m in range(1,width+1))
    measure_root += 2766309*y*sum(data_roots)
    sums = [*data_roots,measure_root]
    low_rows = overlap_rows()
    for e in range(2,17):
        direct,high,low = [],[],[]
        for v,values in enumerate(profile[e]):
            multiplier = (x[0]/ref[0])**v*(x[1]/ref[1])**(e-v)
            direct.append(multiplier*two_upper(max(values)))
            high.append(multiplier*two_upper(max(values[17:])))
            low.append((x[0]/old.X)**v*(x[1]/old.X)**(e-v)*two_upper(max(values[:17])))
        for i in range(3):
            sums[i] += min(old.U[i][e]*max(direct),
                           old.U[i][e]*max(high)+low_rows[i][e]*max(low))
    for i,limit in enumerate(cert.sources):
        sums[i] += old.U[i][17]*ratio**17
        assert sums[i] < limit, (i,'fuentes',float(sums[i]),float(limit))


def matrices():
    base,_,_ = old.prior.geometria()
    order = (0,1,3,0,1,2)
    m = [[base[i][j] for j in order] for i in order]
    n = [row[:3] for row in m[:3]]
    n[2][2] = 0
    return m,n


def check_trees(cert,m):
    kr,kc,kd = cert.direction
    lam = 289*kr+289*kc+83232*kd
    assert 0 < lam < 1
    h = [[f*l for l in READERS] for f in WIDTHS]
    path = [[entry/(1-lam) for entry in row] for row in h]
    s,p,c,near,near_child = [cert.branches[j:j+3] for j in range(0,15,3)]
    both = [a+b for a,b in zip(s,near,strict=True)]
    w = [dot(row,(*cert.sources,*both)) for row in m]
    children = [dot(row[3:],p) for row in path]
    z = [b+l for b,l in zip(w,children,strict=True)]
    for i in (0,1):
        z[i] += dot(h[i][3:],near_child)
    assert all(a<b for a,b in zip(z[:3],Z_CAP,strict=True))
    rates = (max(kr,TOP/17),max(kc,TOP/18),max(kr,kc,kd,TOP/33))
    for j in range(3):
        i = j+3
        before = dot(path[i],(*cert.sources,*(a+b for a,b in zip(c,near,strict=True))))
        straight = dot(h[i][:2],cert.sources[:2])
        attached = WIDTHS[i]*XI
        many = TOP*attached**3*exp(z[i]+attached)
        g = exp(z[i]+2*attached)
        assert s[j] > many+rates[j]*before*g
        assert p[j] > rates[j]*(g-1-children[i])
        assert c[j] > many+rates[j]*before*(g-1)
        assert near[j] > DIRECT*straight*exp(z[i])
        assert near_child[j] > DIRECT*(exp(z[i])-1)
    core = max(cert.sources[0]+both[0]+17*both[2],cert.sources[1]+both[1]+18*both[2])
    return core,cert.sources[2]+35*core


def check_exterior(cert,core,measured,data_price):
    lv,lm = F('1.045')*cert.rate,F('2.055e-11')
    vv,vm,mv,mm = 46870021,1420146,2766309,83810
    tilt = F('.025')
    assert lv > odds(cert.rate)*exp(tilt+vv*(lv+core)+vm*(lm+measured))
    assert lm > odds(Q)*exp(tilt+mv*(lv+core)+mm*(lm+measured))
    assert XI > odds(Q)*exp(tilt+mv*lv+mm*lm)
    assert data_price > odds(cert.rate)*exp(tilt+135*(vv*lv+vm*lm))
    rv,rm = core+18*lv,measured+lm
    assert max(135*(vv*rv+vm*rm),mv*rv+mm*rm) < F('.01')
    assert rv < F('2e-16') and rm < F('2.057e-11')
    return rv,rm


def first_use_rows(m,n):
    """Criterio agregado por primer uso, incluso si hay recaídas posteriores."""
    direction = fractions('.00265 .000225 .0000008')
    areas = {17*x+18*y-2*x*y for x in range(19) for y in range(18)
             if x+y >= 2 and 2*(17*x+18*y-2*x*y) >= 17*(x+y)}
    types = [(17,direction[0]),(18,direction[1])]+[(a,direction[2]) for a in areas]
    actions = [(a,h,2*h-a,bits(h*price)) for a,price in types
               for h in range(a//2+1,a+1) if 100*(2*h-a) >= a]
    assert sum(a>18 for a,h,drop,w in actions) == 9967
    kappa = max(w/drop for a,h,drop,w in actions if a<=18 or drop>=3)
    assert kappa == F(345,64)
    deficits = {1:F(0),2:F(0)}
    ranks = {width:[F(0)]*(width+1) for width in (17,18)}
    for a,h,drop,w in actions:
        prepaid = min(w,kappa*drop)
        if w>prepaid:
            assert a>18 and drop in deficits
            deficits[drop] = max(deficits[drop],w-prepaid)
        for width,row in ranks.items():
            index = width*h//a
            row[index] = max(row[index],(kappa*(a-h)+prepaid)/h)
    assert deficits == {1:F(345,32),2:F(85,16)}
    total = []
    for width,row in ranks.items():
        for j in range(width-1,0,-1):
            row[j] = max(row[j],row[j+1])
        assert all(0<r<=kappa for r in row[1:])
        total.append(sum(row[1:]))
    assert total[0]<F('66.114') and total[1]<F('70.254')

    lam = dot((289,290,83810),direction)
    assert 0<lam<1
    path = [[F(f*g)/(1-lam) for g in READERS] for f in WIDTHS]
    prices = fractions('3.42777e-10 6.04273e-9 2.033e-11')
    clean,ae,am,budget = fractions('1e-25 .25 .00060 1.70')
    s = fractions('9.12647299934e-14 1.62466321632e-12 2.04235963738e-11 '
                  '2.89987253411e-11 2.60850931441e-12 1.86496982662e-13')
    p = fractions('2.40080729654e-6 2.15851458905e-7 2.51888072083e-8')
    c = fractions('2.87519951208e-13 2.73687632682e-14 3.04813760155e-14')
    u = fractions('4.10735503212e-10 7.31177163132e-9 2.06366984074e-11')
    children = [dot(row[3:],p) for row in path]
    z = [dot(row,s)+l for row,l in zip(m,children,strict=True)]
    v,d = max(u[:2]),u[2]
    readers = 10133*v+306*d
    witness = 28200*v+40*d+5508*readers**2/2
    from_measure = 800*v+2766309*v*witness+306*10133*v*readers
    for i,(bound,single) in enumerate(zip((witness,witness,from_measure),
                                          (clean,clean,prices[2]),strict=True)):
        assert u[i] > prices[i]*exp(z[i]+dot(n[i],u))
        assert s[i] > u[i]*bound+single*exp(z[i])
    for j,price in enumerate(direction):
        i = 3+j
        before = dot(path[i],(*s[:3],*c))
        growth = exp(z[i])
        assert s[i] > price*before*growth
        assert p[j] > price*(growth-1-children[i])
        assert c[j] > price*before*(growth-1)
    rv,rm = max(18*s[0]+s[3]+17*s[5],18*s[1]+s[4]+18*s[5]),s[2]
    assert 135*(46870021*rv+1420146*rm)<ae and 2766309*rv+83810*rm<am
    price_v = odds(TARGET)*exp(F('.015')+ae)
    for alpha,price in zip(total,prices[:2],strict=True):
        assert price_v*two_upper(ceil(SCALE*(alpha+budget))) < price
    assert price_v*two_upper(ceil(SCALE*budget))<clean
    assert odds(Q)*exp(F('.015')+am)<prices[2]
    assert rv<F('3.53e-11') and rm<F('2.05e-11')
    return rv,rm


def verify():
    m,n = matrices()
    old.group_rows(n)   # Revalida la supersolución de 18 coordenadas reutilizada.
    for cert in (GENERAL,PARTIAL):
        debt = check_local(cert)
        data_price = F('1.032')*cert.rate
        check_groups(cert,n,data_price)
        rv,rm = check_exterior(cert,*check_trees(cert,m),data_price)
        if cert == GENERAL:
            assert not any(debt.values())
            ref = odds(cert.rate)/odds(F('1e-6'))*F('1.907e-6')
            assert rv+ref < F('2e-16') and rm-odds(Q)*(1+F('.015'))+ref < F('2.6e-13')
        else:
            assert odds(cert.rate) > 2*odds(TARGET)
            assert debt[1] < F('2.606') and debt[2] < F('3.290')
            assert F('.3')*max(debt.values()) < 1
            older = earlier.previous
            small_reference = F('1.907e-6')*odds(TARGET)/odds(F('1e-6'))
            sectors = (older.two_data(),older.activity_counts(),earlier.directional_trees(),
                       earlier.spatial_trees(),(rv,rm),first_use_rows(m,n))
            sums = [sum(row)+small_reference for row in zip(*sectors,strict=True)]
            assert sums[0] < F('1.94e-5') and sums[1] < F('.00070')
            assert max(sums)+F('.004298') < F('.004998')
            print('Recargos máximos para descensos 1 y 2 sin luz inicial propia:',
                  *(f'{float(cost):.7g}' for cost in debt.values()),'bits.')
    # Una corrección física mixta de descenso uno que pasa el filtro.
    edge = {(0,j) for j in range(1,17)} | {(i,0) for i in range(1,18)}
    lights = {(0,j) for j in range(1,9)} | {(i,0) for i in range(1,10)}
    measures = edge-lights
    assert edge ^ measures == lights and len(lights)-len(measures) == 1
    assert min(a for a,d in shapes() if a>18) == 33
    for state in (lights,measures):
        assert max(sum(i==k for i,j in state) for k in range(18)) <= 8
        assert max(sum(j==k for i,j in state) for k in range(17)) <= 9
    assert 2*len(measures)<33 and F(1,33)>=F(1,100)
    print('PASS general: p <= 1.2e-34, q <= 2e-11. Resto de la referencia '
          '< 2e-16 por dato y < 2.6e-13 por medida en esa esquina.')
    print('PASS parcial en p=1e-30, q=2e-11 con un bit por dato en la cuenta de '
          'ocurrencias o 1.70 bits por dato en la cuenta de primeros usos. '
          'La suma completa en esa esquina sigue abierta.')


if __name__ == '__main__':
    verify()
