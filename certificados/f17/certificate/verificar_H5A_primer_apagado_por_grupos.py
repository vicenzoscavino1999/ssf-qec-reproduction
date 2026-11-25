"""Primer apagado por grupos. Certificado racional del apéndice H5-A, Python 3.10+."""

from dataclasses import replace
from fractions import Fraction as F
from functools import cache
from math import ceil, comb
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verificar_H5A_contactos_entre_ejecuciones as prior

P = F('2e-33')
Q = prior.Q
X = F('1e-8')
Y = prior.XI * prior.exp(prior.Z_CAP[2])
PRICE = F('1.032') * P * prior.exp(prior.Z_CAP[0])
N = 28
SCALE = 10**105
DV, DM, MD = 46870021, 1420146, 2766309
DUAL = prior.GENERAL.dual
RANK = prior.GENERAL.rebates
MICRO = prior.SCALE
CREDIT, EXCESS = 6990000, 421000
FIELDS = [(u,w) for u,w in zip(prior.LEVELS,DUAL[1:],strict=True) if w]


def ceiling_ratio(a,b):
    return -(-a//b)


def check_credit():
    bound = 562500
    assert all(max(r[:-1]) <= bound and r[-1]-bound <= EXCESS for r in RANK)
    # w(0) y w(h) en una corrección perfecta. Para h>=17 valen estas dos cotas.
    start = MICRO*prior.bits(prior.TOP)
    end = MICRO*prior.bits(prior.TOP+17*prior.DIRECT)
    assert max(start-17*bound,end) <= 172000
    assert (DUAL[0]-bound)*F(33,2)-172000 >= CREDIT
    # En el caso de un solo dato perfectamente retirado se usa menos de un bit.
    assert end < MICRO


@cache
def recurrence(widths):
    """Una fila inicial por anchura y las primeras bajadas de puntuación."""
    levels = [(tuple(ceil(n*(1+u)/2)-1 for n in widths),ceil(17*(1+u)/2),w)
              for u,w in FIELDS]

    @cache
    def visit(index,k,b):
        if index == len(levels):
            return 0
        caps,minimum,w = levels[index]
        best = 0
        for j in range(min(k,max(caps))+1):
            spent = max(minimum,ceiling_ratio(17*k,18)) if j<k else 0
            remain = min(sum(min(j,cap) for cap in caps),max(0,b-spent))
            best = max(best,w*remain+visit(index+1,j,remain))
        return best
    return visit


@cache
def fees(e):
    # Usar CC para los límites de las fases y el peor precio inicial de ambos tipos.
    visit = recurrence((18,)*e)
    initial = [e*max(DUAL[0]*min(k,n)-sum(row[:k]) for n,row in
                    zip((17,18),RANK,strict=True)) for k in range(19)]
    bounds = [max(initial[k]+visit(0,k,e*k) for k in range(19))]
    for mark in range(1,2*e+1):
        used = ceiling_ratio(33*mark,2)
        pre = max(initial[k]+visit(0,j,min(e*j,max(0,e*k-used)))
                  for k in (17,18) for j in range(min(k,17)+1))
        bounds.append(pre-max(0,CREDIT*mark-EXCESS*e))
    return bounds


def power(bits):
    """Majora 2 elevado a un exponente en millonésimas, también negativo."""
    whole,remainder = divmod(bits,MICRO)
    return F(2)**whole * prior.two_upper(remainder)


def fixed(value):
    return ceil(value*SCALE)


def mul_sum(left,right,degree,orders,weights=False):
    """Coeficiente de grado dado. Cada casilla se redondea una sola vez hacia arriba."""
    out = [0]*(2*degree+1)
    for j in orders:
        for i,l in enumerate(left[j]):
            if l:
                factor = j*l if weights else l
                for k,r in enumerate(right[degree-j]):
                    if r:
                        out[i+k] += factor*r
    divisor = degree*SCALE if weights else SCALE
    return [ceiling_ratio(v,divisor) for v in out]


def block_rows(n,y):
    def series():
        return [[0]*(2*j+1) for j in range(n+1)]
    data,meas,site,inside,bridges,internal,parts = [series() for _ in range(7)]
    exp_bridge,exp_inside,exp_meas,exp_single,exp_multi = [series() for _ in range(5)]
    lines,rectangles = series(),series()
    meas[0][0] = fixed(y)
    exp_bridge[0][0] = fixed(prior.exp(DM*F(meas[0][0],SCALE)))
    for unit in (exp_inside,exp_meas,exp_single,exp_multi):
        unit[0][0] = SCALE
    for e in range(1,n+1):
        site[e][:-2] = [ceiling_ratio(X.numerator*v,X.denominator) for v in exp_bridge[e-1]]
        extra = mul_sum(site,exp_inside,e,range(1,e))
        inside[e] = [a+b for a,b in zip(site[e],extra,strict=True)]
        internal[e] = [612*t for t in inside[e]]
        lines[e] = [630*t for t in site[e]]
        rectangles[e] = [10133*t for t in site[e]]
        exp_inside[e] = mul_sum(internal,exp_inside,e,range(1,e+1),True)
        exp_single[e] = mul_sum(lines,exp_single,e,range(1,e+1),True)
        exp_multi[e] = mul_sum(rectangles,exp_multi,e,range(1,e+1),True)
        data[e] = [a+ceil(18*y*b) for a,b in zip(extra,site[e],strict=True)]
        for mark,count,series_exp in ((1,613,exp_single),(2,5508,exp_multi)):
            companions = mul_sum(site,series_exp,e,range(1,e))
            for i,(solo,others) in enumerate(zip(site[e],companions,strict=True)):
                if i+mark<=2*e:
                    data[e][i+mark] += solo+count*others
                else:
                    assert not solo and not others
        parts[e] = [MD*v for v in data[e]]
        exp_meas[e] = mul_sum(parts,exp_meas,e,range(1,e+1),True)
        meas[e] = [ceil(y*v) for v in exp_meas[e]]
        bridges[e] = [DV*a+DM*d for a,d in zip(data[e],meas[e],strict=True)]
        exp_bridge[e] = mul_sum(bridges,exp_bridge,e,range(1,e+1),True)
    rooted_meas = [[d+ceil(35*y*b) for d,b in zip(di,si,strict=True)]
                   for di,si in zip(meas,site,strict=True)]
    return data,rooted_meas


def pair_fee(n,m,k1,k2,after=None):
    initial = DUAL[0]*(k1+k2)-sum(RANK[n-17][:k1])-sum(RANK[m-17][:k2])
    k,b = after if after is not None else (max(k1,k2),k1+k2)
    return initial+recurrence((n,m))(0,k,b)


def pure_pair(root_width):
    n = root_width
    result = F()
    for m in (17,18):
        for j in range(1,min(n,m)+1):
            if j==1 or n==m:
                ways = 630 if j==1 else 324//j
                result += ways*power(pair_fee(n,m,n-j,m-j))
        remain = max(0,n+m-33)
        worst = max(pair_fee(n,m,n,m,(k,remain)) for k in range(min(17,remain)+1))
        result += DV*power(worst)
        for width in (17,18):
            for i in range(1,width):
                j = width-i
                if (i>1 and n!=width) or (j>1 and m!=width):
                    continue
                ways = (630 if i==1 else 324//i)*(630 if j==1 else 324//j)
                fee = pair_fee(n,m,n,m,(max(n-i,m-j),n+m-width))
                result += ways*power(fee-CREDIT+2*EXCESS)
    return result*PRICE**2


def tail():
    # Medio bit de las filas de mayor puntuación va a los primeros pasos perfectos.
    left = 500000
    offsets = [0,0]
    for u,w in FIELDS:
        held = min(w,left)
        left -= held
        for i,n in enumerate((17,18)):
            cap = ceil(n*(1+u)/2)-1
            offsets[i] += w*cap+held*(n-cap)
    assert left==0
    alpha = max(DUAL[0]*n-sum(r)+v+EXCESS for n,r,v in
                zip((17,18),RANK,offsets,strict=True))
    x = 4*PRICE*power(alpha)
    mark = power(-CREDIT-8250000)
    a,t,d,site = map(F,('3e-10','6.6e-7','2.07e-11','6.5e-7'))
    assert x*prior.exp(DV*a+DM*d)<site and site*prior.exp(612*t)<t
    one = 1+613*(prior.exp(630*site)-1)
    two = 1+5508*(prior.exp(10133*site)-1)
    assert site*(prior.exp(612*t)-1+18*Y+mark*one+mark**2*two)<a
    assert Y*prior.exp(MD*a)<d
    return a/4**(N+1),(d+35*Y*site)/4**(N+1)


def source_rows():
    data,meas = block_rows(N,Y)
    inflated = F('1e-7')
    high,_ = block_rows(3,inflated)
    _,n = prior.matrices()
    measured = F()
    roots = []
    for i,width in enumerate((17,18)):
        factor = PRICE*prior.exp(n[i][2]*Y)
        row = RANK[i]
        by_k = [DUAL[0]*k-sum(row[:k])+recurrence((width,))(0,k,k) for k in range(width+1)]
        roots.append(factor*sum(power(by_k[width-m])*comb(width,m)*Y**m for m in range(width+1)))
        measured += 35*factor*sum(power(by_k[width-m])*comb(width-1,m-1)*Y**m for m in range(1,width+1))
    measured += MD*Y*sum(roots)
    for e in range(2,N+1):
        cost = fees(e)
        rates = [power(c)*(PRICE/X)**e/SCALE for c in cost]
        sums = [v*w for v,w in zip(data[e],rates,strict=True)]
        measured += sum(v*w for v,w in zip(meas[e],rates,strict=True))
        if e==2:
            with_measure = (Y/inflated)*sum(v*w for v,w in zip(high[e],rates,strict=True))
            for i,width in enumerate((17,18)):
                roots[i] += min(sum(sums),with_measure+pure_pair(width))
            continue
        if e==3:
            # Un único centro singleton. Si en él sólo participa su dato, se retira
            # ese apagado perfecto antes de aplicar el perfil conservador del par.
            layouts = 4*612*613*630+613*630**2
            paired = 4*612*DV*PRICE**3*power(max(fees(2))+MICRO)
            restricted = (Y/inflated)*high[e][1]*rates[1]+layouts*PRICE**e*power(cost[1])+paired
            sums[1] = min(sums[1],restricted)
        value = sum(sums)
        for i in range(2):
            roots[i] += value
    rest_data,rest_meas = tail()
    return tuple(v+rest_data for v in roots)+(measured+rest_meas,)


def verify():
    check_credit()
    prior.verify()  # También reconstruye el presupuesto del resto en la meta.
    actual = source_rows()
    sources = prior.fractions('8.85e-24 1.19e-23 1.49e-26')
    assert all(a<b for a,b in zip(actual,sources,strict=True)), [float(a) for a in actual]
    branches = prior.fractions('8.73004832254e-19 8.73005304843e-19 1.22776475412e-17 '
                              '1.26412979288e-11 1.26096663684e-11 1.93016493031e-10 '
                              '8.03421364794e-24 8.5068024667e-24 2.03362045923e-21 '
                              '3.12895359339e-22 3.31301148043e-22 5.63299759173e-21 '
                              '4.68398440671e-7 4.95950578392e-7 8.43158697052e-6')
    cert = replace(prior.GENERAL,rate=P,sources=sources,branches=branches)
    matrix,_ = prior.matrices()
    v,m = prior.check_trees(cert,matrix)
    rv,rm = v+18*F('1.045')*P,m+prior.XI
    lv,lm = F('1.045')*P,prior.XI
    assert lv>prior.odds(P)*prior.exp(F('.025')+DV*(lv+v)+DM*(lm+m))
    assert lm>prior.odds(Q)*prior.exp(F('.025')+MD*(lv+v)+83810*(lm+m))
    assert prior.XI>prior.odds(Q)*prior.exp(F('.025')+MD*lv+83810*lm)
    assert F('1.032')*P>prior.odds(P)*prior.exp(F('.025')+135*(DV*lv+DM*lm))
    assert max(135*(DV*rv+DM*rm),MD*rv+83810*rm)<F('.01')
    reference = F('1.907e-6')*prior.odds(P)/prior.odds(F('1e-6'))
    assert rv+reference<F('2.4e-16')
    assert rm-prior.odds(Q)*(1+F('.015'))+reference<F('2.6e-13')
    # Para un enlace directo en el rectángulo de un múltiple, la fila
    # geométrica usa ya sus 306 celdas, también si el check estaba apagado.
    assert prior.WIDTHS[5]==prior.READERS[5]==18*17
    debt = prior.check_local(prior.PARTIAL)
    assert debt[1]<F('2.606') and debt[2]<F('3.290')
    print('Filas fuente VV, CC y medida:',*(f'{float(v):.8g}' for v in actual))
    print('PASS general: p <= 2e-33, q <= 2e-11. Resto completo en la esquina '
          '< 2.4e-16 por dato y < 2.6e-13 por medida.')
    print('PASS local en p=1e-30: los recargos sólo se necesitan en múltiples '
          'de descenso 1 o 2 cuyo rectángulo no lee ningún dato original del núcleo.')


if __name__=='__main__':
    verify()
