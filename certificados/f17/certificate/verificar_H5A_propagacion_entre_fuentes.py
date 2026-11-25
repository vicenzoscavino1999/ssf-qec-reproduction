"""Cotas de H5A_propagacion_entre_fuentes.txt. Python 3.10 o posterior."""

from fractions import Fraction as F
from functools import cache
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verificar_H5A_comparacion_de_ejecuciones as previous

exp = previous.exp_upper
dot = previous.prior.escalar
P, Q = previous.P, previous.Q
odds = previous.odds
LAMBDA = tuple(map(F, ('.003', '.00016', '.00000016')))


@cache
def bits(prob):
    """Menor múltiplo w de 1/64 con 2**(-w) <= prob, por enteros."""
    assert 0 < prob < 1
    num, den = prob.numerator**64, prob.denominator**64
    shift = max(0, den.bit_length()-num.bit_length())
    if (num << shift) < den:
        shift += 1
    assert (num << shift) >= den and (not shift or (num << (shift-1)) < den)
    return F(shift, 64)


def score_profile():
    singles = [(a,h,bits(rate*h)) for a,rate in zip((17,18),LAMBDA[:2],strict=True)
               for h in range(a//2+1,a+1)]
    kappa = max(w/(2*h-a) for a,h,w in singles)
    assert kappa == F(167,32)
    ranks = {n: [max((kappa*(a-h)+w)/h for a,h,w in singles if n*h >= k*a)
                 for k in range(1,n+1)] for n in (17,18)}
    assert all(0 < v <= kappa for row in ranks.values() for v in row)
    assert sum(ranks[17]) < F('63.4') and sum(ranks[18]) < F('68.5')
    # Formas posibles de un generador, sin exigir que todos sus estados
    # se realicen. Sólo interesan los primeros usos de un múltiple.
    areas = {a for x in range(19) for y in range(18)
             if x+y >= 2 and 2*(a := 17*x+18*y-2*x*y) >= 17*(x+y)}
    assert min(areas) == 33
    deficits = {}
    capped_actions = list(singles)
    low_deficits = {}
    for a in areas:
        for h in range(a//2+1,a+1):
            delta = 2*h-a
            if 100*delta < a:
                continue
            cost = bits(h*LAMBDA[2])
            allowed = min([cost,kappa*delta] +
                          [h*ranks[n][n*h//a-1]-kappa*(a-h) for n in (17,18)])
            assert allowed >= 0
            deficits[a,delta] = cost-allowed
            prepaid = min(cost, kappa*delta)
            capped_actions.append((a,h,prepaid))
            if cost > prepaid:
                low_deficits[a,delta] = cost-prepaid
    assert len(deficits) == 9967
    assert max(deficits.values()) == deficits[33,1] == F(425,32)
    assert all(not deficit for (a,d),deficit in deficits.items() if a >= 62 and d >= 4)
    # Al incluir todos los múltiples en la tabla por rangos, sólo queda
    # por pagar la parte de su primer uso que excede kappa*delta.
    capped_total = {}
    for n in (17,18):
        row = [max((kappa*(a-h)+prepaid)/h for a,h,prepaid in capped_actions
                   if n*h >= k*a) for k in range(1,n+1)]
        assert all(0 < b <= kappa for b in row)
        capped_total[n] = sum(row)
    assert capped_total[17] == F(664916069,10250240)
    assert capped_total[18] == F(920606373,13404160)
    assert capped_total[17]+F('1.5') < F('66.4')
    assert capped_total[18]+F('2.8') < F('71.5')
    assert {d for a,d in low_deficits} == {1,2,3}
    for delta, bound in ((1,F(425,32)), (2,F(255,32)), (3,F(11,4))):
        assert max(r for (a,d),r in low_deficits.items() if d == delta) == bound
    # Saldo para los pasos disjuntos de las huellas originales. Se incluyen
    # también las áreas 17 y 18 y todos los h admitidos por el umbral.
    sig = [(a,h) for a in areas | {17,18} for h in range(a//2+1,a+1)
           if 100*(2*h-a) >= a]
    remote = [sum(max(F(a-h,h) for a,h in sig if n*h>=a*k) for k in range(1,n+1))
              for n in (17,18)]
    assert remote[0] < F('11.110') and remote[1] < F('11.811')
    return len(deficits)


def directional_trees():
    matrix, vv, vc, cv, cc = previous.geometry()
    order = (0,1,3,0,1,2)
    m = [[matrix[i][j] for j in order] for i in order]
    n = [row[:3] for row in m[:3]]
    n[2][2] = 0
    lam = dot((289,290,83810),LAMBDA)
    assert lam == F('.9268096') and lam < 1
    path = [[F(i*j)/(1-lam) for j in (18,17,1,18,17,306)]
            for i in (17,18,1,17,18,306)]
    xv, xc, xi, clean = map(F, ('1.4e-10','4.8e-9','2.032e-11','1e-25'))
    s = tuple(map(F, ('2.678568e-14','9.254615e-13','2.037923e-11',
                     '3.553184e-11','2.008307e-12','4.421444e-14')))
    p = tuple(map(F, ('3.121236e-6','1.769189e-7','9.146489e-9')))
    c = tuple(map(F, ('5.401492e-13','3.230496e-14','1.062242e-14')))
    u = tuple(map(F, ('1.611985e-10','5.569505e-9','2.053659e-11')))
    bridges = [dot(row,s) for row in m]
    children = [dot(row[3:],p) for row in path]
    z = [w+l for w,l in zip(bridges,children,strict=True)]
    for i,price in enumerate((xv,xc,xi)):
        assert u[i] > price*exp(z[i]+dot(n[i],u))
    v, d = max(u[:2]), u[2]
    readers = 10133*v+306*d
    witness_v = 28200*v+40*d+5508*readers**2/2
    witness_m = 800*v+2766309*v*witness_v+306*10133*v*readers
    for i,(w,single) in enumerate(zip((witness_v,witness_v,witness_m),
                                      (clean,clean,xi),strict=True)):
        assert s[i] > u[i]*w+single*exp(z[i])
    for j,price in enumerate(LAMBDA):
        i = 3+j
        predecessor = dot(path[i],(*s[:3],*c))
        growth = exp(z[i])
        assert s[i] > price*predecessor*growth
        assert p[j] > price*(growth-1-children[i])
        assert c[j] > price*predecessor*(growth-1)
    rv = max(18*s[0]+s[3]+17*s[5],18*s[1]+s[4]+18*s[5])
    rm = s[2]
    ae, am = F('.28'), F('.00030')
    assert 135*(vv*rv+vc*rm) < ae and cv*rv+cc*rm < am
    assert rv < F('3.68e-11') and rm < F('2.04e-11')
    price_v = odds(P)*exp(F('.015')+ae)
    assert price_v**10*2**664 < xv**10
    assert price_v**10*2**715 < xc**10
    assert price_v*8 < clean and odds(Q)*exp(F('.015')+am) < xi
    return rv, rm


def ball_sizes(r):
    series = (272**r-1)//271
    next_series = 1+272*series
    return 1+289*series,1+288*series,17*next_series,18*next_series


def spatial_trees():
    # Datos unidos a cinco pasos de L. Una medida activa se puede adscribir
    # a dos pasos y una pasiva en contacto se busca en sentido inverso a tres.
    bv5,bc5,vc5,cv5 = ball_sizes(5)
    bv6,bc6,_,_ = ball_sizes(6)
    vv5 = (bv5*bv6+vc5**2,cv5**2+bc5*bc6)
    bv2,bc2,vc2,cv2 = ball_sizes(2)
    dm2 = (bv2*vc2,bc2*cv2)
    md2 = bv2*cv2+bc2*vc2
    bv3,bc3,vc3,cv3 = ball_sizes(3)
    md3 = bv3*cv3+bc3*vc3
    _, _, vm1, _, mm1 = previous.geometry()
    assert vm1 == 1420146 and mm1 == 83810
    sv, sm = (odds(r)*exp(F('.015')) for r in (P,Q))
    cd, tm = F('1.1e-29'), F('2.031e-11')
    bd, bm = F('.025'), F('.000002')
    for ng,nm in zip(vv5,dm2,strict=True):
        assert ng*cd+135*vm1*tm < bd
        assert sv*exp(2*bd+nm*sm+ng*cd) < cd
    assert md3*cd+mm1*tm < bm and sm*exp(2*bm) < tm
    assert tm+md2*cd < F('2.04e-11')
    return cd, tm+md2*cd


def avalanche():
    # El recuento de capas usa enteros para profundidad arbitraria. La nota demuestra
    # que el orden global siempre termina volteando exactamente estas capas.
    for depth in (3,6,11,29,83):
        layers = [2*8**j*9**(depth-j) for j in range(depth+1)]
        assert all(8*a == 9*b for a,b in zip(layers,layers[1:]))
        assert F(sum(layers[1:]),layers[0]) == 8*(1-F(8,9)**depth)
        if depth >= 6:
            assert sum(layers[1:]) > 4*layers[0]

    # Ejemplo pequeño con redondeo. Cada vértice de una capa interior
    # tiene nueve incidencias en la capa anterior. Una arista es un check.
    layers = []
    offset = 0
    for size in (100,88,78,68,60,52,46,40,34):
        layers.append(list(range(offset,offset+size)))
        offset += size
    neighbors = [[] for _ in range(offset)]

    def edge(a,b):
        assert a != b
        neighbors[a].append(b)
        neighbors[b].append(a)

    for ring,odd in ((layers[0],True),(layers[-1],False)):
        for i,v in enumerate(ring):
            for shift in range(1,5):
                edge(v,ring[(i+shift) % len(ring)])
        if odd:
            for i in range(len(ring)//2):
                edge(ring[i],ring[i+len(ring)//2])
    for sources,targets in zip(layers,layers[1:]):
        spare = 8*len(sources)-9*len(targets)
        assert 0 <= spare <= len(sources) and spare % 2 == 0
        for j in range(0,spare,2):
            edge(sources[j],sources[j+1])
        outgoing = [v for j,v in enumerate(sources) for _ in range(8-int(j<spare))]
        # Las incidencias paralelas son checks con etiquetas distintas.
        for j,v in enumerate(outgoing):
            edge(v,targets[j % len(targets)])
    assert all(len(row) == 17 for row in neighbors)
    spins = [i < len(layers[0]) for i in range(offset)]
    hits = [sum(spins[v] != spins[w] for w in row) for v,row in enumerate(neighbors)]
    steps = 0
    while max(hits) >= 9:
        v = hits.index(max(hits))
        assert not spins[v]
        for w in neighbors[v]:
            hits[w] += -1 if spins[w] else 1
        spins[v] = True
        hits[v] = 17-hits[v]
        steps += 1
    assert steps == 466 and all(spins) and not any(hits)
    print('Ejemplo físico por capas: 100 datos originales, 466 acciones, síndrome final cero.')


def verify():
    signatures = score_profile()
    dynamic = directional_trees()
    spatial = spatial_trees()
    low = previous.two_data()
    old_actions = previous.activity_counts()
    ref = F('1.907e-6')*odds(P)/odds(F('1e-6'))
    total = [sum(row)+ref for row in zip(dynamic,spatial,low,old_actions,strict=True)]
    assert total[0] < F('1.94e-5') and total[1] < F('.00070')
    assert max(total)+F('.004298') < F('.004998')
    avalanche()
    print(f'{signatures} firmas de acción múltiple. Máximo recargo 425/32 bits.')
    print('Segundo cómputo: sólo hay recargo para descensos de 1, 2 o 3 en el primer uso.')
    print('A p=1e-30 y q=2e-11, las filas adicionales por criterio son <3.68e-11.')
    print('Total cubierto por dato <1.94e-5 y por medida <0.00070. '
          'La cota completa en la meta sigue abierta.')


if __name__ == '__main__':
    verify()
