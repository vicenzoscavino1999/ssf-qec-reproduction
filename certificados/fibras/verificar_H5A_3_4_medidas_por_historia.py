"""HGP (3,4). Certificado uniforme en (1.2e-7, .001) con la biblioteca estándar."""

from collections import Counter
from fractions import Fraction as Q
from math import ceil, factorial, floor
from pathlib import Path
from runpy import run_path


WEAK = ((486, 324), (576, 384))
STRONG = ((154, 117), (208, 129))
NEIGHBORS = ((640, 441), (784, 513))


def far_credits(profiles, stable):
    # Unidades de 1/6 del descenso. Sólo el rango tercero del CC aislado
    # en el grafo fuerte recibe un incremento respecto del crédito ordinario.
    generic = ((2, 3, 6), (2, 2, 4, 6))
    isolated_cc = (2, 2, 6, 6)
    for lit, h, v, options in profiles:
        limit = 6*(2*h-v)
        maximum = max(max(row[:n*h//v]) for n, row in zip((3, 4), generic))
        choices = [(mask, count, blocked) for mask, count, blocked in options[1] if mask]
        def place(start, used, blocked, credit, pristine):
            if not pristine or stable[lit & blocked & ~used]:
                assert credit + maximum*(lit & ~blocked).bit_count() <= limit
            for j in range(start, len(choices)):
                mask, count, hit = choices[j]
                if mask & blocked:
                    continue
                for rank in range(count, 4*h//v+1):
                    share = sum(isolated_cc[rank-count:rank])
                    place(j+1, used|mask, blocked|hit, credit+share, pristine or rank == 4)
        place(0, 0, 0, 0, False)
    for row, expected in zip(generic, (7, 10)):
        assert max(6*r-sum(row[:r]) for r in range(len(row)+1)) == expected
    assert max(6*(2*r-4)-sum(isolated_cc[:r]) for r in range(5)) == 8


def data_coefficients(old_tables):
    result = [Counter(), Counter()]
    for root, table in enumerate(old_tables):
        for (nv, nc, old_bulk, old_small), multiplicity in table.items():
            # El histograma anterior guarda el coste a granel en sextos.
            # En grados hasta cuatro identifica unívocamente los dos conteos.
            fits = [(v, c) for v in range(nv+1) for c in range(nc+1)
                    if 84*nv+118*nc-11*v-28*c == old_bulk]
            assert len(fits) == 1
            _, isolated_cc = fits[0]
            bulk = 7*nv+10*nc-2*isolated_cc
            small = min(bulk//6, max(0, (old_small-1)//66)) if nv+nc >= 3 else 0
            result[root][nv, nc, bulk, small] += multiplicity
    return result


def ball(variable, radius):
    if radius % 2:
        return (3 if variable else 4)*(6**((radius+1)//2)-1)//5
    return 1+(9 if variable else 8)*(6**(radius//2)-1)//5


def geometric_rows(Y):
    def kernel(v1, r1, v2, r2):
        def size(k):
            return ball(v1, r1+2*k)*ball(v2, r2+2*k)
        a0, a1, a2 = (size(k) for k in range(3))
        c36 = Q(a2-7*a1+6*a0, 1050)
        c6 = Q(a1-a0-35*c36, 5)
        c1 = a0-c6-c36
        for k in range(9):
            assert size(k) == c36*36**k+c6*6**k+c1
            assert size(k) <= c36*36**k
        # Cada factor de la bola tiene forma A*6^k-B, con 0 <= B < A.
        assert c6 <= 0 and c6+c1 <= 0
        weighted = c1+(1-Y)*(c6/(1-6*Y)+c36/(1-36*Y))
        return ceil(weighted), ceil(c36)
    dd = ((kernel(True,4,True,6)[0], kernel(True,5,True,5)[0]),
          (kernel(False,5,False,5)[0], kernel(False,6,False,4)[0]))
    dg = kernel(True,3,True,4)[0], kernel(False,4,False,3)[0]
    gv, gc, gg = kernel(False,3,True,4), kernel(False,4,True,3), kernel(False,2,True,2)
    return dd, dg, (gv[0], gc[0]), gg[0], (gv[1], gc[1], gg[1])


def obstruction():
    """Condición necesaria de las dos ramas a granel, para todo 0 < Y < 1/36."""
    p, q = Q("8e-7"), Q(".001")
    base = Q("1.817")
    assert base**6 < 36
    def lower_exp(x):
        return sum(x**j/factorial(j) for j in range(7))
    def down(x):
        return Q(floor(x*10**12), 10**12)
    m = [p*lower_exp(2*q*k) for k in (414, 536)]
    d, e = [Q(0)]*2, [Q(0)]*2
    for step in range(1, 61):
        x = [sum(k*s for k,s in zip(row,d)) for row in WEAK]
        z = [sum(k*s for k,s in zip(row,e)) for row in STRONG]
        d = [down(lower_exp(x[i])*m[i]*(base**(7+i)
                 +base**(7+3*i)*(lower_exp(z[i])-1))) for i in range(2)]
        e = [down(m[i]*base**(7+3*i)*lower_exp(x[i]+z[i])) for i in range(2)]
        if e[1] >= Q(".1"):
            break
    assert e[1] >= Q(".1") and step == 25
    assert m[1]*base**10*(129*Q(".1"))**3/6 > Q(".1")
    return step


def small_group_obstruction(b, Y, gd, dg, dd, gg, halos, pointed, pairs_and_triples):
    """Condición necesaria en 2e-7 con pesos fijos, sólo con parejas y tríos."""
    p,q = Q("2e-7"), Q(".001")
    def lower_exp(x):
        return sum(x**j/factorial(j) for j in range(7))
    def down(x):
        return Q(floor(x*10**18), 10**18)
    y = q*(1+Q(".015"))/(1-2*q)
    g, visible, h = Q(0), Q(0), (Q(0), Q(0))
    for step in range(1, 17):
        I = 588*g
        ell = 89*b+I
        lower_ratio = y*(1+b)*(1+ell+3*ell**2/2+7*ell**3/6)/b
        if lower_ratio > Q("1.01"):
            return step
        s = down(y*lower_exp(90*b+2*I))
        assert s < 1
        m = [down(p*lower_exp(b*k+2*kg*g)/(1-s)**k) for k,kg in zip(halos[0],halos[3])]
        x = [down(m[i]*lower_exp(sum(k*v for k,v in zip(dd[i],h))+dg[i]*visible))
             for i in range(2)]
        h = [down(pairs_and_triples(i,x)) for i in range(2)]
        J = sum(k*v for k,v in zip(gd,h))+gg*visible
        g = down(Y*(pointed[0]*h[0]+pointed[1]*h[1]+pointed[2]*visible))
        visible = down(Y*(lower_exp(J)-1-J))
    raise AssertionError("No se alcanzó la obstrucción de las filas pequeñas.")


def certify(tables, halos, exp_upper):
    p,q,b,alpha,eta = map(Q, ("1.2e-7", ".001", ".00118", ".001", ".0001"))
    r = Q(".9889")
    a = r**55
    Y = a**6
    assert 0 < 36*Y < 1
    def upper(x):
        return Q(ceil(exp_upper(x)*10**24), 10**24)
    dd,dg,gd,gg,pointed = geometric_rows(Y)
    assert dd == ((257772,171947),(305683,203723)) and dg == (14280,16928)
    assert gd == (19040,12696) and gg == 1047 and pointed == (1867,1245,104)
    assert halos == [(414,536),(1264,1652),(Q(1719,2),1095),(Q(5247,2),3374)]
    initial = [(210,252),(640,784),(441,513),(1344,1596)]

    # Supersolución racional. h cuenta ramas del árbol y physical sólo
    # núcleos completos no simples con un original fijado.
    s,theta = Q(".0011503968748685456"), Q(".001290508639749")
    C = tuple(map(Q, ("2.1438646442233673e-7","2.529513595658138e-7")))
    u = tuple(map(Q, ("1.7912020068799636e-7","2.8603865073831506e-7")))
    h = tuple(map(Q, ("4.206720379938035e-8","2.494542132730112e-7")))
    physical = tuple(map(Q, ("1.615402295793722e-8","1.1336892907247649e-7")))
    d = tuple(map(Q, ("3.086522248467517e-5","8.864348485088505e-5")))
    e = tuple(map(Q, ("3.086522248467517e-5",".0002753185073549987")))
    visible, rooted = Q("2.22084926794114e-7"), Q("1.0405602401557564e-5")
    g = 40*physical[0]+27*physical[1]+rooted
    row = max(sum(k*ui for k,ui in zip(line,u)) for line in NEIGHBORS)
    I = 280*u[0]+189*u[1]+588*g
    ell = 89*b+I
    y = q*upper(Q(".015"))/(1-2*q)
    z = p*upper(Q(".015")+alpha)/(1-2*p)
    eps = upper(ell)*(upper(ell)-1)
    assert eps < 1 and max(u) < eta and row < Q(".00028717") < alpha
    assert row*Q(3,8)**3 < Q(".00001515")
    price = [y*(1+b)*(1+eps)/b, y*upper(90*b+2*I)/s]
    assert 0 < s <= theta*(1-theta)**89
    factors = (b,2*u[0],2*u[1],2*g)
    def data_price(table, i):
        return z*upper(eta+2*row+sum(w*line[i] for w,line in zip(factors,table)))
    m = tuple(Q(ceil(10**24*C[i]/(1-theta)**halos[0][i]),10**24) for i in range(2))
    x = tuple(m[i]*upper(sum(k*v for k,v in zip(dd[i],h))+dg[i]*visible) for i in range(2))
    pair = ((576+64/Y+9/Y**2,441/Y+9/Y**2),
            (784/Y+16/Y**2,384+72/Y+65/Y**2))
    tree = [mi*Q("1.003") for mi in m]
    def order(i,w,size):
        return sum(n*w[0]**nv*w[1]**nc/Y**k
                   for (nv,nc,_,k),n in tables[i].items() if nv+nc == size)
    def pairs_and_triples(i,w):
        return w[i]*sum(k*v for k,v in zip(pair[i],w))+order(i,w,3)
    def first_four(i,w):
        return Y*w[i]+pairs_and_triples(i,w)+order(i,w,4)
    tails = []
    for i in range(2):
        price.append(data_price(halos,i)/C[i])
        X = sum(k*v for k,v in zip(WEAK[i],d))
        Z = sum(k*v for k,v in zip(STRONG[i],e))
        good, high = x[i]/a**(7+i), x[i]/a**(7+3*i)
        full = upper(X)*(good+high*(upper(Z)-1))
        price.extend((full/d[i], high*upper(X+Z)/e[i]))
        initial_four = sum(n*x[0]**nv*x[1]**nc/a**bulk
                           for (nv,nc,bulk,_),n in tables[i].items())
        tail = full-initial_four
        assert tail >= Q(125,24)*WEAK[i][i]**4*good**5
        tails.append(tail)
        all_groups = first_four(i,x)+tail
        connected_exp = upper(sum(k*ui for k,ui in zip(NEIGHBORS[i],tree)))
        assert m[i]*connected_exp <= tree[i]
        physical_groups = all_groups-first_four(i,m)+m[i]*(connected_exp-1)
        simple = data_price(initial,i)*(1+y)**(3+i)
        price.extend((all_groups/h[i], physical_groups/physical[i], (simple+physical[i])/u[i]))
    J = sum(k*v for k,v in zip(gd,h))+gg*visible
    one_edge = pointed[0]*h[0]+pointed[1]*h[1]+pointed[2]*visible
    price.extend((Y*(upper(J)-1-J)/visible, Y*upper(J)*one_edge/rooted))
    assert 0 < min(price) and max(price) < 1
    print(f"p=1.2e-7, q=.001. Fila completa {float(row):.12g} < .00028717 frente a .001.")
    print(f"Cociente pasivo {float(price[0]):.12f}. Mayor cociente comprobado {float(max(price)):.12f}.")
    print(f"Colas desde cinco datos por raíz VV y CC: {float(tails[0]):.10g}, {float(tails[1]):.10g}.")
    step = small_group_obstruction(b,Y,gd,dg,dd,gg,halos,pointed,pairs_and_triples)
    print(f"A p=2e-7 con estos pesos, sólo parejas y tríos imponen un cociente pasivo >1.01 "
          f"en {step} rondas de cotas inferiores.")
    print(f"A p=4.5e-8, tres o más componentes < .00001515. Obstrucción desde p=8e-7 "
          f"en {obstruction()} iteraciones para cualquier Y admisible de estas filas.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para comprobar las desigualdades.")
    folder = Path(__file__).resolve().parent
    def previous(name):
        return run_path(str(folder/f"verificar_H5A_3_4_{name}.py"))
    contacts, striped = previous("contacto_de_checks"), previous("franjas")
    profiles, stable, representatives = contacts["geometry"]()
    far_credits(profiles,stable)
    _, caps = striped["cost_tables"](striped["light_profiles"](),previous("fondo_pasivo")["step_costs"],
                                     previous("pares_adyacentes")["stripe_profiles"])
    tables = contacts["positions"](profiles,stable,caps,striped["coefficients"](4))
    base = previous("medidas_y_un_dato")
    certify(data_coefficients(tables),contacts["exposure"](representatives,base["budgets"]),base["exp_upper"])
