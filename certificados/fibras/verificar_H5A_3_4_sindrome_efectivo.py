"""Síndrome efectivo y acceso uniforme en HGP (3,4). Python 3.10+, biblioteca estándar."""

from collections import Counter
from fractions import Fraction as Q
from math import ceil, factorial
from pathlib import Path
from runpy import run_path


def previous(stem):
    return run_path(str(Path(__file__).with_name(f"verificar_H5A_3_4_{stem}.py")))


propagation = previous("propagacion_lejana")
dot, N = propagation["dot"], propagation["N"]
W, B = propagation["W"], propagation["B"]
ROWS, COLS = propagation["ROWS"], propagation["COLS"]
exp_upper = previous("medidas_y_un_dato")["exp_upper"]


def upper(x):
    return Q(ceil(10**24*exp_upper(x)), 10**24)


def effective_budgets(winners, stable):
    """Dos costes por subconjunto de una huella. Se conservan todos los empates."""
    cuts = set()
    for bits in range(1, 128):
        mask = 0
        for index, stripe in enumerate(ROWS+COLS):
            if bits >> index & 1:
                mask ^= stripe
        if mask and 2*mask.bit_count() >= 3*bits.bit_count():
            cuts.add(mask)
    assert len(cuts) == 63
    denominator = 27720  # Múltiplo de 1,...,12, para comparar h/v sin redondeo.
    prices = [(mask, denominator//mask.bit_count()) for mask in cuts]
    scores = [max((mask & pending).bit_count()*price for mask, price in prices)
              for pending in range(4096)]
    no_move = 10000  # Un coste realizable es <= 12 pasos de coste <= 36.
    costs = {}
    for score in set(scores):
        discharge, excess = [no_move]*4096, [no_move]*4096
        for lit, footprint in winners:
            h, size = lit.bit_count(), footprint.bit_count()
            descent = 2*h-size
            assert 1 <= descent and h <= 3*descent
            if h*denominator < score*size:
                continue
            firsts = lit
            while firsts:
                excess[firsts] = min(excess[firsts], 3*descent-h)
                if score < denominator or stable[lit ^ firsts]:
                    discharge[firsts] = min(discharge[firsts], descent)
                firsts = (firsts-1) & lit
        costs[score] = discharge, excess
    min_discharge, min_excess = [0]*4096, [0]*4096
    for pending in range(1, 4096):
        discharge, excess = costs[scores[pending]]
        firsts = pending
        best_discharge = best_excess = no_move
        while firsts:
            rest = pending ^ firsts
            best_discharge = min(best_discharge, discharge[firsts]+min_discharge[rest])
            best_excess = min(best_excess, excess[firsts]+min_excess[rest])
            firsts = (firsts-1) & pending
        min_discharge[pending], min_excess[pending] = best_discharge, best_excess
    caps, credits = [], []
    for fragment in range(4096):
        n = fragment.bit_count()
        pending = fragment
        cap, credit = 0, 5*n  # Incluye el pendiente vacío.
        while pending:
            r = pending.bit_count()
            credit = min(credit, 5*(n-r)+min_excess[pending])
            if stable[fragment ^ pending]:
                cap = max(cap, 2*r-n-min_discharge[pending])
            pending = (pending-1) & fragment
        caps.append(cap)
        credits.append(credit)
    assert Counter(caps) == {0:208, 1:2274, 2:1614}
    assert all(caps[mask] == 0 for mask in cuts)
    cross = ROWS[0] ^ COLS[0]
    assert cross.bit_count() == 5 and cross in cuts
    assert [credits[mask] for mask in (ROWS[0], COLS[0], cross)] == [6,8,9]
    assert caps[ROWS[0] ^ ROWS[1] ^ ROWS[2]] == 0
    return caps


def pair_classes(winners, stable, caps):
    old = propagation["pair_counts"](winners, stable)
    assert old == (((576,54,9),(0,279,9)),((0,496,16),(384,72,20)))
    assert all(caps[row ^ col] == 0 for row in ROWS for col in COLS)
    # Las clases de solapamiento mixto tienen 9 y 16 destinos, respectivamente.
    return ((old[0][0], (9,279,0)), ((16,496,0), old[1][1]))


def discount(saving):
    price = Q(saving, 1000)
    return 1/sum(price**j/factorial(j) for j in range(6))


def sector_polynomials(Y, ball, trees, pairs):
    saving = ((12,118),(118,74))
    regular, physical = [Counter(), Counter()], [Counter(), Counter()]
    for i in range(2):
        for marked, table in enumerate((regular[i], physical[i])):
            price = ((lambda k: Y**(-k)), (lambda k: propagation["anchor"](ball,i,k)))[marked]
            table[1-i,i] = (0 if marked else Y)*discount((108,64)[i])
            assert 2*price(0) <= price(1) and 2*price(1) <= price(2)
            for j in range(2):
                table[2-i-j,i+j] = dot(pairs[i][j],[price(k) for k in range(3)])*discount(saving[i][j])
            for (nv,nc), count in trees[i].items():
                if nv+nc >= 3:
                    table[nv,nc] = count*price(2*(nv+nc)//3)
    return regular, physical


def evaluate(table, values):
    return sum(count*values[0]**nv*values[1]**nc for (nv,nc),count in table.items())


def check_old_obstruction(a, kernel, poly):
    """Sólo el subsistema necesario de §6/v7, con las dos clases mixtas corregidas."""
    def pair_terms(table, w):
        return evaluate({key:val for key,val in table.items() if sum(key)<=2},w)
    h = tuple(map(Q,("1.27875235e-7","2.38691211e-7")))
    d = tuple(map(Q,(".000222429809",".000797740688")))
    e = tuple(map(Q,(".000222429809",".001714820104")))
    m = [Q("1e-6")*upper(Q(rate,1000)) for rate in (504,628)]
    w = [m[i]*upper(dot([k[0] for k in kernel[0][i]],h)) for i in range(2)]
    ratios = []
    for i in range(2):
        X, Z = dot(W[i],d), dot(B[i],e)
        low, high = w[i]/a**(7+i), w[i]/a**(7+3*i)
        ratios.extend((pair_terms(poly[i],w)/h[i], upper(X)*(low+high*(upper(Z)-1))/d[i],
                       high*upper(X+Z)/e[i]))
    assert max(ratios) < 1 and e[1] < Q(".00172")
    print(f"G3, subsistema necesario de v7: mayor cociente {float(max(ratios)):.10f}, e_C<.00172.")


def full_tail(series, a, w, i, d, e):
    X, Z = dot(W[i],d), dot(B[i],e)
    low, high = w[i]/a**(7+i), w[i]/a**(7+3*i)
    total = upper(X)*(low+high*(upper(Z)-1))
    assert total < d[i] and high*upper(X+Z) < e[i]
    ten = sum(n*w[0]**nv*w[1]**nc/a**nu for (nv,nc,nu),n in series[i].items())
    tail = total-ten
    assert tail >= Q(11**9,factorial(10))*W[i][i]**10*low**11
    return tail


def sector_tail(trees, Y, w, i, d):
    scale = Q("13.58")
    assert scale**3 >= Y**-2
    total = scale*w[i]*upper(dot(N[i],d))
    assert total < d[i]
    ten = evaluate(trees[i], [scale*wj for wj in w])
    tail = total-ten
    assert tail >= Q(11**9,factorial(10))*N[i][i]**10*(scale*w[i])**11
    return tail


def certify(Y, s, halos, regular, pure, kernel, tail_for, witness):
    p, b, alpha, row_cap, tree_factor = map(Q, witness["rates"])
    q, eta, tilt = Q(".001"), Q(".0001"), Q(".015")
    m, u, h, physical = (tuple(map(Q, witness[k])) for k in ("m","u","h","physical"))
    theta, visible, rooted = map(Q, witness["outer"])
    dd,dg,gd,gg = kernel
    g = 40*physical[0]+27*physical[1]+rooted
    R = max(dot(row,u) for row in N)
    I = 280*u[0]+189*u[1]+588*g
    y, z = q*upper(tilt)/(1-2*q), p*upper(tilt+alpha)/(1-2*p)
    ell = 89*b+I
    eps = upper(ell)*(upper(ell)-1)
    assert eps < 1 and max(u) < eta and R < row_cap < alpha
    ratios = [y*(1+b)*(1+eps)/b]
    def log_bound(t):
        return t+t*t/(2*(1-t))
    c = Q(89,2)*log_bound(theta**2)
    guard, ordinary = log_bound(theta)-c, b+2*c-log_bound(theta)
    assert guard >= q and ordinary >= 0
    ratios.append(y*(1-theta)*upper(99*guard+90*ordinary+2*I)/theta)
    for i in range(2):
        interaction = 2*dot([row[i] for row in halos[1:]],(u[0],u[1],g))
        ratios.append(z*upper(eta+2*R+guard*(504,628)[i]+ordinary*halos[0][i]+interaction)/m[i])
    w = tuple(m[i]*upper(dot([k[0] for k in dd[i]],h)+dg[i][0]*visible) for i in range(2))
    depart = tuple(dot([k[1] for k in dd[i]],h)+dg[i][1]*visible for i in range(2))
    initial = [(210,252),(640,784),(441,513),(1344,1596)]
    trees = tuple(mi*tree_factor for mi in m)
    removed, tails = [], []
    for i in range(2):
        tail = tail_for(w,i)
        tails.append(tail)
        branches = evaluate(regular[i],w)+tail
        removed.append(evaluate(regular[i],m))
        pointed = sum(count*w[0]**nv*w[1]**nc*(nv*depart[0]+nc*depart[1])
                      for (nv,nc),count in regular[i].items())
        connected = upper(dot(N[i],trees))
        assert m[i]*connected < trees[i]
        fixed = min(branches-removed[i],pointed+tail)+m[i]*(connected-1)
        simple = z*upper(eta+2*R+dot((b,2*u[0],2*u[1],2*g),[row[i] for row in initial]))*(1+y)**(3+i)
        ratios.extend((branches/h[i],fixed/physical[i],(simple+physical[i])/u[i]))
    J = dot([k[0] for k in gd],h)+gg[0]*visible
    one = dot([k[1] for k in gd],h)+gg[1]*visible
    without_exterior = dot([k[1] for k in gd],removed)
    direct = sum(evaluate(tab,m) for tab in pure)
    outer = Y/s**2
    root = outer*(upper(J)*one-without_exterior)+direct
    ratios.extend((outer*(upper(J)-1-J)/visible,root/rooted))
    assert all(0 < ratio < 1 for ratio in ratios), [(j,float(v)) for j,v in enumerate(ratios) if not 0<v<1]
    print(f"{witness['name']}: p={float(p):g}, q=.001, fila {float(R):.12g} < {float(row_cap):g}.")
    print(f"  Cociente exterior máximo {float(max(ratios)):.12f}. "
          f"Colas VV, CC {float(tails[0]):.8g}, {float(tails[1]):.8g}.")


FULL = dict(name="Todos los conjuntos", rates=("3.25e-7",".001155",".001",".000589","1.003"),
            m=("6.40434396083e-7","7.57050851175e-7"), u=("4.38346037129e-7","4.77809404376e-7"),
            h=("5.47853759670e-8","2.89690689400e-7"), physical=("2.63443321284e-9","1.73708936335e-8"),
            outer=(".0011599599554","3.87235375512e-7","1.02865359088e-5"))
SECTOR = dict(name="Componentes con saldo acotado", rates=("1e-6",".00118",".004",".001992","1.02"),
              m=("2.29589564260e-6","2.82859013003e-6"), u=("1.48008203174e-6","1.62023462830e-6"),
              h=("4.08042787401e-7","9.07019646249e-7"), physical=("2.97813663199e-8","6.64226615640e-8"),
              outer=(".00119819220697","5.38214437031e-7","3.21382017618e-5"))


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para verificar las desigualdades.")
    contacts, striped, base, history, old = map(previous,("contacto_de_checks","franjas","medidas_y_un_dato",
                                                            "medidas_por_historia","medidas_aisladas"))
    profiles, stable, reps = contacts["geometry"]()
    winners = propagation["winners"](reps)
    caps = effective_budgets(winners, stable)
    pairs = pair_classes(winners, stable, caps)
    print("4096 síndromes: 208 con saldo <=0, 2274 con techo 1 y 1614 con techo 2. Las 63 huellas tienen techo cero.")
    halos = contacts["exposure"](reps,base["budgets"])
    ball = history["ball"]
    trees = striped["coefficients"](10)
    Y, s = Q("0.02"), Q("0.9")
    regular, pure = sector_polynomials(Y,ball,trees,pairs)
    sector_d = tuple(map(Q,("3.61827330841e-5","4.58463611394e-5")))
    certify(Y,s,halos,regular,pure,propagation["kernels"](ball,Y,s),
            lambda w,i: sector_tail(trees,Y,w,i,sector_d),SECTOR)
    _, prior = striped["cost_tables"](striped["light_profiles"](),previous("fondo_pasivo")["step_costs"],
                                        previous("pares_adyacentes")["stripe_profiles"])
    hist = contacts["positions"](profiles,stable,prior,striped["coefficients"](4))
    tables = history["data_coefficients"](hist)
    triples = old["count_triples"](profiles,stable,contacts["rooted_trees"])
    series = old["group_series"](10)
    a, s = Q(".5488178589126"), Q(".86")
    regular, pure = propagation["group_polynomials"](a**6,ball,tables,triples,series,pairs,
                         propagation["cover_profiles"](reps),striped["coefficients"](propagation["LIMIT"]))
    kernel = propagation["kernels"](ball,a**6,s)
    check_old_obstruction(a,kernel,regular)
    full_d = tuple(map(Q,("5.28010437671e-5",".000134329789883")))
    full_e = tuple(map(Q,("5.28010437671e-5",".000391752404757")))
    certify(a**6,s,halos,regular,pure,kernel,lambda w,i: full_tail(series,a,w,i,full_d,full_e),FULL)
