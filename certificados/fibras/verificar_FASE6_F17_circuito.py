"""Transferencia al circuito de una ronda (17,18). Python 3.10+, biblioteca estándar."""
from dataclasses import replace
from fractions import Fraction as F
from hashlib import sha256
from math import ceil
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile
import sys


PHEN_P, RM, RC = map(F, ("5e-33", "2e-11", "2e-35"))
BASE = 34*RC
FULL, DEFECT, MANY = F(".704"), F(".0052"), F("1e-30")
ARCHIVE = "SOURCE_CUMULANT_CERTIFICATE.zip"
DIGEST = "a83972e07afc351d39cab0046d7a81c225996ba0e4c77d4413532621d906c328"


def finite_f17():
    # Dos lados del HGP de la matriz de 289 filas por 306 columnas.
    a, b = 306, 289
    n, m = a*a+b*b, a*b
    layers, edges = 35+18-1, 35*m
    locations = (layers+2)*(n+m)-edges
    assert (n, m, locations) == (177157, 88434, 11246724)
    op = PHEN_P/(1-PHEN_P)
    read = 1/((1-PHEN_P)*(1-RM))-1
    total = (locations-m)*op+m*read
    error = 2*total/(1-2*total)
    assert total < F("1.769e-6") and error < 1
    norm = error/(1-error)
    assert norm < F("3.538e-6")
    return total, norm


def maximum_timing_weight(n, u, v, price):
    # Mejor prefijo para cada número de CNOT ya leídas. La puerta h-ésima
    # por capa está a más tardar en 35+h-1. Las dos esperas de borde son fijas.
    states = {0:u*price[0]}
    for layer in range(1, 53):
        following = {j:weight+u*price[j] for j, weight in states.items()}
        for j, weight in states.items():
            if j < n and layer <= 35+j:
                gate = weight+v*(price[j]+price[j+1])/2
                following[j+1] = max(following.get(j+1, F(0)), gate)
        states = following
    return states[n]+u*price[n]


def check_fibre_laws(first):
    u, v = 2*RC/(3-2*RC), 8*RC/(15-8*RC)
    target = 4*RC/(15-8*RC)
    # La factorización de la CNOT reproduce los cuatro pesos marginales X.
    a = 8*RC/15
    assert ((1-a)*(1-target), (1-a)*target, a/2, a/2) == (1-4*RC/5, 4*RC/15, 4*RC/15, 4*RC/15)
    discount = first.power(-sum(first.DUAL))
    for n in (17, 18):
        total = (1+u)**(54-n)*(1+v)**n-1
        many = total-(54-n)*u-n*v
        assert 0 <= many <= MANY*BASE and total < BASE
        full = maximum_timing_weight(n, u, v, [1]+[0]*n)
        assert full <= 35*u+v/2 < FULL*BASE
        # A lo sumo 35 esperas anteriores a la última interacción. El
        # último color del cable es <=35+n-1. Se incluye la espera inicial.
        absent = (35*u+v/(1-discount))*discount+((54-n)*u+v/2)*discount**n
        weighted = maximum_timing_weight(n, u, v, [0]+[discount**h for h in range(1, n+1)])
        assert weighted <= absent and many+absent < DEFECT*BASE
    return (1-2*RC/3)**(-19)/(1-RM)/(1-target)**35-1


def group_blocks(first, limit, y, full=FULL, defect=DEFECT):
    """Coeficientes por número de datos y marcas, a la escala racional heredada."""
    def empty():
        return [[0]*(2*e+1) for e in range(limit+1)]

    data, checks, site, overlap, bridge, inner, check_children = [empty() for _ in range(7)]
    e_bridge, e_inner, e_check, e_line, e_rect = [empty() for _ in range(5)]
    lines, rectangles = empty(), empty()
    checks[0][0] = first.fixed(y)
    e_bridge[0][0] = first.fixed(first.prior.exp(first.DM*F(checks[0][0], first.SCALE)))
    for series in (e_inner, e_check, e_line, e_rect):
        series[0][0] = first.SCALE

    for e in range(1, limit+1):
        site[e][:-2] = [ceil(first.X*v) for v in e_bridge[e-1]]
        extra = first.mul_sum(site, e_inner, e, range(1, e))
        overlap[e] = [x+t for x, t in zip(site[e], extra)]
        inner[e] = [612*t for t in overlap[e]]
        lines[e] = [ceil(630*full*t) for t in site[e]]
        rectangles[e] = [ceil(10133*full*t) for t in site[e]]
        for source, forest in ((inner, e_inner), (lines, e_line), (rectangles, e_rect)):
            forest[e] = first.mul_sum(source, forest, e, range(1, e+1), True)
        data[e] = [z+ceil((18*y+defect)*x) for z, x in zip(extra, site[e])]
        for marks, count, forest in ((1, 613, e_line), (2, 5508, e_rect)):
            members = first.mul_sum(site, forest, e, range(1, e))
            for j, (one, more) in enumerate(zip(site[e], members)):
                if j+marks <= 2*e:
                    data[e][j+marks] += ceil(full*(one+count*more))
                else:
                    assert one == more == 0
        check_children[e] = [first.MD*v for v in data[e]]
        e_check[e] = first.mul_sum(check_children, e_check, e, range(1, e+1), True)
        checks[e] = [ceil(y*v) for v in e_check[e]]
        bridge[e] = [first.DV*a+first.DM*b for a, b in zip(data[e], checks[e])]
        e_bridge[e] = first.mul_sum(bridge, e_bridge, e, range(1, e+1), True)
    rooted = [[d+ceil(35*y*s) for d, s in zip(dr, sr)] for dr, sr in zip(checks, site)]
    return data, rooted


def two_data(first, top, n):
    debt = sum(first.DUAL)
    result = 2*MANY*first.DV*first.power(max(first.fees(2)))
    for m in (17, 18):
        def fee(k, l, after=None, common=0):
            missing = n+m-2*common-k-l
            assert missing >= 0
            return first.pair_fee(n, m, k, l, after)-debt*missing

        shared = []
        for i in top.overlaps(n, m):
            largest = max(fee(k, l, common=i) for k in range(n-i+1) for l in range(m-i+1))
            shared.append(first.power(largest)/i)
        result += top.incidence(n, m)*max(shared)
        incomplete = max(fee(k, l) for k in range(n) for l in range(m))
        big = max(fee(k, l, (j, k+l-33)) for k in range(n+1) for l in range(m+1)
                  if k+l >= 33 and (k == n or l == m) for j in range(max(k, l)+1))
        # Si la primera acción retira un dato completo, quedan <1+18 bits.
        result += first.DV*(first.power(incomplete)+first.power(big)+2**19)
        for h in (17, 18):
            by_incidence = []
            for i in top.overlaps(n, h):
                other = h-i
                if other <= 0 or (h != m and other > 1):
                    continue
                largest = max(fee(k, l, (max(k-i, l-other), k+l-h))
                              for k in range(i, n+1) for l in range(other, m+1)
                              if k == n or l == m)
                readers = 18 if m == 17 else 17
                by_incidence.append(readers*first.power(largest)/i)
            result += top.incidence(n, h)*max(by_incidence, default=F(0))
    return result


def tail(first, price):
    held = 500000
    offsets = [0, 0]
    for threshold, weight in first.FIELDS:
        share = min(held, weight)
        held -= share
        for i, n in enumerate((17, 18)):
            cap = ceil(n*(1+threshold)/2)-1
            offsets[i] += weight*cap+share*(n-cap)
    alpha = max(n*first.DUAL[0]-sum(row)+extra+first.EXCESS
                for n, row, extra in zip((17, 18), first.RANK, offsets))
    assert held == 0 and alpha == 85857226
    x = 4*price*first.power(alpha)
    mark = first.power(-first.CREDIT-8250000)
    site, tree, root, check = map(F, ("3e-7", "3.1e-7", "2.1e-9", "2.071e-11"))
    exp = first.prior.exp
    assert x*exp(first.DV*root+first.DM*check) < site and site*exp(612*tree) < tree
    line = FULL*(1+613*(exp(630*FULL*site)-1))
    rect = FULL*(1+5508*(exp(10133*FULL*site)-1))
    assert site*(exp(612*tree)-1+18*first.Y+DEFECT+mark*line+mark**2*rect) < root
    assert first.Y*exp(first.MD*root) < check
    return root/4**29, (check+35*first.Y*site)/4**29


def source_rows(first, top):
    prior = first.prior
    price = F("1.032")*BASE*prior.exp(prior.Z_CAP[0])
    y = first.Y
    data, meas = group_blocks(first, 28, y)
    y_hi = F("1e-7")
    pair_hi, _ = group_blocks(first, 2, y_hi)
    assert group_blocks(first, 4, y, F(1), F(0)) == first.block_rows(4, y)

    def row(arr, e):
        return (price/first.X)**e / first.SCALE * sum(
            count*first.power(fee) for count, fee in zip(arr[e], first.fees(e)))

    t_data, t_meas = tail(first, price)
    data_rest = sum(row(data, e) for e in range(3, 29))+t_data
    one = price*prior.exp(first.DM*y)*2**18
    for n in (17, 18):
        cost = [first.DUAL[0]*k-sum(first.RANK[n-17][:k])
                + first.recurrence((n,))(0, k, k)-sum(first.DUAL)*(n-k) for k in range(n+1)]
        assert max(cost) <= 18*first.MICRO
    roots = [one+data_rest+y/y_hi*row(pair_hi, 2)+price**2*two_data(first, top, n)
             for n in (17, 18)]
    measured = first.MD*y*one+sum(row(meas, e) for e in range(2, 29))+t_meas
    return (*roots, measured), (t_data, t_meas)


def verify_uniform(first, top):
    prior = first.prior
    first.check_credit()
    assert not any(prior.check_local(prior.GENERAL).values())
    gm = check_fibre_laws(first)
    actual, rest = source_rows(first, top)
    limits = prior.fractions("6e-24 1.3e-23 1.5e-25")
    assert all(a < b for a, b in zip(actual, limits)), [float(x) for x in actual]
    branches = prior.fractions("7.268855917008e-19 7.268859154131e-19 1.022239188961e-17 "
                "1.066266327011e-11 1.063571738677e-11 1.587554982135e-10 "
                "5.503252902255e-24 5.826965119896e-24 1.392966706244e-21 "
                "2.813300297532e-22 2.978789876782e-22 5.064592264202e-21 "
                "3.80789064072e-7 4.031878144832e-7 6.854443070565e-6")
    cert = replace(prior.GENERAL, sources=limits, branches=branches)
    mat, _ = prior.matrices()
    rv, core_check = prior.check_trees(cert, mat)
    passive = prior.XI
    rm = core_check+passive
    tilt = F(".025")
    assert gm*prior.exp(tilt+first.MD*rv+83811*rm) < passive
    assert gm*prior.exp(tilt+83811*passive) < prior.XI
    assert BASE*prior.exp(tilt+199*first.DM*passive) < F("1.032")*BASE
    load_v = 199*((first.DV+1)*rv+first.DM*rm)
    load_c = first.MD*rv+83811*rm
    assert max(load_v, load_c) < F(".01")
    print("Familia (17,18), rc <= 2e-35, rm <= 2e-11. Filas VV, CC y medida:",
          *(f"{float(v):.10g}" for v in actual))
    print("Cola desde 29 datos:", *(f"{float(v):.10g}" for v in rest))
    print("Cargas exteriores por dato y medida:", *(f"{float(v):.12g}" for v in (load_v, load_c)),
          "frente a 0.01.")


def main():
    total, norm = finite_f17()
    print("F17 original, rc <= 5e-33, rm <= 2e-11. Suma de odds:",
          f"{float(total):.12g}. Norma logarítmica en radio 2: <= {float(norm):.12g}.")
    here = Path(__file__).resolve().parent
    refs = (here/ARCHIVE, here.parent/ARCHIVE, here.parents[1]/"qec"/ARCHIVE)
    archive = next((p for p in refs if p.is_file()), None)
    if archive is None or sha256(archive.read_bytes()).hexdigest() != DIGEST:
        raise RuntimeError(f"Se requiere {ARCHIVE} del manuscrito con SHA-256 {DIGEST}.")
    with TemporaryDirectory() as td, ZipFile(archive) as source:
        tmp = Path(td)
        for name in source.namelist():
            if name.startswith("certificate/"):
                dest = tmp/Path(name).name
                dest.write_bytes(source.read(name))
        sys.path.insert(0, str(tmp))
        import verificar_H5A_primer_apagado_por_grupos as first
        import verificar_H5A_fusiones_sin_repetir as top
        verify_uniform(first, top)


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para verificar las desigualdades.")
    main()
