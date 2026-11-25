"""Cadenas de fusiones y primer apagado. Certificado aritmético H5-A."""

from dataclasses import replace
from fractions import Fraction as F
from functools import cache
from math import ceil, comb, factorial
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verificar_H5A_encuentro_de_dos_frentes as fronts
import verificar_H5A_primer_apagado_por_grupos as first

pools = first.prior
P = F('5e-33')
Q = pools.Q
PRICE = F('1.032') * P * pools.exp(pools.Z_CAP[0])
Y = first.Y
J = 83810
WIDTHS = (17, 18)


def incidence(n, m):
    readers = 18 if m == 17 else 17
    return n * (readers - int(n == m))


def overlaps(n, m):
    return range(1, n + 1 if n == m else 2)


@cache
def profile(widths, lights):
    entry = sum(first.DUAL[0]*k - sum(first.RANK[n-17][:k])
                for n, k in zip(widths, lights, strict=True))
    future = first.recurrence(tuple(sorted(widths)))
    return entry + future(0, max(lights), sum(lights))


@cache
def weight(fee):
    return first.power(fee)


def pair(n):
    total = F(0)
    for m in WIDTHS:
        best = max(weight(profile((n, m), (n-i, m-i))) / i for i in overlaps(n, m))
        total += incidence(n, m) * best
        left = max(0, n+m-33)
        best_fee = max(first.pair_fee(n, m, n, m, (k, left)) for k in range(left+1))
        total += first.DV * weight(best_fee)
        for h in WIDTHS:
            choices = []
            for i in overlaps(n, h):
                j = h-i
                if j <= 0 or (j > 1 and m != h):
                    continue
                fee = first.pair_fee(n, m, n, m, (max(n-i, m-j), n+m-h))
                readers = 18 if m == 17 else 17
                choices.append(readers * weight(fee-first.CREDIT+2*first.EXCESS) / i)
            total += incidence(n, h) * max(choices, default=F(0))
    return total * PRICE**2


def four_unmarked(n, matrix):
    connected = F(4**2, factorial(3)) * 612**3 * weight(first.fees(4)[0])
    total = connected
    for m in WIDTHS:
        for v in WIDTHS:
            bridges = matrix[n-17][v-17] + matrix[m-17][v-17]
            for w in WIDTHS:
                widths = (n, m, v, w)
                largest = max(weight(profile(widths, (n-i, m-i, v-j, w-j))) / (i*j)
                              for i in overlaps(n, m) for j in overlaps(v, w))
                total += incidence(n, m) * bridges * incidence(v, w) * largest
    return total * PRICE**4


@cache
def cover_fee(e, r):
    """Saldo de e datos una vez consumida una línea cubierta por r de ellos."""
    levels = [(ceil(18*(1+u)/2)-1, ceil(17*(1+u)/2), w) for u, w in first.FIELDS]

    def capacity(k):
        return (e-r)*k + min(r*k, 18*r-17)

    @cache
    def visit(index, k, balance):
        if index == len(levels):
            return 0
        cap, minimum, multiplier = levels[index]
        best = 0
        for j in range(min(cap, k)+1):
            spent = max(minimum, ceil(F(17*k, 18))) if j < k else 0
            remaining = min(capacity(j), max(0, balance-spent))
            best = max(best, multiplier*remaining + visit(index+1, j, remaining))
        return best

    values = []
    for k in (17, 18):
        entry = e*max(first.DUAL[0]*min(n, k)-sum(row[:k])
                      for n, row in zip(WIDTHS, first.RANK, strict=True))
        for j in range(min(k, 17)+1):
            values.append(entry + visit(0, j, min(capacity(j), e*k-17)))
    return max(values) - max(0, first.CREDIT-first.EXCESS*e)


def multiple_cover(e):
    total = F(0)
    for r in range(2, min(e, 18)+1):
        m = e-r
        b = m//2
        layouts = (613 * 35**(r-1) * F(e**m, factorial(m)) *
                   comb(m, b) * 612**(m-b) * first.DV**b)
        total += layouts * weight(cover_fee(e, r))
    return PRICE**e * total


def tail():
    # Misma inclinación de la serie infinita que en primer_apagado, a la nueva p.
    reserved = 500000
    offsets = [0, 0]
    for u, w in first.FIELDS:
        held = min(reserved, w)
        reserved -= held
        for i, n in enumerate(WIDTHS):
            cap = ceil(n*(1+u)/2)-1
            offsets[i] += w*cap + held*(n-cap)
    assert reserved == 0
    alpha = max(first.DUAL[0]*n-sum(r)+o+first.EXCESS
                for n, r, o in zip(WIDTHS, first.RANK, offsets, strict=True))
    assert alpha == 85857226
    x = 4*PRICE*weight(alpha)
    mark = weight(-first.CREDIT-8250000)
    a, t, d, site = map(F, ('2.1e-9', '1.80e-6', '2.071e-11', '1.78e-6'))
    exp = pools.exp
    assert x*exp(first.DV*a+first.DM*d) < site and site*exp(612*t) < t
    one = 1 + 613*(exp(630*site)-1)
    two = 1 + 5508*(exp(10133*site)-1)
    assert site*(exp(612*t)-1+18*Y+mark*one+mark**2*two) < a
    assert Y*exp(first.MD*a) < d
    return a/4**(first.N+1), (d+35*Y*site)/4**(first.N+1)


def source_rows():
    data, meas = first.block_rows(first.N, Y)
    tilt = F('1e-7')
    high, _ = first.block_rows(12, tilt)
    pure, _ = first.block_rows(12, F(0))
    matrix, neighbours = pools.matrices()
    roots = []
    measured = F(0)
    for i, n in enumerate(WIDTHS):
        factor = PRICE*pools.exp(neighbours[i][2]*Y)
        ranks = [weight(profile((n,), (k,))) for k in range(n+1)]
        roots.append(factor*sum(ranks[n-m]*comb(n, m)*Y**m for m in range(n+1)))
        measured += 35*factor*sum(ranks[n-m]*comb(n-1, m-1)*Y**m for m in range(1, n+1))
    measured += first.MD*Y*sum(roots)
    for e in range(2, first.N+1):
        fees = first.fees(e)
        unit = (PRICE/first.X)**e / first.SCALE
        rates = [unit*weight(c) for c in fees]
        sums = [v*w for v, w in zip(data[e], rates, strict=True)]
        measured += sum(v*w for v, w in zip(meas[e], rates, strict=True))
        if e == 2:
            mixed = Y/tilt * sum(v*w for v, w in zip(high[e], rates, strict=True))
            for i, n in enumerate(WIDTHS):
                roots[i] += min(sum(sums), mixed+pair(n))
            continue
        if e <= 12:
            mixed = Y/tilt * high[e][1] * rates[1]
            one = pure[e][1] * unit * weight(max(first.fees(e-1))+first.MICRO)
            sums[1] = min(sums[1], mixed+one+multiple_cover(e))
        for i, n in enumerate(WIDTHS):
            unmarked = sums[0]
            if e == 4:
                mixed = Y/tilt * high[e][0] * rates[0]
                unmarked = min(unmarked, mixed+four_unmarked(n, matrix))
            roots[i] += unmarked + sum(sums[1:])
    rest_data, rest_meas = tail()
    return tuple(v+rest_data for v in roots) + (measured+rest_meas,)


def general():
    first.check_credit()
    assert not any(pools.check_local(pools.GENERAL).values())
    actual = source_rows()
    print('Filas fuente generales:', *(f'{float(v):.10g}' for v in actual), flush=True)
    sources = pools.fractions('6e-24 1.3e-23 1.5e-25')
    assert all(v < limit for v, limit in zip(actual, sources, strict=True))
    branches = pools.fractions('7.268855917008e-19 7.268859154131e-19 1.022239188961e-17 '
                               '1.066266327011e-11 1.063571738677e-11 1.587554982135e-10 '
                               '5.503252902255e-24 5.826965119896e-24 1.392966706244e-21 '
                               '2.813300297532e-22 2.978789876782e-22 5.064592264202e-21 '
                               '3.80789064072e-7 4.031878144832e-7 6.854443070565e-6')
    cert = replace(pools.GENERAL, rate=P, sources=sources, branches=branches)
    matrix, _ = pools.matrices()
    core, measured = pools.check_trees(cert, matrix)
    rv, rm = pools.check_exterior(cert, core, measured, F('1.032')*P)
    ref = F('1.907e-6')*pools.odds(P)/pools.odds(F('1e-6'))
    residual = rv+ref, rm-pools.odds(Q)*(1+F('.015'))+ref
    assert residual[0] < F('2e-16') and residual[1] < F('2.6e-13')
    print('PASS general p <= 5e-33, q <= 2e-11. Filas del resto en la esquina:',
          *(f'{float(v):.10g}' for v in residual), flush=True)


def distance_rows():
    sp, sq = (pools.odds(r)*pools.exp(F('.015')) for r in (F('1e-30'), Q))
    growth, base = 273**2, 2**20
    sizes = fronts.spatial.ball_sizes
    offsets = [v-272*u for u, v in zip(sizes(1), sizes(2), strict=True)]
    assert offsets == [18, 17, 17, 18] and all(o <= u for o, u in zip(offsets, sizes(1)))
    assert growth < base
    kernel = {r: F(fronts.data_from_check(r), 1)/(1-F(growth, base)) for r in (3, 5, 6)}
    y, passive = F('1.2e-29'), F('1.634e-10')
    bd = 135*first.DM*(passive+kernel[3]*y)
    bm = J*passive+kernel[5]*y
    # Cada check fresco fuera de las huellas es una luz creada o una medida de D.
    fresh = 18 + 99*18//2
    reads = 2*fresh*J+18
    assert fresh == 909 and reads == 152366598
    assert 8*sq*pools.exp(2*bm) < passive
    assert 8*sq*pools.exp(kernel[5]*y+2*J*passive) < passive
    closed = 4*sp*pools.exp(2*bd+fresh*kernel[5]*y+reads*passive)
    assert closed < y
    measured = passive*(1+kernel[6]*y)
    assert measured < F('1.32e-8')
    print('Fusiones con 20 saltos de exceso <= 2e + 3d. Consultas frescas <=',
          f'{reads}e + {2*J}d. Filas <= 1.2e-29 por dato, < 1.32e-8 por medida.')
    return y, measured


def target():
    rows = distance_rows()
    prev = fronts.spatial.previous
    stored = ((F('3.68e-11'), F('2.04e-11')), (F('1.1e-29'), F('2.04e-11')),
              (F('2e-16'), F('2.057e-11')), (F('3.53e-11'), F('2.05e-11')))
    all_rows = (prev.two_data(), prev.activity_counts(), *stored, fronts.rows(), rows)
    reference = F('1.907e-6')*pools.odds(F('1e-30'))/pools.odds(F('1e-6'))
    totals = [sum(row)+reference for row in zip(*all_rows, strict=True)]
    assert totals[0] < F('1.94e-5') and totals[1] < F('.00070')
    assert max(totals)+F('.004298') < F('.004998')
    print('PASS parcial p = 1e-30. Filas acumuladas:',
          *(f'{float(v):.10g}' for v in totals), 'Resto conjunto abierto.')


def compare_chain(g):
    """Cierres exactos del ejemplo de g grupos sobre su partición por capas."""
    degrees = {}

    def edge(a, b, forward, backward):
        degrees.setdefault(a, {})[b] = forward
        degrees.setdefault(b, {})[a] = backward

    def trunk(start, name):
        for j in range(3):
            dest = (*name, j)
            edge(start, dest, 8, 9)
            start = dest
        return start

    seeds = [('fuente', i) for i in range(g)]
    ends = []
    for i, seed in enumerate(seeds):
        edge(seed, seed, 9, 9)
        ends.append(trunk(seed, ('brazo', i)))
    previous = ends[0]
    gates = []
    for i in range(1, g):
        gate = ('unión', i)
        gates.append(gate)
        if i < g-1:
            edge(previous, gate, 8, 5)
            edge(ends[i], gate, 8, 4)
            previous = trunk(gate, ('tramo', i))
        else:
            edge(previous, gate, 8, 8)
            edge(ends[i], gate, 8, 8)
            edge(gate, gate, 1, 1)
    assert all(sum(row.values()) == 17 for row in degrees.values())
    sizes = {seeds[0]: F(1)}
    pending = [seeds[0]]
    while pending:
        a = pending.pop()
        for b, count in degrees[a].items():
            size = sizes[a]*count/degrees[b][a]
            if b in sizes:
                assert sizes[b] == size
            else:
                sizes[b] = size
                pending.append(b)
    assert sizes.keys() == degrees.keys() and min(sizes.values()) > 0

    def closure(initial):
        reached = set(initial)
        while True:
            added = {a for a, row in degrees.items() if a not in reached
                     and sum(n for b, n in row.items() if b in reached) >= 9}
            if not added:
                return reached
            reached.update(added)

    for count in range(1, g+1):
        reached = closure(seeds[:count])
        assert set(gates) & reached == set(gates[:count-1])
    assert reached == set(degrees)
    # Por monotonía, excluir cada fuente de una en una domina los subconjuntos propios.
    for missing in seeds:
        reached = closure(set(seeds)-{missing})
        assert any(b not in reached for a in reached for b in degrees[a])
    return len(gates)


if __name__ == '__main__':
    general()
    target()
    for g in (2, 3, 8):
        assert compare_chain(g) == g-1
    print('PASS comparación física de los prefijos en cadenas de 2, 3 y 8 grupos.')
