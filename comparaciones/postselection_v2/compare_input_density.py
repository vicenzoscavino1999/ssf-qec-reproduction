"""Check a seven-lit-check discard rule on the affine 225-qubit Bernoulli example.

Run with Python 3.12, passing the root of the reproduction package.
The extension from reduced two-data patterns to arbitrary readouts uses Lemma 29.
The sector counts and union bounds are those of Appendix E.2.1 of the manuscript.
"""

import argparse
import json
from collections import Counter
from decimal import Decimal, localcontext, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
from itertools import combinations
from math import comb
from pathlib import Path
from runpy import run_path


def decimal_bound(value, rounding):
    with localcontext() as context:
        context.prec = 12
        context.rounding = rounding
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def readout_tail(physical_weight, m, q):
    # Independently erase lit checks and light previously unlit checks.
    def binomial(n, k):
        return comb(n, k) * q**k * (1 - q)**(n - k)

    available = m - physical_weight
    return sum(binomial(physical_weight, erased)
               * sum(binomial(available, added)
                     for added in range(max(0, 7 - physical_weight + erased), available + 1))
               for erased in range(physical_weight + 1))


def compare(release):
    audit = run_path(str(release / 'comparaciones' / 'auditar_esquinas_HGP_225.py'))
    _, checks, generators = audit['incidence']()
    tags, decode, _ = audit['decoder'](checks, generators)
    geometry = json.loads((release / 'sectores/referencia/resultado_225.json').read_bytes())
    census = geometry['interval']['positive_census']
    archived = json.loads((release / 'sectores/flag_rates/pilot-results.json').read_bytes())
    n, m = len(checks), geometry['geometry']['checks']
    assert (n, m) == (225, 108) and max(x.bit_count() for x in checks) == 4

    silent_pairs, silent_triples = Counter(), Counter()
    for row in range(12):
        sites = range(12 * row, 12 * row + 12)
        for u, v in combinations(sites, 2):
            physical = checks[u] ^ checks[v]
            for bit in range(9 * row, 9 * row + 9):
                readout = 1 << bit
                syndrome = physical ^ readout
                end, correction = decode(syndrome)
                if end == readout and correction != (tags[u] ^ tags[v]):
                    silent_pairs[syndrome.bit_count()] += 1
        for u, v, w in combinations(sites, 3):
            syndrome = checks[u] ^ checks[v] ^ checks[w]
            end, correction = decode(syndrome)
            if end == 0 and correction != (tags[u] ^ tags[v] ^ tags[w]):
                silent_triples[syndrome.bit_count()] += 1
    assert silent_pairs == {7: census['2,1']} == {7: 324}
    assert silent_triples == {9: census['3,0']} == {9: 36}

    p, q = Fraction(1, 10**6), Fraction(1, 10**4)
    def weight(a, b):
        return p**a * (1-p)**(n-a) * q**b * (1-q)**(m-b)

    spectra = (Counter({0: 1}), Counter(x.bit_count() for x in checks),
               Counter((x ^ y).bit_count() for x, y in combinations(checks, 2)))
    tails = {h: readout_tail(h, m, q) for spectrum in spectra for h in spectrum}
    abort_lo = sum(p**a * (1-p)**(n-a) * sum(count * tails[h] for h, count in spec.items())
                   for a, spec in enumerate(spectra))
    abort_hi = abort_lo + 1 - sum(comb(n, a)*p**a*(1-p)**(n-a) for a in range(3))
    elementary_abort_hi = comb(m, 7)*q**7 + n*p*comb(m, 3)*q**3 + comb(n, 2)*p**2
    assert abort_lo <= abort_hi < elementary_abort_hi

    # All two-data silent inputs have seven lit checks in their faulty row.
    # Exterior flipped readouts only add lights. All 36 silent (3,0) inputs have nine.
    orbit = archived['orbit_bounds']['225']
    counts = {(3,1): orbit['3,1']['count_upper'], (3,2): orbit['3,2']['count_upper'],
              (3,3): 139_914_784_512, (4,0): archived['C40_direct225']['silent'],
              (4,1): orbit['4,1']['count_upper']}
    remaining = sum(count*weight(a, b) for (a, b), count in counts.items())
    remaining += comb(n, 3)*p**3*comb(m, 4)*q**4
    remaining += comb(n, 4)*p**4*comb(m, 2)*q**2 + comb(n, 5)*p**5
    density_risk_hi = remaining / (1 - abort_hi)
    ssf_risk_lo = census['3,0'] * weight(3, 0)  # D=0 implies the original flag keeps them.
    assert density_risk_hi < ssf_risk_lo
    print('Silent reduced (2,1) and all (3,0) by input weight:', silent_pairs, silent_triples)
    print('Discards per million:', decimal_bound(10**6 * abort_lo, ROUND_FLOOR),
          decimal_bound(10**6 * abort_hi, ROUND_CEILING))
    print('Elementary upper bound on discards per million:', decimal_bound(10**6 * elementary_abort_hi, ROUND_CEILING))
    print('Density-rule risk per retained round <=', decimal_bound(density_risk_hi, ROUND_CEILING))
    print('SSF-flag risk per retained round >=', decimal_bound(ssf_risk_lo, ROUND_FLOOR))
    print('SSF-flag risk / density-rule risk >=', decimal_bound(ssf_risk_lo / density_risk_hi, ROUND_FLOOR))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('release', type=Path)
    compare(parser.parse_args().release)
