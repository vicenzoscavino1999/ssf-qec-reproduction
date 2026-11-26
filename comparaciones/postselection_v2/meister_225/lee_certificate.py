"""Finite Lee acceptance test and a Bernoulli upper bound at the archived cutoff."""

from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import comb
from operator import xor

from extra_silent import independent_decoder
from inputs import BITS, LEE_CUT, M, N, P, sector_weight
from logical_gap import set_bits


def lee_silent_bounds(compare, seed_rows, checks, tags, ssf, features):
    # Log-odds ratio t=w_p/w_q is between 3/2 and 8/5. The total weight
    # in units of w_q is < 468. Any three reads, two data, or a mixture
    # of at least one of each has squared cluster norm >=3*w_q^2.
    assert 999999**2 > 9999**3 and 999999**5 < 9999**8
    assert N*F(8,5)+M == 468 and 3*10**12 > (468*LEE_CUT)**2
    simple = {0,*BITS,*checks}
    simple.update(a^b for a,b in combinations(BITS,2))
    assert len(simple) == 6112
    assert all(lee or word in simple for word,(_,lee) in features.items())
    # Order-zero LSD assembles a syndrome-consistent solution supported on
    # its final valid clusters. Every accepted input is therefore in simple.
    independent_tags, independent_ssf = independent_decoder(compare, checks)
    for label, decode in ((tags,ssf),(independent_tags,independent_ssf)):
        singles = {(0,0),*zip(checks, label, strict=True)}
        for s in simple:
            stop, correction_label = decode(s)
            assert (s^stop,correction_label) in singles

    # At silent failure the effective correction can thus be represented
    # by at most one site. With five data errors its completion has weight
    # exactly six. Every weight-six logical is a single seed-code row.
    minimum_words = [w for w in range(1 << 12) if w.bit_count()==6
                     and all((w & row).bit_count()%2 == 0 for row in seed_rows)]
    logicals = [w << (12*r) for r in range(12) for w in minimum_words]
    fives, proven = set(), set()
    for logical in logicals:
        indexes = tuple(set_bits(logical))
        assert reduce(xor,(checks[i] for i in indexes)) == 0
        assert reduce(xor,(tags[i] for i in indexes)) != 0
        for v in indexes:
            e = logical ^ (1 << v)
            fives.add(e)
            s = checks[v]  # observed syndrome when the readout is perfect
            stop, correction = ssf(s)
            if not features[s][1] and stop == 0 and correction != reduce(
                    xor,(tags[j] for j in set_bits(e)),0):
                proven.add(e)
    assert len(logicals)==72 and len(fives)==432 and proven==fives
    five = len(fives)*P**5*(1-P)**(N-5)
    tail6 = 1-sum(comb(N,i)*P**i*(1-P)**(N-i) for i in range(6))
    # The first bound uses only d_X=6. The refined bound also classifies
    # all possible weight-five supports and pays the entire six-data tail.
    crude_upper = comb(N,5)*P**5
    assert five+tail6 <= crude_upper
    return len(simple), len(logicals), len(fives), len(proven), (
        len(proven)*sector_weight(5,0), five+tail6, crude_upper)
