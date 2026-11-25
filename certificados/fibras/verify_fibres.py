"""Exact examples for the one-round data/check fibre proposition. Standard Python."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import prod


def xor_product(laws):
    result = {0: Q(1)}
    for law in laws:
        updated = Counter()
        for previous, p in result.items():
            for value, q in law.items():
                updated[previous ^ value] += p*q
        result = updated
    return result


def relative(law, baseline):
    return {x: p/baseline for x, p in law.items()}


def bits(probability, effect):
    return {0: 1-probability, effect: probability}


# Two data share two Z checks. Both physical columns are 11, so a filtered
# residual may be 11 and two competing actions have exactly the same score.
GATES = ((0, 0), (1, 1), (0, 1), (1, 0))
SINGLE_PAULI_X = (0, 1, 1, 0)  # I, X, Y, Z


def laws(present, r, rm):
    # A key packs the two data bits below the two syndrome (or discrepancy) bits.
    primitive = []
    data = [[], []]
    checks = [[], []]
    baselines = []
    for k, (v, j) in enumerate(GATES):
        if not (present >> k) & 1:
            continue
        e, check = 1 << v, 1 << (j+2)
        future = sum(1 << (b+2) for a, b in GATES[k+1:] if a == v)
        past = 12 ^ future  # includes the current ideal gate, regardless of J
        direct = Counter()
        for i, h in product(range(4), repeat=2):
            a, b = SINGLE_PAULI_X[i], SINGLE_PAULI_X[h]
            direct[a*e ^ a*future ^ b*check] += r/15 if i or h else 1-r
        primitive.append(relative(direct, 1-r))
        baselines.append(1-r)
        pa, pt = 8*r/15, 4*r/(15-8*r)
        data[v].append({0: 1-pa, e ^ past: pa/2, e ^ past ^ check: pa/2})
        checks[j].append(bits(pt, check))

    if present & (1 << 4):  # final data wait on wire zero
        p = 2*r/3
        primitive.append(relative(bits(p, 1), 1-r))
        data[0].append(bits(p, 1 ^ 12))
        baselines.append(1-r)
    for j in range(2):  # two joint Pauli/readout locations
        if present & (1 << (5+j)):
            flip = 2*r/3 + rm - 2*(2*r/3)*rm
            check_law = bits(flip, 1 << (j+2))
            baseline = (1-r)*(1-rm)
            primitive.append(relative(check_law, baseline))
            checks[j].append(check_law)
            baselines.append(baseline)
    fibres = [xor_product(wire) for wire in (*data, *checks)]
    scale = prod(law[0] for law in fibres)/prod(baselines)
    effective = xor_product(relative(law, law[0]) for law in fibres)
    # (E,D) -> (E,s); H has two identical columns of weight two.
    joint = {key ^ (12 if (key & 3).bit_count() % 2 else 0): value
             for key, value in effective.items()}
    return xor_product(primitive), joint, scale


def partitions(joint, allowed):
    result = [Q(0)]*4
    for key, mass in joint.items():
        error, syndrome = key & 3, key >> 2
        # A strict descent on column 11 is possible only from syndrome 11.
        # The full competition then picks wire zero under the fixed tie-break.
        correction = int(syndrome == 3)
        residual = error ^ correction
        if correction & ~allowed or residual.bit_count() % 2:
            continue
        for u in range(4):
            result[u] += (-1)**((u & residual).bit_count()) * mass
    return result


def check_subhosts():
    detects_sign = False
    for present in range(1 << 7):
        primitive, effective, scale = laws(present, Q(1, 20), Q(1, 17))
        assert primitive == {x: scale*p for x, p in effective.items()}
        for allowed in range(4):
            direct = partitions(primitive, allowed)
            assert direct == [scale*p for p in partitions(effective, allowed)]
            detects_sign |= direct[0] != direct[1]
    assert detects_sign


def check_rates():
    r, rm, p0 = map(Q, ("4.2e-8", ".001", "3.2e-7"))
    for h in (3, 4):
        data = (1-2*r/3)**-(12-h) * (1-8*r/15)**-h - 1
        assert data < Q("3.192001e-7") < p0/(1-2*p0)
    pt = 4*r/(15-8*r)
    check = (1-2*r/3)**-5/(1-rm)/(1-pt)**7 - 1
    assert check < Q(".001001220") < rm/(1-2*rm)

    r, rm = map(Q, ("2e-7", "1e-5"))
    pt = 4*r/(15-8*r)
    parity = (1-4*r/3)**5*(1-2*rm)*(1-2*pt)**7
    odds = (1-parity)/(1+parity)
    assert odds > Q("1.1040e-5") > Q("1.00003e-5") > rm/(1-2*rm)
    print(f"Corner 3 full-check odds {float(odds):.12g}, budget {float(rm/(1-2*rm)):.12g}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run without -O to enable the rational checks.")
    check_subhosts()
    check_rates()
    print("Joint response and signed normalizations agree on all 512 primitive/action subhosts.")
