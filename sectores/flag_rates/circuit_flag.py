"""Deterministic circuit postselection bounds for the eight affine 225-qubit points."""
import argparse
from array import array
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction as F
from itertools import combinations
import json
from math import comb, lcm, prod
from pathlib import Path
from runpy import run_path


def outward(value, upper=False):
    with localcontext() as ctx:
        ctx.prec = 13
        ctx.rounding = ROUND_CEILING if upper else ROUND_FLOOR
        return format(Decimal(value.numerator) / Decimal(value.denominator), ".12E")


def subset_masks(support):
    subset = support
    while subset:
        yield subset
        subset = (subset - 1) & support


def flag_catalogue(checks, run, audit_run):
    words = sorted({word for support in checks for word in subset_masks(support)})
    index = {s: j for j, s in enumerate(words)}

    def flagged(s):
        end = run(s)[0]
        if audit_run(s)[0] != end:
            raise ArithmeticError(f"Different stopping syndromes for input {s}")
        return 0 < end.bit_count() < s.bit_count()

    singles = [flagged(s) for s in words]
    neighbours = [array("H") for s in words]
    for i, s in enumerate(words):
        for j in range(i + 1, len(words)):
            if flagged(s ^ words[j]):
                neighbours[i].append(j)
                neighbours[j].append(i)
    bits = [1 << j for j in range(108)]
    assert not flagged(0) and not any(singles[index[b]] for b in bits)
    assert sum(flagged(a ^ b) for a, b in combinations(bits, 2)) == 432
    triples = sum(flagged(a ^ b ^ c) for a, b, c in combinations(bits, 3))
    return index, singles, neighbours, triples


def silent_triples(checks, tags, run, audit_tags, audit_run, expected):
    result = []
    for row in range(12):
        for a, b, c in combinations(range(12*row, 12*row+12), 3):
            s = checks[a] ^ checks[b] ^ checks[c]
            end, correction = run(s)
            end2, correction2 = audit_run(s)
            silent = end == 0 and correction != (tags[a] ^ tags[b] ^ tags[c])
            assert silent == (end2 == 0 and correction2 != (audit_tags[a] ^ audit_tags[b] ^ audit_tags[c]))
            if silent:
                assert not (checks[a] & checks[b] or checks[a] & checks[c] or checks[b] & checks[c])
                result.append((a, b, c))
    # The archived census covers every physical triple, including triples outside one row.
    assert len(result) == expected == 36
    assert [sum(s[0]//12 == row for s in result) for row in range(12)] == [3]*12
    return result


def discard_bounds(checks, laws, norm, h, catalogue):
    index, singles, neighbours, triple_flags = catalogue
    m = 108
    # Each entry is a distribution on OBSERVED syndrome increments, including invisible faults.
    fibres = [{s: zero[s] + one[s ^ support] for s in zero}
              for support, (zero, one) in zip(checks, laws)]
    assert all(sum(row.values()) == norm for row in fibres)
    no_syndrome = F(prod(row[0] for row in fibres), norm**len(checks)) * (1-h)**m
    # Independent count of the covered mass by the number of active data wires.
    data0, data1, data2 = 1, 0, 0
    for row in fibres:
        empty, active = row[0], norm-row[0]
        data0, data1, data2 = data0*empty, data1*empty+data0*active, data2*empty+data1*active
    fibres += [{0: h.denominator-h.numerator, 1 << j: h.numerator} for j in range(m)]
    # Integer odds with one common denominator keep the pair census exact.
    scale = lcm(*(row[0] for row in fibres))
    prior = [0] * len(index)
    flagged_one = flagged_two = sum_one = sum_two = 0
    for row in fibres:
        multiplier = scale // row[0]
        weights = [(index[s], mass*multiplier) for s, mass in row.items() if s]
        for i, weight in weights:
            if singles[i]:
                flagged_one += weight
            flagged_two += weight * sum(prior[j] for j in neighbours[i])
        # Update only after the entire fibre, so no pair uses one fibre twice.
        for i, weight in weights:
            prior[i] += weight
        total = sum(weight for _, weight in weights)
        sum_two += sum_one * total
        sum_one += total
    assert 0 <= flagged_one <= sum_one and 0 <= flagged_two <= sum_two
    check_odds = h/(1-h)
    one_two = F(flagged_one*scale+flagged_two, scale**2)
    all_one_two = F(sum_one*scale+sum_two, scale**2)
    lower = no_syndrome * (one_two + triple_flags*check_odds**3)
    remainder = 1 - no_syndrome * (1 + all_one_two + comb(m, 3)*check_odds**3)
    readouts = [comb(m, k)*h**k*(1-h)**(m-k) for k in range(4)]
    covered = (data0*sum(readouts) + data1*sum(readouts[:2]) + data2*readouts[0]) / norm**len(checks)
    assert remainder == 1-covered
    assert 0 < lower < lower+remainder < 1
    return lower, lower+remainder


def unflagged_bounds(laws, norm, h, triples, four_count):
    m = 108
    n = len(laws)
    zero = [z[0] for z, _ in laws]
    data = [F(sum(one.values()), norm) for _, one in laws]
    background = [F(sum(z.values())-z[0], norm) for z, _ in laws]
    burden = sum(background)
    odds = h/(1-h)
    # On each of the 36 supports the three readout neighbourhoods are disjoint.
    # The independent check bits match their combined offset with probability
    # (1-h)^m times odds raised to the sum of the three offset weights.
    ratios = [sum(mass*odds**d.bit_count() for d, mass in one.items())/zero[i]
              for i, (_, one) in enumerate(laws)]
    clean = F(prod(zero), norm**n) * (1-h)**m * sum(prod(ratios[v] for v in S) for S in triples)
    outside = sum(prod(data[v] for v in S) * (burden - sum(background[v] for v in S))
                  for S in triples)
    pmax = max(data)
    upper = clean + outside + four_count*pmax**4 + comb(n, 5)*pmax**5
    assert 0 < clean <= upper
    return clean, upper


def compute(repo):
    comp = repo/"comparaciones"
    old = json.loads((comp/"piloto_HGP_225_esquinas_20261124.json").read_bytes())
    total_four, four = map(int, (repo/"sectores/flag_rates/counts/four_1.txt").read_text().split())
    assert (total_four, four) == (comb(225, 4), 14904)
    circ = run_path(str(comp/"circuito_HGP_225.py"))
    audit = run_path(str(comp/"auditar_esquinas_HGP_225.py"))
    _, checks, tags, _, _, run = circ["pilot"]["build"]()
    _, check2, generators = audit["incidence"]()
    assert checks == check2
    tags2, run2, _ = audit["decoder"](check2, generators)
    print("Decoding projected fibre labels and pure readout triples with both SSF implementations.", flush=True)
    catalogue = flag_catalogue(checks, run, run2)
    triples = silent_triples(checks, tags, run, tags2, run2, old["trios_Bernoulli_cero_lecturas"])
    index, singles, neighbours, flagged_readout_triples = catalogue
    result = {
        "code": "[[225,9,6]]", "observed_masks": len(index),
        "nonzero_fibre_labels": sum(2**s.bit_count()-1 for s in checks)+108,
        "distinct_observed_masks_flagged": sum(singles),
        "distinct_mask_pairs_tested": comb(len(index), 2),
        "distinct_mask_pairs_flagged": sum(map(len, neighbours))//2,
        "pure_readout_triples": comb(108, 3),
        "pure_readout_triples_flagged": flagged_readout_triples,
        "silent_data_triples_without_offset": len(triples), "rows": [],
    }
    for row in old["circuito"]:
        cc_first = row["orden_ejecutado"] == "CC primero"
        rc, rm = F(row["rc"]), F(row["rm"])
        schedule, _ = circ["schedule"](checks, cc_first=cc_first)
        laws, norm = circ["joint_fibres"](schedule, rc)
        h = circ["ancilla_rate"](rc, rm)
        gates, _ = audit["greedy_calendar"](check2, cc_first)
        other, denom, h2 = audit["circuit_laws"](gates, rc, rm)
        assert len(laws) == len(other) == len(checks) == 225
        assert norm == denom and h == h2 and all(tuple(b) == a for a, b in zip(laws, other))
        a_low, a_high = discard_bounds(checks, laws, norm, h, catalogue)
        u_low, u_high = unflagged_bounds(laws, norm, h, triples, four)
        risk_low, risk_high = u_low/(1-a_low), u_high/(1-a_high)
        lower_unfiltered = F(row["determinista"]["L"])
        gain = lower_unfiltered/risk_high
        assert gain > 1
        result["rows"].append({
            "order": "CC first" if cc_first else "VV first", "rc": row["rc"], "rm": row["rm"],
            "unfiltered_deterministic_lower": row["determinista"]["L"],
            "discards_per_million": [outward(10**6*a_low), outward(10**6*a_high, True)],
            "unflagged_joint": [outward(u_low), outward(u_high, True)],
            "risk_per_retained_round": [outward(risk_low), outward(risk_high, True)],
            "gain_lower": outward(gain),
        })
    assert len(result["rows"]) == 8
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repo", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = compute(args.repo.resolve())
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(f"Verified {len(result['rows'])} deterministic circuit intervals.", flush=True)
