"""Postselection decisions on the same one-round affine [[225,9,6]] SSF experiment."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import json
from math import ceil, comb, lcm, prod
from pathlib import Path
from runpy import run_path

from cluster_score import check_circuit_priors, detector_priors, scorer

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
COMP = REPO / "comparaciones"
CIRC = run_path(str(COMP / "circuito_HGP_225.py"))
FLAG = REPO / "sectores" / "flag_rates"
N, M = 225, 108
BITS = tuple(1 << j for j in range(M))
GRID = 10**6


def load(path):
    return json.loads(path.read_bytes())


def endpoints(lo, hi):
    assert 0 <= lo <= hi
    with localcontext() as ctx:
        ctx.prec = 12
        def rounded(x, direction):
            ctx.rounding = direction
            return format(Decimal(x.numerator)/Decimal(x.denominator), "E")
        return [rounded(lo, ROUND_FLOOR), rounded(hi, ROUND_CEILING)]


def projected_wires(checks, laws):
    wires = []
    for support, (zero, one) in zip(checks, laws, strict=True):
        row = Counter(zero)
        for offset, mass in one.items():
            row[support ^ offset] += mass
        wires.append(row)
    return wires


def wire_odds(wires, norm):
    scale = lcm(*(row[0] for row in wires))
    no_wire = F(prod(row[0] for row in wires), norm**len(wires))
    odds = [[(s, value*(scale//row[0])) for s, value in row.items() if s and value]
            for row in wires]
    return no_wire, scale, odds


def density_discards(h, no_wire, scale, odds):
    # Condition on zero, one or two wire increments. The 108 check bits remain free.
    prior = Counter()
    ones = Counter()
    twos = Counter()
    for wire in odds:
        for word, weight in wire:
            previous = Counter()
            for other, mass in prior.items():
                previous[(word ^ other).bit_count()] += mass
            for size, mass in previous.items():
                twos[size] += weight*mass
            ones[word.bit_count()] += weight
        prior.update(dict(wire))

    @lru_cache(None)
    def lit_tail(size):
        result = 0
        for erased in range(size+1):
            need = max(0, 7-size+erased)
            new_lights = 1-sum(comb(M-size, j)*h**j*(1-h)**(M-size-j) for j in range(need))
            result += comb(size, erased)*h**erased*(1-h)**(size-erased)*new_lights
        return result

    low = no_wire*(lit_tail(0) +
                   sum(mass*lit_tail(s) for s, mass in ones.items())/scale +
                   sum(mass*lit_tail(s) for s, mass in twos.items())/scale**2)
    covered = no_wire*(1+F(sum(ones.values()), scale)+F(sum(twos.values()), scale**2))
    return low, low+1-covered


def cluster_cut(score, h, no_wire, scale, odds, budget, include_triples):
    # Calibrate solely from input probabilities; SSF outcomes play no part in the cutoff.
    t = h/(1-h)
    last = 3 if include_triples else 2
    unit = lcm(scale, t.denominator**last)
    no_input = no_wire*(1-h)**M
    word_mass = Counter()
    for wire in odds:
        for word, mass in wire:
            word_mass[word] += mass*(unit//scale)
    for size in range(1, last+1):
        weight = int(unit*t**size)
        assert weight == unit*t**size
        for group in combinations(BITS, size):
            word_mass[sum(group)] += weight
    hist = Counter()
    for word, mass in word_mass.items():
        hist[ceil(GRID*score(word))] += mass
    total = sum(hist.values())
    rest = 1-no_input*(1+F(total, unit))
    assert 0 <= rest < budget
    allowance = (budget-rest)*unit/no_input
    rejected = total
    cutoff = 0
    for cutoff, mass in sorted(hist.items()):
        rejected -= mass
        if rejected <= allowance:
            break
    low = no_input*F(rejected, unit)
    assert low+rest <= budget
    return cutoff, (low, low+rest)


def refine_discards(rejects, h, no_wire, scale, odds, bounds, bern):
    # Additional disjoint classes at the two headline points, after fixing the cutoff.
    t = h/(1-h)
    z = no_wire*(1-h)**M
    weight = Counter()
    for row in odds:
        weight.update(dict(row))
    two = sum(mass*sum(rejects(word ^ bit) for bit in BITS) for word,mass in weight.items())
    three = sum(rejects(a ^ b ^ c) for a,b,c in combinations(BITS,3))
    accepted_mass = F(sum(weight.values())*M-two,scale)*t + (comb(M,3)-three)*t**3
    rejected_mass = F(two,scale)*t+three*t**3
    if bern:
        assert all(len(row)==1 for row in odds)
        positives = [row[0] for row in odds]
        total = marked = 0
        for (a,wa),(b,wb) in combinations(positives,2):
            total += wa*wb
            if rejects(a ^ b):
                marked += wa*wb
        accepted_mass += F(total-marked,scale**2)
        rejected_mass += F(marked,scale**2)
    lower,upper = bounds
    tightened = lower+z*rejected_mass, upper-z*accepted_mass
    assert lower<=tightened[0]<=tightened[1]<=upper
    return tightened


def geometry():
    H, checks, tags, _, _, decode = CIRC["pilot"]["build"]()
    pairs = []
    assert CIRC["pilot"]["row_sectors"](checks, tags, decode, patterns=pairs) == 324
    triples, _ = CIRC["physical_triples"](H, checks, tags, decode)
    assert len(triples) == 1656
    weights = Counter((checks[u]^checks[v]^checks[w]^d).bit_count() for (u,v,w),d in triples)
    assert weights == {4:324, 8:324, 9:1008} and ((0,3,4),1 << 3) in triples
    return checks, pairs, triples


def pair_decisions(pairs, checks, reject):
    result = {}
    for row, u, v, d in pairs:
        observed = checks[u] ^ checks[v] ^ d
        assert observed.bit_count() == 7
        result[row,u,v,d] = (int(reject(observed)),
                       sum(reject(observed ^ bit) for j, bit in enumerate(BITS) if j//9 != row))
    return result


def circuit_risks(laws, norm, h, checks, pairs, triples, decisions, reject, row):
    restricted = {}
    spectra = CIRC["pair_spectra"](laws, pairs, norm, restricted_patterns=restricted)
    low_two, upper_two = CIRC["pair_bounds"](spectra, h)
    grid = CIRC["GRID"]
    assert set(restricted) == set(decisions)
    assert [sum(v[j] for v in restricted.values()) for j in range(10)] == spectra[0]
    mass = [F((1-h)**99), h*(1-h)**98]
    reject_two = keep_two = 0
    for key, moments in restricted.items():
        inner = sum(F(c, 512*grid)*(1-2*h)**j for j, c in enumerate(moments))
        # At most 512/grid absolute error per pattern, as in the released pair inversion.
        inner = max(F(0), inner-F(512,grid))
        single, exterior = decisions[key]
        reject_two += inner*(single*mass[0]+exterior*mass[1])
        keep_two += inner*((1-single)*mass[0]+(99-exterior)*mass[1])
    assert abs(keep_two+reject_two-low_two) < F(1,10**32)
    kept = {"density": [], "cluster": []}
    for sites, d in triples:
        u,v,w = sites
        s = checks[u] ^ checks[v] ^ checks[w] ^ d
        if s.bit_count() < 7:
            kept["density"].append((sites,d))
        if not reject(s):
            kept["cluster"].append((sites,d))
    clean = {name: CIRC["clean_triple_mass"](laws, norm, keep)*(1-h)**M
             for name, keep in kept.items()}
    model = CIRC["trial_inputs"](laws,norm)
    burden = sum(model["prob_bg"]) + M*h
    m3, m4 = model["mass"][3:5]
    det = model["tail5"] + m3*min(1,burden)+m4
    b3,b4 = (F(row["estratos"][str(k)]["limite_exitos"],
               row["estratos"][str(k)]["ensayos"]) for k in (3,4))
    stat = min(det, model["tail5"]+m3*burden*b3+m4*b4)
    return {"density": (clean["density"], clean["density"]+det, clean["density"]+stat),
            "cluster": (keep_two+clean["cluster"], upper_two-reject_two+clean["cluster"]+det,
                        upper_two-reject_two+clean["cluster"]+stat)}


def bernoulli_risks(p, q, checks, triples, decisions, reject):
    archive = load(FLAG/"pilot-results.json")
    orbit = archive["orbit_bounds"]["225"]
    def w(a,b):
        return p**a*(1-p)**(N-a)*q**b*(1-q)**(M-b)
    remainder = sum(count*w(a,b) for a,b,count in (
        (3,1,orbit["3,1"]["count_upper"]), (3,2,orbit["3,2"]["count_upper"]),
        (3,3,139914784512), (4,0,archive["C40_direct225"]["silent"]),
        (4,1,orbit["4,1"]["count_upper"])))
    remainder += (comb(N,3)*p**3*comb(M,4)*q**4 +
                  comb(N,4)*p**4*comb(M,2)*q**2 + comb(N,5)*p**5)
    rejected1,rejected2 = map(sum, zip(*decisions.values()))
    keep_pairs = (324-rejected1)*w(2,1)+(324*99-rejected2)*w(2,2)
    tail_pairs = 324*p**2*(1-p)**223*q*(1-q)**8-324*(w(2,1)+99*w(2,2))
    triple_words = [checks[u]^checks[v]^checks[w_] for (u,v,w_),d in triples if d==0]
    assert len(triple_words)==36 and all(s.bit_count()==9 for s in triple_words)
    keep = keep_pairs + sum(not reject(s) for s in triple_words)*w(3,0)
    return {"density": (F(0), remainder, remainder),
            "cluster": (keep, keep+tail_pairs+remainder, keep+tail_pairs+remainder)}


def one(index):
    checks, pairs, triples = geometry()
    bern = index == -1
    if bern:
        p,h = F(1,10**6), F(1,10**4)
        norm = p.denominator
        laws = [({0:norm-p.numerator},{0:p.numerator}) for _ in checks]
        priors = detector_priors(checks,h,rate=p)
        base = load(FLAG/"pilot-results.json")["flag225_direct"]
        budget = F(base["aborted"][0])
        label = {"model":"Bernoulli","p":"1e-6","q":"1e-4"}
    else:
        row = load(COMP/"piloto_HGP_225_esquinas_20261124.json")["circuito"][index]
        orig = load(FLAG/"circuit-flag.json")["rows"][index]
        rc,rm = F(row["rc"]),F(row["rm"])
        cc = row["orden_ejecutado"] == "CC primero"
        assert (orig["rc"],orig["rm"],orig["order"]) == (row["rc"],row["rm"],"CC first" if cc else "VV first")
        schedule,_ = CIRC["schedule"](checks,cc_first=cc)
        laws,norm = CIRC["joint_fibres"](schedule,rc)
        h = CIRC["ancilla_rate"](rc,rm)
        priors = detector_priors(checks,h,schedule,rc)
        check_circuit_priors(priors,schedule,rc,rm)
        budget = F(orig["discards_per_million"][0])/10**6
        label = {"model":"circuit","order":orig["order"],"rc":row["rc"],"rm":row["rm"]}
    wires = projected_wires(checks,laws)
    assert all(sum(row.values()) == norm for row in wires)
    # Independent check of the projected graph against the archived joint fibre laws.
    for t in (*BITS, sum(BITS[::2]), sum(BITS[::3]), sum(BITS)):
        wire_char = prod(sum(n*(-1)**((s&t).bit_count()%2) for s,n in w.items())/norm for w in wires)
        target = wire_char*(1-2*float(h))**t.bit_count()
        actual = prod(1-2*v for s,v in priors.items() if (s&t).bit_count()%2)
        assert abs(target-actual)<5e-13
    for parity, t in ((1, 9), ((1<<N)-1,0),
                      (sum(1<<j for j in range(0,N,3)),sum(BITS[::5]))):
        joint = prod(sum(m*(-1)**(((t & (d^(checks[v] if e else 0))).bit_count()
                                    + e*((parity >> v)&1))%2)
                         for e,table in enumerate(pair) for d,m in table.items())/norm
                     for v,pair in enumerate(laws))
        actual = prod(1-2*v for s,v in priors.items() if (s & (t | (parity<<M))).bit_count()%2)
        assert abs(joint*(1-2*float(h))**t.bit_count()-actual)<5e-13
    score = scorer(priors)
    z,scale,odds = wire_odds(wires,norm)
    density = density_discards(h,z,scale,odds)
    cut,cluster = cluster_cut(score,h,z,scale,odds,budget,not bern and rm==F(1,1000))
    def rejects(word):
        return ceil(GRID*score(word))>cut
    if index in (-1,6):
        cluster = refine_discards(rejects,h,z,scale,odds,cluster,bern)
    decisions = pair_decisions(pairs,checks,rejects)
    risks = (bernoulli_risks(p,h,checks,triples,decisions,rejects) if bern else
             circuit_risks(laws,norm,h,checks,pairs,triples,decisions,rejects,row))
    output = {**label, "observed_graph_columns":len(priors), "cluster_cutoff_millionths":cut,
              "cluster_rejected_21_and_22":list(map(sum,zip(*decisions.values())))}
    for name, (lo,hi) in (("density", density),("cluster",cluster)):
        assert hi <= budget
        lower,det,stat = risks[name]
        assert 0 <= lower <= stat <= det
        output[name] = {"discard_per_million":endpoints(10**6*lo,10**6*hi),
                        "risk_deterministic":endpoints(lower/(1-lo),det/(1-hi)),
                        "risk_stratified":endpoints(lower/(1-lo),stat/(1-hi))}
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--point", type=int, choices=range(-1,8))
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    indices = list(range(-1,8)) if args.point is None else [args.point]
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        result = []
        for row in pool.map(one,indices):
            result.append(row)
            print(row["model"], row.get("order",""), row.get("rm",""), row["cluster_cutoff_millionths"], flush=True)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
