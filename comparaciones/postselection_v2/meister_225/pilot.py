"""Run the one-point cluster-geometry pilot at the native Lee discard rate."""

from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
import argparse
import json
import platform
import random

PROCESS_START = perf_counter()

from inputs import (LEE_CUT, P, Q, collect, input_catalogues, load_comparison,
                    read_features, sector_weight, words_for)
from extra_silent import two_data_three_readouts
from lee_certificate import lee_silent_bounds
from logical_gap import LogicalGap, set_bits


def calibrated_boundary(gap_mass, lee_mass, missing):
    # The true Lee rate, mass below the boundary, and boundary mass all
    # include their as yet unclassified shares of "missing".
    below = F(0)
    for boundary, equal in enumerate(gap_mass):
        if equal and below + missing <= lee_mass and lee_mass + missing <= below + equal:
            alpha_lo = (lee_mass - below - missing) / (equal + missing)
            alpha_hi = (lee_mass + missing - below) / equal
            assert 0 <= alpha_lo <= alpha_hi <= 1
            return boundary, alpha_lo, alpha_hi
        below += equal
    raise ValueError("The marginal tail does not isolate one cutoff; refine the input catalogue.")


def check_metric(metric, erasures):
    # Independent global distance bound on subsets of a known minimum logical.
    w = list(set_bits(metric.witness))
    for size in range(7):
        assert metric(sum(1 << v for v in w[:size])) == 6-size
    assert metric(1 << 144) == 6
    # For erasures containing some VV site the verified row witness has
    # cost <=5. Compare geometric shortcuts with an unrestricted quotient
    # search on a repeatable sample from all observed erasure sizes.
    rng = random.Random(22507433)
    checked = 0
    for size in sorted({e.bit_count() for e in erasures}):
        usable = sorted(e for e in erasures if e.bit_count() == size and e % (1 << 144))
        for erasure in rng.sample(usable, min(12, len(usable))):
            assert metric(erasure) == metric.quotient_gap(erasure, 5)
            checked += 1
    return checked


def analyse(repo, cache):
    started = perf_counter()
    compare = load_comparison(repo)
    h, checks, tags, marginal, silent, triples, ssf, patterns = input_catalogues(compare)
    words = words_for(marginal, silent)
    features, meta = read_features(cache, words)
    began_gap = perf_counter()
    metric = LogicalGap(h, checks, tags)
    erasures = {entry[0] for entry in features.values()}
    gaps = {e:metric(e) for e in erasures}
    metric_stats = dict(metric.counts)
    gap_seconds = perf_counter()-began_gap
    checked = check_metric(metric, erasures)
    projected = {word:(gaps[e], lee) for word,(e,lee) in features.items()}

    gaps_mass = [F(0) for _ in range(7)]
    lee_mass = F(0)
    for sector, catalogue in marginal.items():
        counts = Counter()
        lee_count = 0
        for word, n in catalogue.items():
            gap, lee = projected[word]
            counts[gap] += n
            lee_count += n * lee
        w = sector_weight(*sector)
        for gap, count in counts.items():
            gaps_mass[gap] += count*w
        lee_mass += w*lee_count
    missing = 1-sum(gaps_mass)
    assert 0 < missing < F(2, 10**8)
    cap_lo, cap_hi = lee_mass, lee_mass+missing
    boundary, alpha_lo, alpha_hi = calibrated_boundary(gaps_mass, lee_mass, missing)
    ppm = 10**6
    endpoints = compare.endpoints

    # Use precisely the archived deterministic upper bound on all remaining
    # silent sectors, by telling the existing routine that every enumerated
    # silent history was discarded. Its contribution on the known part is zero.
    rejected_all = {j:(1,99) for j in range(324)}
    original = compare.bernoulli_risks(P,Q,checks,triples,rejected_all,lambda _:True)
    zero, original_remainder, _ = original["cluster"]
    assert zero == 0
    cap31 = compare.load(compare.FLAG/"pilot-results.json")["orbit_bounds"]["225"]["3,1"]["count_upper"]
    assert cap31 == 101088 and silent[3,1].total() <= cap31
    began_census = perf_counter()
    candidates23, silent23 = two_data_three_readouts(patterns, checks, tags, ssf)
    census_seconds = perf_counter()-began_census
    # Pay every silent (2,3) history conservatively in the remaining tail.
    remainder = (original_remainder - cap31*sector_weight(3,1)
                 - (candidates23-silent23)*sector_weight(2,3))
    assert remainder > 0
    known_lee = known_keep = known_ties = F(0)
    old_known_lee = F(0)
    reports = {}
    for sector,catalogue in silent.items():
        hist, joint = Counter(), Counter()
        kept_lee = 0
        for word,n in catalogue.items():
            gap,lee = projected[word]
            hist[gap] += n
            joint[gap,lee] += n
            kept_lee += n * (1-lee)
        kept = sum(n for gap,n in hist.items() if gap > boundary)
        ties = hist[boundary]
        w = sector_weight(*sector)
        known_lee += w*kept_lee
        if sector != (3,1):
            old_known_lee += w*kept_lee
        known_keep += w*kept
        known_ties += w*ties
        reports[",".join(map(str,sector))] = {
            "total":catalogue.total(), "lee_retains":kept_lee, "gap_histogram":dict(sorted(hist.items())),
            "always_retained":kept, "on_randomized_boundary":ties,
            "counts_by_gap_and_lee_discards":{f"{g},{lee}":m for (g,lee),m in sorted(joint.items())}}
    known_m_lo = known_keep + (1-alpha_hi)*known_ties
    known_m_hi = known_keep + (1-alpha_lo)*known_ties
    simple, logicals, five_candidates, witnessed, lee_joint = lee_silent_bounds(
        compare,h,checks,tags,ssf,features)
    lee_lo = max(known_lee,lee_joint[0])/(1-cap_lo)
    lee_hi = min(known_lee+remainder,lee_joint[1])/(1-cap_hi)
    geom_lo, geom_hi = known_m_lo/(1-cap_lo), (known_m_hi+remainder)/(1-cap_hi)
    reference = compare.load(compare.HERE/"comparison.json")[0]
    assert reference["cluster_cutoff_millionths"] == LEE_CUT
    assert endpoints(ppm*cap_lo, ppm*cap_hi) == reference["cluster"]["discard_per_million"]
    assert endpoints(old_known_lee/(1-cap_lo),(old_known_lee+original_remainder)/(1-cap_hi)) == reference["cluster"]["risk_deterministic"]
    # The original density remainder charged every tie-relaxed (3,1)
    # candidate. Replace just that term with the exact fixed-order count
    # having fewer than seven lit checks.
    density_candidates = Counter()
    for word, count in silent[3,1].items():
        density_candidates[word.bit_count()] += count
    density_count = sum(n for size,n in density_candidates.items() if size < 7)
    density_known = density_count*sector_weight(3,1)
    density_upper = original["density"][1]-(cap31-density_count)*sector_weight(3,1)
    density_discard = [F(x)/ppm for x in reference["density"]["discard_per_million"]]
    density_risk = endpoints(density_known/(1-density_discard[0]),
                             density_upper/(1-density_discard[1]))
    # With two clustered data bits, the Lee norm is already above cutoff,
    # regardless of cluster partition. A geometry gap below 5 requires at
    # least two such bits since the code distance is six. Bound the weight
    # ratio rigorously using 1 < log(999999)/log(9999) < 2:
    # sqrt(2)*w_p / (225*w_p+108*w_q) > sqrt(2)/333 > 0.003368.
    assert boundary == 5 and 2*10**12 > (333*LEE_CUT)**2
    # If Lee accepts an unclassified silent history, geometry can discard
    # it with probability at most alpha_hi (only the gap-five boundary).
    delta_lo = known_m_lo-known_lee-alpha_hi*remainder
    delta_hi = known_m_hi-known_lee+remainder
    paired = [n/(1-a) for n in (delta_lo,delta_hi) for a in (cap_lo,cap_hi)]
    diff_lo = max(geom_lo-lee_hi,min(paired))
    diff_hi = min(geom_hi-lee_lo,max(paired))
    assert diff_lo <= diff_hi
    ratio_lower = geom_lo/lee_hi
    if diff_lo > 0:
        verdict = "Lee has strictly lower retained risk at the same discard probability."
        ratio_lower = max(ratio_lower,1+diff_lo/lee_hi)
    elif diff_hi < 0:
        verdict = "The cluster-geometry rule has strictly lower retained risk at equal discard."
    else:
        verdict = "The certified comparison does not order the retained risks."
    example = checks[0]^checks[3]^checks[4]^(1 << 3)
    assert projected[example] == (6,1) and features[example][0] == 0
    assert ssf(example)[0] == 1 << 3 and ssf(example)[1] != tags[0]^tags[3]^tags[4]
    wall = perf_counter()-started
    report = {
        "code":"affine [[225,9,6]]", "noise":"one-round independent Bernoulli, X sector",
        "p":str(P),"q":str(Q),"lee_cutoff_millionths":LEE_CUT,
        "discard_per_million_both":endpoints(ppm*cap_lo,ppm*cap_hi),
        "discard_calibration":{
            "unclassified_input_mass":endpoints(missing,missing),
            "known_gap_mass_per_million":{
                str(i):endpoints(ppm*v,ppm*v) for i,v in enumerate(gaps_mass)},
            "randomized_gap":boundary,
            "boundary_discard_probability":endpoints(alpha_lo,alpha_hi)},
        "enumerated_silent_histories":reports,
        "two_data_three_readouts":{"candidates":candidates23,"silent":silent23},
        "lee_certificate":{"admissible_input_superset":simple,"weight_six_logicals":logicals,
            "possible_weight_five_data_supports":five_candidates,"confirmed_at_zero_readout":witnessed,
            "upper_using_only_distance":endpoints(0,lee_joint[2]/(1-cap_hi))[1],
            "risk_factor_using_only_distance_lower":endpoints(geom_lo/(lee_joint[2]/(1-cap_hi)),
                geom_lo/(lee_joint[2]/(1-cap_hi)))[0]},
        "disagreeing_example":{"data_sites":[0,3,4],"readout_check":3,
            "lit_input_checks":list(set_bits(example)),"clustered_data_sites":[],
            "gap":6,"lee_discards":True},
        "unclassified_silent_mass_upper":endpoints(remainder,remainder)[1],
        "seven_lit_checks":{"silent_31_by_input_weight":dict(sorted(density_candidates.items())),
                            "retained_31":density_count,
                            "discard_per_million":reference["density"]["discard_per_million"],
                            "risk_per_retained_round":density_risk},
        "risk_per_retained_round":{"lee":endpoints(lee_lo,lee_hi),
            "cluster_geometry":endpoints(geom_lo,geom_hi)},
        "geometry_minus_lee_risk":endpoints_signed(diff_lo,diff_hi,endpoints),
        "paired_difference_without_lee_support_certificate":endpoints_signed(min(paired),max(paired),endpoints),
        "geometry_to_lee_risk_lower":endpoints(ratio_lower,ratio_lower)[0],
        "verdict":verdict,
        "work":{"unique_syndromes":len(words),"unique_data_erasure_sets":len(erasures),
            "metric_solver_calls":metric_stats,"crosschecked_erasure_patterns":checked},
        "runtime":{"gap_calculation_seconds":gap_seconds,
            "two_data_three_readouts_seconds":census_seconds,
            "input_geometry_seconds":meta["geometry_seconds"],
            "common_features_seconds":meta["feature_seconds"],
            "feature_stage_seconds":meta["elapsed_seconds"],"analysis_stage_seconds":wall,
            "total_measured_seconds":meta["elapsed_seconds"]+wall},
    }
    return report


def endpoints_signed(lo, hi, endpoints):
    def one(x, lower):
        ends = endpoints(abs(x), abs(x))
        if x < 0:
            return "-" + ends[1 if lower else 0]
        return ends[0 if lower else 1]
    return [one(lo,True),one(hi,False)]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reuse-cache", action="store_true")
    args = parser.parse_args()
    if not args.reuse_cache:
        collect(args.repo, args.cache)
    result = analyse(args.repo, args.cache)
    result["runtime"].update(command_seconds=perf_counter()-PROCESS_START,
                          reused_features=args.reuse_cache, python=platform.python_version(),
                          processor=platform.processor())
    args.output.write_text(json.dumps(result, indent=2)+'\n',encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2))
