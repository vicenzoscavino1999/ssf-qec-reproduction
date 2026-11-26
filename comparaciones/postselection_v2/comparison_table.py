"""Typeset the numerical rows of the postselection table with outward rounding."""
import json
from decimal import Decimal as D, ROUND_CEILING, ROUND_FLOOR
from pathlib import Path

HERE = Path(__file__).resolve().parent
FLAGS = HERE.parents[1]/"sectores/flag_rates"


def rounded(x, upper=False):
    if not x:
        return "0"
    unit = D(10)**(x.adjusted()-2)
    return format(x.quantize(unit, rounding=ROUND_CEILING if upper else ROUND_FLOOR), "f")


def interval(values, scale=1):
    lo,hi = (D(x)*scale for x in values)
    return f"[{rounded(lo)},{rounded(hi,True)}]"


def rows():
    inputs = json.loads((HERE/"comparison.json").read_bytes())
    bernoulli = json.loads((HERE/"meister_225/results.json").read_bytes())
    base = json.loads((FLAGS/"circuit-flag.json").read_bytes())["rows"]
    bflag = json.loads((FLAGS/"pilot-results.json").read_bytes())["flag225_direct"]
    result = []
    for i, point in enumerate(inputs):
        if i:
            ssf = base[i-1]
            assert all(ssf[k]==point[k] for k in ("order","rc","rm"))
            data = format(D(point["rc"])*D(10)**7, ".3f")
            readout = D(point["rm"])
            label = "VV" if point["order"]=="VV first" else "CC"
            discard, risk = ssf["discards_per_million"], ssf["risk_per_retained_round"]
        else:
            label, data, readout = "Bern.", "10", D(point["q"])
            discard = [D(x)*10**6 for x in bflag["aborted"]]
            risk = [bflag["silent_unflagged"][0], bflag["silent_per_accepted_upper"]]
        density,cluster = point["density"],point["cluster"]
        assert D(cluster["discard_per_million"][1])<=D(discard[0])
        assert D(density["discard_per_million"][1])<=D(discard[0])
        assert point["cluster_rejected_21_and_22"]==[324,32076]
        density_risk = density["risk_stratified"]
        cluster_upper = rounded(D(cluster["risk_stratified"][1])*10**18, True)
        if not i:
            assert (bernoulli["discard_per_million_both"] == cluster["discard_per_million"]
                    and bernoulli["seven_lit_checks"]["discard_per_million"] == density["discard_per_million"])
            density_risk = bernoulli["seven_lit_checks"]["risk_per_retained_round"]
            bound = D(bernoulli["risk_per_retained_round"]["lee"][1])*10**18
            exponent = bound.adjusted()
            cluster_upper = "$"+rounded(bound.scaleb(-exponent), True)+r"\cdot10^{"+str(exponent)+"}$"
        digits = (label, "$10^{"+str(readout.adjusted())+"}$", data,
                  interval(discard), interval(risk,10**18),
                  interval(density["discard_per_million"]), interval(density_risk,10**18),
                  str(point["cluster_cutoff_millionths"]), interval(cluster["discard_per_million"]),
                  cluster_upper)
        result.append(" & ".join(digits)+r" \\")
    return "\n".join(result)


if __name__=="__main__":
    print(rows())
