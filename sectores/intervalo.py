"""Exact rational bounds on the total Bernoulli probability, including every source weight."""
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import json

if __package__:
    from .mayorante import false_perfect, no_perfect
else:
    from mayorante import false_perfect, no_perfect


def render(x,direction=ROUND_CEILING):
    with localcontext() as ctx:
        ctx.prec=15
        ctx.rounding=direction
        return str(Decimal(x.numerator)/Decimal(x.denominator))


def regional_interval(qubits, p, q):
    """Corollary 6 for the paper's two Bernoulli codes, using rational rates."""
    remainder = {225: Q(".0180"), 900: Q(".1728")}
    if qubits not in remainder:
        raise ValueError("Choose the 225- or 900-qubit code.")
    if not 0 <= q <= Q(1, 10_000) or not 0 <= p < Q(1, 2):
        raise ValueError("The query requires 0 <= q <= 1e-4 and 0 <= p < 1/2.")
    if p == 0:
        return {"lower": "0", "upper": "0", "estimate": "0", "max_relative_percent": None}
    # Cross-multiplication also excludes q=0 when p>0.
    if p * (1 - q) > Q(9999, 999999) * q * (1 - p):
        raise ValueError("Outside Corollary 6. For p > 0, require q > 0 and "
                         "(p/(1-p))/(q/(1-q)) <= 9999/999999.")

    reference = Path(__file__).resolve().parent / "referencia" / f"resultado_{qubits}.json"
    data = json.loads(reference.read_text(encoding="utf-8-sig"))
    checks = data["geometry"]["checks"]
    count = data["interval"]["positive_census"]["2,1"]
    rho = remainder[qubits]
    t21 = count * p**2 * q * (1 - p)**(qubits - 2) * (1 - q)**(checks - 1)
    lower = render(t21, ROUND_FLOOR)
    upper = render(t21 * (1 + rho), ROUND_CEILING)
    estimate = render(t21 * 2 * (1 + rho) / (2 + rho), ROUND_HALF_EVEN)
    # Certify the printed estimate against both outward-rounded endpoints.
    point = Q(estimate)
    error = max(point / Q(lower) - 1, 1 - point / Q(upper))
    return {"lower": lower, "upper": upper, "estimate": estimate,
            "max_relative_percent": render(100 * error, ROUND_CEILING)}


def certificate(meta, counts, exceptional):
    n,m=meta['qubits'],meta['checks']
    A,B=12*meta['sheets'],9*meta['sheets']
    p,q=Q(1,10**6),Q(1,10**4)
    def w(a,b):
        return p**a*(1-p)**(n-a)*q**b*(1-q)**(m-b)
    size=lambda a,b:comb(n,a)*comb(m,b)
    # True initial removal to a failing smaller physical configuration.
    upper=dict(counts)
    upper[3,1]=min(size(3,1),exceptional+(n-2)*counts[2,1])
    for a,bb in ((3,(2,3)),(4,(0,1))):
        for b in bb:
            danger=false_perfect(n,m,a,b)+no_perfect(n,m,a,b)*size(a,b)
            pairs=upper.get((2,b),meta['marked_pairs']*comb(m,b))
            smaller=comb(n-2,a-2)*pairs
            if a==4:
                smaller+=(n-3)*upper[3,b]
            upper[a,b]=min(size(a,b),danger+smaller)
    assert counts[2,0]==0
    # For >=3 readout errors, at least one lies in the data row or two share
    # one of its 3A rectangles. The remaining readouts are unrestricted.
    pair_tail=meta['marked_pairs']*p**2*(1-p)**(n-2)*q**3*(
        B*comb(m-1,2)+3*A*comb(12,2)*(m-2))
    other_tails=comb(n,3)*p**3*comb(m,4)*q**4+comb(n,4)*p**4*comb(m,2)*q**2+comb(n,5)*p**5
    lo=sum(w(a,b)*v for (a,b),v in counts.items())
    hi=sum(w(a,b)*v for (a,b),v in upper.items())+pair_tail+other_tails
    mid=2*lo*hi/(lo+hi)
    err=(hi-lo)/(hi+lo)
    assert 0<lo<=mid<=hi<1 and err<Q(1,10)
    # Half the log-partition budget. The positive Taylor partial sum proves
    # exp(x) >= 2N(1+zeta)/zeta, using log(1+zeta) >= zeta/(1+zeta).
    zeta=Q(1,1280*(1+999999**6))
    source_order=6436 if n==225 else 6529
    x=Q(15,1000)*(source_order+1)
    term=total=Q(1)
    for j in range(1,257):
        term*=x/j
        total+=term
    assert total*zeta>=2*(n+m)*(1+zeta)
    return dict(lower=render(lo,ROUND_FLOOR),upper=render(hi),estimate=render(mid),
                max_relative_percent=render(100*err),pair_tail_upper=render(pair_tail),
                other_tails_upper=render(other_tails),first_step_exceptional_triples=exceptional,
                weighted_upper={f'{a},{b}':render(w(a,b)*v) for (a,b),v in sorted(upper.items())},
                uniform_sufficient_source_order=source_order,
                positive_census={f'{a},{b}':v for (a,b),v in counts.items()})

