"""Verify the twelve inequalities printed in Appendix B.3 from its finite count inputs."""
from fractions import Fraction as Q
from math import ceil


def round_up(value, digits):
    grid = 10**digits
    return Q(ceil(value * grid), grid)


def exp_up(x):
    assert 0 <= x < 26
    term = total = Q(1)
    for k in range(1, 25):
        term *= x / k
        total += term
    return round_up(total + term * x / (25 * (1 - x/26)), 26)


def dot(row, vector):
    return sum(a*b for a, b in zip(row, vector, strict=True))


W = tuple(tuple(map(Q, row)) for row in (
    ("479.413463293782", "319.712053370976"), ("568.376983770624", "379.040194911360")))
S = tuple(tuple(map(Q, row)) for row in (
    ("146.995636562163", "112.725959011272"), ("200.401704908928", "123.178993446112")))
a = Q(".568")
separated = 1 + round_up(a**-4 - 1, 13)
strong = [round_up(a**-(7+3*i), 13) for i in range(2)]
surcharge = [round_up((a**-(7+i) - 1)/2, 13) for i in range(2)]


def free_supersolution(weights):
    state = [[Q(0), Q(0)] for _ in range(5)]
    for _ in range(500):
        free, attached, core, satellite, premium = state
        new = [[] for _ in state]
        for i, x in enumerate(weights):
            g = dot(W[i], free)/2
            b = g + dot(W[i], premium)/2
            forest = exp_up(g+b)
            strong_root = x * strong[i] * forest * exp_up(dot(W[i], satellite))
            strong_child = exp_up(dot(S[i], attached))
            new_core = strong_root * (strong_child-1)
            new_satellite = x * surcharge[i] * forest
            new_free = x*exp_up(g)*(1 + separated*(exp_up(b)-1)) + new_core + new_satellite*dot(W[i], core)
            updates = (new_free, strong_root*strong_child, new_core, new_satellite, x*(separated-1)*exp_up(g))
            for row, value in zip(new, updates, strict=True):
                row.append(round_up(value, 24))
        assert all(lo <= hi for old, nxt in zip(state, new) for lo, hi in zip(old, nxt))
        if new == state:
            return free
        state = new
    raise ArithmeticError("Positive iteration did not produce a supersolution.")


def free_coefficients(weights, degree):
    """Compute the coefficients of U(t weights) by total degree."""
    def series():
        return [[Q(0), Q(0)] for _ in range(degree+1)]

    free, attached, core, satellite, premium = (series() for _ in range(5))
    logs = [series() for _ in range(4)]
    exps = [series() for _ in logs]
    for s in exps:
        s[0] = [Q(1), Q(1)]
    for n in range(1, degree+1):
        for i, x in enumerate(weights):
            unit, forest, with_satellites, with_strong = (s[n-1][i] for s in exps)
            core[n][i] = x*strong[i]*(with_strong-with_satellites)
            attached[n][i] = x*strong[i]*with_strong
            satellite[n][i] = x*surcharge[i]*forest
            premium[n][i] = x*(separated-1)*unit
            chosen = sum(satellite[k][i]*dot(W[i], core[n-k]) for k in range(1, n))
            free[n][i] = x*(unit + separated*(forest-unit)) + core[n][i] + chosen
        for i in range(2):
            weak = dot(W[i], free[n])
            logs[0][n][i] = weak/2
            logs[1][n][i] = weak + dot(W[i], premium[n])/2
            logs[2][n][i] = logs[1][n][i] + dot(W[i], satellite[n])
            logs[3][n][i] = logs[2][n][i] + dot(S[i], attached[n])
            for log, exponential in zip(logs, exps):
                exponential[n][i] = sum(k*log[k][i]*exponential[n-k][i] for k in range(1, n+1))/n
    return free


def check():
    m, u, h, H = ([Q(v)/10**6 for v in row] for row in (
        ("1.103585668", "1.125499296"), ("1.064424176", "1.087972930"),
        (".04753509019", ".1387925312"), (".006387261997", ".02275491635")))
    visible, root = Q(".6862867109e-6"), Q("2.337412703e-6")
    b, theta, eta, alpha = map(Q, (".0001085", ".000108", ".0001", ".004"))
    q, p = Q("1e-4"), Q("1e-6")
    y = q*exp_up(Q(".015"))/(1-2*q)
    z = p*exp_up(Q(".015")+alpha)/(1-2*p)
    def log_bound(x):
        return x + x*x/(2*(1-x))
    A = log_bound(theta) - Q(89,2)*log_bound(theta**2)
    B = b + 89*log_bound(theta**2) - log_bound(theta)
    g = 40*H[0]+27*H[1]+root
    I = 280*u[0]+189*u[1]+588*g
    N = ((640,441), (784,513))
    R = max(dot(row, u) for row in N)
    phi = exp_up(89*b+I)
    epsilon = phi*(phi-1)
    assert A > 0 and B > 0 and y/(1+y) <= theta and R < alpha and max(u) < eta and epsilon < 1
    K = ((1028316,685678,21163), (1218982,812566,25086), (28217,18815,547))
    pointed = ((174179,116136,3657), (206463,137631,4335), (4876,3252,90))
    children = (*h, visible)
    w = [round_up(m[i]*exp_up(dot(K[i],children)), 26) for i in range(2)]
    gross = free_supersolution(w)
    head = free_coefficients(w, 10)
    naked = free_coefficients(m, 13)
    tail = [gross[i] - sum(head[k][i] for k in range(1,11)) for i in range(2)]
    naked_tail = [sum(naked[k][i] for k in range(11,14)) for i in range(2)]
    def counts(s):
        return [Q(v)/10**8 for v in s.split()]
    error = Q("1e-17")  # One last printed unit in the finite count inputs.
    for actual, printed in ((tail, counts(".016370937 .094225744")),
                            (naked_tail, counts(".002031421 .011364058"))):
        assert all(abs(v-w) <= error for v, w in zip(actual, printed))
    group = counts("4.732393468 13.765875275")
    bare = counts("2.726968511 6.229548460")
    distinguished = counts(".489174459 2.023387217")
    halo = ((414,1264,Q("859.5"),Q("2623.5")), (536,1652,1095,3374))
    simple = ((210,640,441,1344), (252,784,513,1596))
    rows = {"readouts": [y*(1+b)*(1+epsilon)/b, y*(1-theta)*exp_up(99*A+90*B+2*I)/theta],
            "source": [], "ordinary": [], "fixed": [], "site": []}
    for i in range(2):
        exposure = dot(halo[i], (B,2*u[0],2*u[1],2*g))
        M = z*exp_up(eta+2*R+(504,628)[i]*A+exposure)
        simple_mass = z*(1+y)**(3+i)*exp_up(eta+2*R+dot(simple[i],(b,2*u[0],2*u[1],2*g)))
        group_at_w, group_at_m = group[i]+error, bare[i]-error
        connected = exp_up(dot(N[i], [Q("1.01")*v for v in m]))
        assert connected <= Q("1.01")
        fixed = min(group_at_w-group_at_m, distinguished[i]+error) + tail[i]-naked_tail[i] + m[i]*(connected-1)
        rows["source"].append(M/m[i])
        rows["ordinary"].append((group_at_w+tail[i])/h[i])
        rows["fixed"].append(fixed/H[i])
        rows["site"].append((simple_mass+H[i])/u[i])
    J, d = dot(K[2], children), dot(pointed[2], children)
    omega = a**6/Q(".678")**2
    remaining = exp_up(J)*d - sum(pointed[2][i]*(bare[i]-error+naked_tail[i]) for i in range(2))
    assert remaining > 0
    generator = omega*(omega*remaining + (exp_up(J)-1)*d/2) + Q("1.6902840e-8")
    rows["generator"] = [omega*(exp_up(J)-1-J)/visible, generator/root]
    caps = ((".947357",".954746"), (".99952024",".99952592"), (".999002",".998621"),
            (".998307",".997744"), (".999509",".999491"), (".998418",".998614"))
    for (label, pair), bounds in zip(rows.items(), caps, strict=True):
        assert all(0 < x <= Q(bnd) < Q(".999526") for x,bnd in zip(pair,bounds))
        print(label, *(f"{float(x):.12f}" for x in pair))


if __name__ == "__main__":
    check()
