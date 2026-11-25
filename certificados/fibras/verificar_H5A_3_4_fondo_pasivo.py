"""Certificado uniforme HGP (3,4) con el fondo pasivo, Python 3.10+."""

from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
from runpy import run_path


def step_costs(profiles, far):
    # B_n(j) es la carga máxima por un check en su primer apagado.
    ranks = [sum(max(Q(far*(v-h)+1,h) for v,h,_,_ in profiles if n*h >= j*v)
                 for j in range(1,n+1)) for n in (3,4)]

    def pair(n, m):
        @lru_cache(None)
        def credit(r, s):
            if r+s == 0:
                return 0
            possible = []
            for v,h,row,col in profiles:
                if n*h < r*v or m*h < s*v:
                    continue
                minimum = (v+1)//2 if r == n or s == m else 1
                for i in range(min(r, row if n == 3 else col, h)+1):
                    for j in range(min(s, row if m == 3 else col, h-i)+1):
                        if i+j >= minimum:
                            possible.append(far*(2*h-v)-1+credit(r-i,s-j))
            return min(possible)

        costs = [[], []]
        for overlap in range(1+(min(n,m) if n == m else 1)):
            for t in range(1+min(n-overlap,n//2)):
                for u in range(1+min(m-overlap,m//2)):
                    r, s = n-overlap-t, m-overlap-u
                    balance = r+s-t-u
                    costs[bool(overlap)].append(far*balance-credit(r,s))
        return tuple(map(max,costs))

    return ranks, [pair(3,3), pair(3,4), pair(4,4)]


def certificate(exp_upper, endpoint_counts, coverage, ranks, costs):
    p, q, alpha, eta, pin = map(Q, ("2.5e-9",".001",".001",".0001",".00117"))
    def radius(r):
        return r*exp_upper(Q(".015"))/(1-2*r)
    data = p*exp_upper(Q(".015")+alpha)/(1-2*p)
    measure = radius(q)
    # Ancla visitada, etiqueta activa, precios y ramas del árbol de un núcleo.
    U,s,C1,C2,d1,d2,u1,u2,un,uf,vn,vf,cn,cf = map(Q, (
        "0.00000873550119132 0.00113134456171 4.40375523367e-9 5.14968894171e-9 "
        "0.00000260340392659 0.0000354377272028 1.33681829991e-8 4.01733813290e-8 "
        "0.00000112668679409 0.00000104496346482 0.00608124643753 0.0000986786301680 "
        "5.31583055616e-8 4.93025102055e-8").split())
    assert U < eta and U < pin < 1
    local_mass = 89*pin+469*U
    error = local_mass*exp_upper(2*local_mass)
    assert error < 1
    passive_ratio = measure*(1+pin)*(1+error)/pin
    assert passive_ratio < Q(".989")

    root = Q(".9729")
    base = root**5
    near, far = base**3, base**33
    assert far == near**11 and (root**11)**15 == far and 36*far < 1
    j,f1,f2,k1,k2,m1,m2,mg = endpoint_counts(root**11)
    assert (j,f1,f2,k1,k2,m1,m2,mg) == (152,2170,2574,2893,1931,352,421,968)
    assert coverage == (414,536) and ranks == [14,Q(59,3)]
    assert costs == [(14,24),(14,24),(25,35)]
    B = 1297*U
    extra = un+uf
    R1,R2 = 640*d1+441*d2,784*d1+513*d2
    x1,x2 = C1/near**14, C2/base**59
    V1 = (1+s)**3*exp_upper(1344*extra+30*vn+near*(f1*vf+m1*s))
    V2 = (1+s)**4*exp_upper(1596*extra+36*vn+near*(f2*vf+m2*s))

    def pair_price(i, ordinary, common):
        k, l = costs[i]
        return ordinary/near**k + common/near**l

    pair1 = V1*C1*C1*pair_price(0,640,9) + V2*C1*C2*pair_price(1,441,9)
    pair2 = V1*C1*C2*pair_price(1,784,16) + V2*C2*C2*pair_price(2,513,8)
    exp1, exp2 = exp_upper(R1), exp_upper(R2)
    child1, child2 = V1*x1*(exp1-1), V2*x2*(exp2-1)
    large1 = exp1-1-R1+640*child1+441*child2
    large2 = exp2-1-R2+784*child1+513*child2
    L = j*vf+far*mg*s
    G = (1+s)**12*exp_upper(1792*u1+1197*u2+3648*extra+L)
    Hn, Hf = 40*u1+27*u2, j*(cn+cf)+near*(k1*u1+k2*u2)
    ratios = (
        max(u1+3*extra,u2+4*extra)/U,
        measure*exp_upper(90*pin)/s,
        data*exp_upper(eta+Q(17,3)*B+pin*coverage[0])/C1,
        data*exp_upper(eta+Q(103,15)*B+pin*coverage[1])/C2,
        V1*x1*exp1/d1, V2*x2*exp2/d2,
        V1*(C1/near+pair1+x1*large1)/u1, V2*(C2/near+pair2+x2*large2)/u2,
        near*Hn*G/un, far*Hf*G/uf,
        near*(G-1-L-12*s)/vn, far*(G-1-L-12*s)/vf,
        near*Hn*(G-1)/cn, far*Hf*(G-1)/cf,
    )
    assert min(ratios) > 0 and max(ratios) < Q(".9991")
    row = max(640*u1+441*u2, 784*u1+513*u2)
    assert row < Q(".00003110") < alpha
    intermediate = Q("2e-9")
    scale = intermediate/p
    assert scale == Q(4,5) and scale**3*row < Q(".00001593")
    # Techo del mismo recuento al precio de 11 por generador lejano,
    # para todo precio de base. e^(1+536q) se acota por debajo.
    v = 1+536*q
    lower = 1+v+v*v/2+v**3/6+v**4/24+v**5/120+v**6/720
    assert (Q("7e-7")*513*lower)**33 * 36**59 > 1
    print(f"p=2.5e-9, q=.001. Respuesta pasiva {float(passive_ratio):.10f}. "
          f"Máximo de los 14 cocientes de núcleos {float(max(ratios)):.10f}.")
    print(f"Fila completa {float(row):.12g} < .00003110, presupuesto .001. "
          f"En p=2e-9, mayorante de tres o más {float(scale**3*row):.12g} < .00001593.")


def check_shared_background():
    """Contraste de la biyección en el producto de la incidencia afín de orden 3."""
    from itertools import combinations

    def bits(mask):
        while mask:
            bit = mask & -mask
            yield bit.bit_length()-1
            mask ^= bit

    def subsets(items):
        return [sum(1 << i for i in subset) for k in range(len(items)+1)
                for subset in combinations(items,k)]

    lines = [tuple(3*x+(m*x+b)%3 for x in range(3)) for m in range(3) for b in range(3)]
    lines += [tuple(3*x+y for y in range(3)) for x in range(3)]
    incident = [[i for i,line in enumerate(lines) if c in line] for c in range(9)]
    psi = [sum(1 << (9*i+c) for c in lines[j]) for i in range(12) for j in range(12)]
    psi += [sum(1 << (9*i+c) for i in incident[r]) for r in range(9) for c in range(9)]
    assert len(psi) == 225 and set(map(int.bit_count,psi)) == {3,4}
    rectangles, candidates = [], {}
    for r in range(9):
        for j in range(12):
            qubits = [12*i+j for i in incident[r]] + [144+9*r+c for c in lines[j]]
            rectangles.append(sum(1 << (9*i+c) for i in incident[r] for c in lines[j]))
            for f in subsets(qubits)[1:]:
                phi = 0
                for v in bits(f):
                    phi ^= psi[v]
                if phi and 2*phi.bit_count() >= 3*f.bit_count():
                    candidates[f] = phi
    # Una prioridad física fija. Entre huellas idénticas sólo gana la primera.
    best_rep = {}
    for f in sorted(candidates,key=lambda f: (-f.bit_count(),f)):
        best_rep.setdefault(candidates[f], f)
    actions = [(phi, phi.bit_count(), f) for phi,f in best_rep.items()]
    readers = [sum(1 << k for k,r in enumerate(rectangles) if r & p) for p in psi]

    @lru_cache(None)
    def run(s):
        best, lit, area = -1, 1, 2
        for j,(phi,n,_) in enumerate(actions):
            h = (s & phi).bit_count()
            if h*area > lit*n:
                best,lit,area = j,h,n
        if best < 0:
            return ()
        assert 100*(2*lit-area) >= area
        return (best,)+run(s ^ actions[best][0])

    def connected(w):
        sites = list(bits(w))
        reached = readers[sites[0]]
        while True:
            old = reached
            for v in sites:
                if reached & readers[v]:
                    reached |= readers[v]
            if reached == old:
                return all(reached & readers[v] for v in sites), reached

    # Tres primeras coordenadas paralelas; hay medidas que leen varios grupos.
    inputs = subsets([0, 3, 12, 24])
    measurements = subsets([0,3,9,18,27,39,51,57])
    stable = {d for d in measurements if not run(d)}
    direct, cores = {}, []
    for e in inputs:
        sigma = original = 0
        for x in bits(e):
            sigma ^= psi[x]
            original |= psi[x]
        for d in measurements:
            trace = run(sigma ^ d)
            residual, active, w = sigma,original,e
            for j in trace:
                phi,_,f = actions[j]
                residual ^= phi
                active |= phi
                w |= f
            if residual:
                continue
            direct[e,d] = trace
            if e and not (d & ~active):
                one, ports = connected(w)
                if one:
                    allowed = {p for p in stable if not (p & active) and run(sigma ^ d ^ p) == trace}
                    assert 0 in allowed
                    assert all((p ^ (1 << i)) in allowed for p in allowed for i in bits(p))
                    halo = 0
                    for rect in rectangles:
                        if rect & active:
                            halo |= rect
                    projections = {p & halo for p in allowed}
                    assert all((p in allowed) == ((p & halo) in projections)
                               for p in stable if not (p & active))
                    cores.append((e,d,w,ports,trace,allowed,halo))
    reconstructed, counts = set(), [0,0]

    def compose(start, chosen, e, d, ports, allowed):
        for p in allowed:
            key = (e, d|p)
            assert key in direct and key not in reconstructed
            reconstructed.add(key)
            actual = direct[key]
            assert len(actual) == sum(len(c[4]) for c in chosen)
            for c in chosen:
                assert tuple(j for j in actual if actions[j][2] & c[2]) == c[4]
            if len(chosen) >= 3:
                counts[0] += 1
                if any(sum(bool(c[6] & (1 << j)) for c in chosen) >= 2 for j in bits(p)):
                    counts[1] += 1
        for k in range(start,len(cores)):
            c = cores[k]
            if not (ports & c[3]):
                compose(k+1, chosen+[c], e|c[0], d|c[1], ports|c[3], allowed & c[5])
    compose(0,[],0,0,0,stable)
    assert reconstructed == direct.keys() and counts[1] > 0
    print(f"Contraste HGP de 225 qubits: {len(inputs)*len(measurements)} entradas, "
          f"{len(direct)} filtradas, {counts[0]} con tres o más núcleos, "
          f"{counts[1]} de ellas con una medida pasiva en varias vecindades de lectura.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para comprobar las desigualdades.")
    folder = Path(__file__).resolve().parent
    first = run_path(str(folder/"verificar_H5A_3_4_medidas_y_un_dato.py"))
    adjacent = run_path(str(folder/"verificar_H5A_3_4_pares_adyacentes.py"))
    corners = run_path(str(folder/"verificar_H5A_3_4_esquinas.py"))
    two = run_path(str(folder/"verificar_H5A_3_4_dos_componentes.py"))
    profiles = adjacent["stripe_profiles"]()
    assert step_costs(profiles,5)[1] == adjacent["component_budgets"](profiles)
    coverage = two["local_checks"](first["budgets"])
    ranks, costs = step_costs(profiles,11)
    certificate(first["exp_upper"], corners["endpoint_counts"], coverage, ranks, costs)
    check_shared_background()
