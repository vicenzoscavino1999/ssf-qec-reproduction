"""H5-A (3,4). Pares adyacentes, rectángulo uniforme y filas de cumulantes.

Python 3.10 o posterior. Todos los certificados usan enteros o Fraction.
"""

from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
from runpy import run_path


def stripe_profiles():
    rows = [7 << (3*i) for i in range(4)]
    cols = [sum(1 << (3*i+j) for i in range(4)) for j in range(3)]
    candidates = []
    for selection in range(1, 128):
        phi = 0
        for i, stripe in enumerate(rows+cols):
            if selection & (1 << i):
                phi ^= stripe
        size = phi.bit_count()
        if size and 2*size >= 3*selection.bit_count():
            assert selection.bit_count() <= 5
            candidates.append((phi, size, max((phi&r).bit_count() for r in rows),
                                max((phi&c).bit_count() for c in cols)))
    assert len(candidates) == 101
    assert min(v for _, v, _, _ in candidates) == 3
    result = set()
    for state in range(4096):
        best_h, best_v = -1, 1
        tied = []
        for phi, v, row_cap, col_cap in candidates:
            h = (phi & state).bit_count()
            cmp = h*best_v-best_h*v
            if cmp > 0:
                best_h, best_v, tied = h, v, []
            if cmp >= 0:
                tied.append((v, h, row_cap, col_cap))
        if 100*(2*best_h-best_v) >= best_v:
            result.update(tied)
            # Una ejecución cuyas luces caben en una franja sigue en esa línea,
            # incluso si se permiten todos los empates del máximo global.
            if any(state & ~row == 0 for row in rows):
                assert all(v == 3 and row_cap == 3 for v, _, row_cap, _ in tied)
            if any(state & ~col == 0 for col in cols):
                assert all(v == 4 and col_cap == 4 for v, _, _, col_cap in tied)
    return result


def component_budgets(profiles):
    # Cinco unidades por paso lejano, una por paso cercano.
    ranks = []
    for width in (3, 4):
        row = [max(Q(5*(v-h)+1, h) for v, h, _, _ in profiles if width*h >= k*v)
               for k in range(1, width+1)]
        ranks.append(row)
    assert ranks == [[Q(11,3), 3, Q(1,3)], [Q(11,3), Q(11,3), 2, Q(1,3)]]
    assert list(map(sum, ranks)) == [Q(7), Q(29,3)]

    def pair(n1, n2):
        @lru_cache(None)
        def credit(r1, r2):
            if r1+r2 == 0:
                return 0
            best = None
            for v, h, row_cap, col_cap in profiles:
                if n1*h < r1*v or n2*h < r2*v:
                    continue
                cap1 = row_cap if n1 == 3 else col_cap
                cap2 = row_cap if n2 == 3 else col_cap
                minimum = (v+1)//2 if r1 == n1 or r2 == n2 else 1
                for a1 in range(min(r1,cap1,h)+1):
                    for a2 in range(min(r2,cap2,h-a1)+1):
                        if a1+a2 >= minimum:
                            value = 5*(2*h-v)-1+credit(r1-a1,r2-a2)
                            best = value if best is None else min(best, value)
            assert best is not None, (n1,n2,r1,r2)
            return best

        budgets = [[], []]
        max_overlap = min(n1,n2) if n1 == n2 else 1
        for overlap in range(max_overlap+1):
            for t1 in range(min(n1-overlap,n1//2)+1):
                for t2 in range(min(n2-overlap,n2//2)+1):
                    r1, r2 = n1-overlap-t1, n2-overlap-t2
                    balance = r1+r2-t1-t2
                    budgets[bool(overlap)].append(5*balance-credit(r1,r2))
        return tuple(map(max,budgets))

    pairs = [pair(3,3), pair(3,4), pair(4,4)]
    assert pairs == [(8,12), (8,12), (13,17)]
    return pairs


def rectangle(exp_upper, endpoint_caps, pairs):
    p, q, tau = map(Q, ("4e-10", ".001", ".01"))
    t = Q(11,15)**3
    y = t**5
    j = endpoint_caps(y)[-1]
    assert j == 142
    U, P, s = map(Q, ("7.7474732e-6", ".001288021", ".0012786947"))
    d1, d2 = map(Q, ("2.2689657e-6", "4.2725422e-5"))
    u1, u2 = map(Q, ("9.9106884e-9", "5.1859267e-8"))
    un, uf = map(Q, ("7.5080672e-7", "1.1673035e-6"))
    vn, vf, cf = map(Q, (".022344282", ".00022914915", "6.2413733e-8"))
    C1, C2 = map(Q, ("1.6716251e-9", "2.2905615e-9"))
    zp, zq = p/(1-p), q/(1-q)
    B, BD = 1297*U, 469*U+90*P
    l1, l2 = C1/t, C2/t
    x1, x2 = C1/t**7, C2/Q(11,15)**29

    R1, R2 = 640*d1+441*d2, 784*d1+513*d2
    V1 = (1+s)**3*exp_upper(1344*(un+uf)+30*vn)
    V2 = (1+s)**4*exp_upper(1596*(un+uf)+36*vn)
    paired = [(a*b/t**i, a*b/t**k) for a,b,(i,k) in
              ((C1,C1,pairs[0]), (C1,C2,pairs[1]), (C2,C2,pairs[2]))]
    (w11,o11),(w12,o12),(w22,o22) = paired
    pair1 = V1*(640*w11+9*o11)+V2*(441*w12+9*o12)
    pair2 = V1*(784*w12+16*o12)+V2*(513*w22+8*o22)
    child1 = V1*x1*(exp_upper(R1)-1)
    child2 = V2*x2*(exp_upper(R2)-1)
    large1 = exp_upper(R1)-1-R1+640*child1+441*child2
    large2 = exp_upper(R2)-1-R2+784*child1+513*child2
    L = j*vf
    G = (1+s)**12*exp_upper(1792*u1+1197*u2+3648*(un+uf)+L)
    Hn, Hf = 40*u1+27*u2, j*(un+cf)
    ratios = {
        "medida pasiva": zq*exp_upper(tau+2*BD)/P,
        "medida activa": zq*exp_upper(tau+180*P)/s,
        "peso de una VV": zp*exp_upper(tau+Q(17,3)*B+2*P*528)/C1,
        "peso de una CC": zp*exp_upper(tau+Q(103,15)*B+2*P*Q(3228,5))/C2,
        "rama VV de un grupo grande": V1*x1*exp_upper(R1)/d1,
        "rama CC de un grupo grande": V2*x2*exp_upper(R2)/d2,
        "entrada VV": V1*(l1+pair1+x1*large1)/u1,
        "entrada CC": V2*(l2+pair2+x2*large2)/u2,
        "cercano, raíz/puente/hijo": t*Hn*G/un,
        "lejano, raíz/puente": y*Hf*G/uf,
        "cercano desde su causa": t*(G-1)/vn,
        "lejano desde su causa": y*(G-1-L)/vf,
        "lejano desde su hijo": y*Hf*(G-1)/cf,
        "ancla de dato": max(u1+3*(un+uf),u2+4*(un+uf))/U,
    }
    largest = max(ratios.values())
    assert all(0 < r < Q(".998") for r in ratios.values()), ratios
    assert max(U,P) < 1 and p/Q("2.4e-11") == Q(50,3)
    print(f"Rectángulo (4e-10, 1e-3). {len(ratios)} cocientes racionales. "
          f"Máximo {float(largest):.9f} ({max(ratios,key=ratios.get)}).")


def tanner_ball(start, distance):
    degree = 3 if start == "V" else 4
    other = 4 if start == "V" else 3
    if distance % 2:
        return degree*(6**((distance+1)//2)-1)//5
    return 1+degree*(other-1)*(6**(distance//2)-1)//5


def pair_halos():
    # Dos coordenadas del centro y radios hasta un dato VV o CC de origen.
    shapes = {"medida": ("VC", ((8,9),(9,8))),
              "VV": ("VV", ((8,10),(9,9))),
              "CC": ("CC", ((9,9),(10,8)))}
    counts = {}
    for key, (center, radii) in shapes.items():
        origins = [tanner_ball(center[0],a)*tanner_ball(center[1],b) for a,b in radii]
        counts[key] = 1080*origins[0]+1296*origins[1]
    assert counts == {"medida":28198469520, "VV":63453609360, "CC":75207562128}
    return counts


def target(exp_upper):
    p, q, tilt = Q("1e-6"), Q(".001"), exp_upper(Q(".015"))
    b, m = p/(1-2*p)*tilt, q*tilt
    V, M = Q(".00000195"), Q(".00128")
    assert b*exp_upper(504*M)<V and m*exp_upper(178*M+938*V)<M
    marked_measurements = M/(1-178*M)
    caps = pair_halos()
    baseline_data = b*(exp_upper(504*M)-1)+2*1296*b*b
    baseline_measure = m*(exp_upper(178*M+938*V)-1)
    linear_data = 4*1296*22*252*b*b*M
    linear_measure = 4*caps["medida"]*b*b*marked_measurements
    triple_data = 6*b**3*(1296*22*1297+max(caps["VV"],caps["CC"]))
    data = baseline_data+linear_data+triple_data
    measure = baseline_measure+linear_measure
    assert data<Q("1.433e-6") and measure<Q(".000455")
    assert max(data,measure)+Q(".004545")<Q(".005")
    print(f"Objetivo (1e-6, 1e-3). Contribuciones controladas por raíz: "
          f"{float(data):.12g} y {float(measure):.12g}. Margen suficiente 0.004545.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para conservar las comprobaciones.")
    folder = Path(__file__).resolve().parent
    measure_proof = run_path(str(folder/"verificar_H5A_3_4_medidas_y_un_dato.py"))
    contact_proof = run_path(str(folder/"verificar_H5A_3_4_primer_contacto.py"))
    original = measure_proof["local_profiles"]()
    assert measure_proof["first_extinctions"](original) == (528, Q(3228,5))
    contact_proof["first_contact"](original)
    profiles = stripe_profiles()
    assert {(v,h) for v,h,_,_ in profiles} == {(v,h) for _,v,h,_ in original}
    prices = component_budgets(profiles)
    rectangle(measure_proof["exp_upper"], contact_proof["endpoint_caps"], prices)
    target(measure_proof["exp_upper"])
