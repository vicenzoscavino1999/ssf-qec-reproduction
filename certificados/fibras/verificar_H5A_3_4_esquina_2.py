"""Esquina de medida 10^-4 del certificado general HGP (3,4). Python estándar."""
from contextlib import redirect_stdout
from fractions import Fraction as Q
from io import StringIO
from pathlib import Path
from runpy import run_path


HERE = Path(__file__).resolve().parent
V29 = run_path(str(HERE / "verificar_H5A_3_4_primer_apagado_de_ternas.py"))
ENGINE = V29["engine"]
F5 = ENGINE["F5"]
V12 = ENGINE["V15"]["v12"]
MEASUREMENT = Q(".0001")
THETA = Q(".000108")
B_PASS = Q(".0001085")
PAIR_PROB = tuple(map(Q, (".013611", ".015012", ".04752")))
original_pair_prices = ENGINE["pair_prices"]


def physical_prices(v15, a, theta, ball, discounts, hist):
    """v18 con el techo CC-CC ya comprobado en v12 para esta esquina."""
    y = a**6
    ordinary, marked = v15["v12"]["prices"](y, Q(1), theta, ball, discounts, hist)
    clean = y*y + PAIR_PROB[2]
    clean_one = y*y + PAIR_PROB[2]/y
    pair = 384*clean + 72*clean_one + 24 + 20/y + 12/y**2
    ordinary[1][0, 2] = min(ordinary[1][0, 2], discounts[0, 2]*pair)
    return ordinary, marked


def pair_prices(i, j, energies, y, discount, regular, mark, circuit):
    """La fila Bernoulli de v22 es afín en el techo local de incidencias."""
    old = original_pair_prices(i, j, energies, y, discount, regular, mark, circuit)
    if circuit:
        return old
    sphere = ENGINE["sphere"]

    def count(xy):
        return sphere(i, xy[0])*sphere(i, xy[1])

    if i == j == 0:
        slope = sum(count(p)*energies[p] for p in energies if p[0] != 0)
    elif i == j == 1:
        slope = 384*energies[4, 2] + sum(count(p)*energies[p]/y for p in ((0, 2), (2, 2)))
    else:
        slope = count((3, 3))*energies[3, 3]*(1 + 1/(2*y))
    previous_cap = (Q(".152"), Q(".167"), Q(".528"))[i+j]
    result = old - discount*slope*(previous_cap - PAIR_PROB[i+j])
    assert 0 < result <= old
    return result


def configure():
    # run_path carga una instancia privada. Todas las funciones activas usan este F5.
    F5.update(RM=MEASUREMENT, THETA=THETA, B_PASS=B_PASS)
    assert F5["TILT"] == Q(".015")
    local = V12["prices"].__globals__
    assert local is V12["check_local"].__globals__
    local.update(THETA_MAX=THETA, PAIR_PROB=PAIR_PROB)
    ENGINE["v18"]["physical_prices"] = physical_prices
    ENGINE["verify_model"].__globals__["pair_prices"] = pair_prices
    guard, other = ENGINE["field_coefficients"]()
    assert 0 < guard < Q(".001") and other > 0


def geometry():
    ENGINE["check_local"]()
    data, circuit = ENGINE["v21"]["previous_tables"]()
    assert data[2] == circuit[3]
    last = V29["last"]
    data, circuit = last["improved_tables"](data, circuit, last["check_local"](), False)
    caps = V29["geometry_caps"]()
    assert caps == {((1, 0, 3), 4): 432, ((0, 1, 2), 3): 5022, ((1, 1, 2), 3): 17904}
    return V29["tighten"](data, circuit, caps)


def coverage(data, circuit):
    # El techo trivial 1 es válido en las dos esquinas. El valor final se
    # calcula con el guard nuevo mediante el Taylor de v27, incluidos los pares.
    trivial = dict.fromkeys(data[2], Q(1))
    discount, pair = V29["last"]["cross"]["coverage_discounts"](ENGINE, circuit[4], trivial)
    assert set(discount) == set(data[2]) == set(circuit[3])
    assert all(0 < value <= 1 for value in (*discount.values(), *pair.values()))
    return (*data[:2], discount, data[3]), (*circuit[:3], discount, *circuit[4:]), pair


def fields(witness):
    _, u, h, physical, visible, root = witness
    u, h, physical = (tuple(map(Q, row)) for row in (u, h, physical))
    g = 40*physical[0] + 27*physical[1] + Q(root)
    return (*ENGINE["field_coefficients"](), 2*u[0], 2*u[1], 2*g, *h, Q(visible))


CASES = (
    ("Bernoulli", "1e-6", None, ".568", ".678",
     (("1.103585668e-6", "1.125499296e-6"), ("1.064424176e-6", "1.087972930e-6"),
      ("4.753509019e-8", "1.387925312e-7"), ("6.387261997e-9", "2.275491635e-8"),
      "6.862867109e-7", "2.337412703e-6")),
    ("Circuito VV primero", "2e-7", 4, ".56", ".735",
     (("1.583803264e-6", "1.586540506e-6"), ("7.009445854e-7", "1.246307714e-6"),
      ("6.208982359e-8", "2.636955071e-7"), ("8.480201465e-9", "4.121954550e-8"),
      "8.863421353e-7", "1.732494362e-6")),
    ("Circuito en cualquier orden admitido", "2e-7", 7, ".558", ".745",
     (("1.462315182e-6", "1.465048528e-6"), ("1.018676803e-6", "1.153293948e-6"),
      ("9.512959359e-8", "2.707254298e-7"), ("1.150919716e-8", "4.441037537e-8"),
      "8.825105554e-7", "1.483819788e-6")),
    ("Circuito VV primero", "1.85e-7", 4, ".56", ".735",
     (("1.556645943e-6", "1.559762524e-6"), ("6.892972725e-7", "1.227067616e-6"),
      ("6.292080239e-8", "2.623605436e-7"), ("8.536645986e-9", "4.072267627e-8"),
      "9.521317178e-7", "1.937809829e-6")),
    ("Circuito en cualquier orden admitido", "1.68e-7", 7, ".556", ".757",
     (("1.410698352e-6", "1.414168190e-6"), ("9.843092660e-7", "1.121928294e-6"),
      ("1.029704577e-7", "2.944270920e-7"), ("1.225033433e-8", "4.757120021e-8"),
      "9.961419869e-7", "1.788974783e-6")),
)


def ordinary_contributions(case, data, circuit, pair, floor):
    """Separar los grados 1–10 y >=11 en las dos filas ordinarias de v22."""
    _, rate, first, a, s, witness = case
    rate, a, s = map(Q, (rate, a, s))
    m, u, h, physical = (tuple(map(Q, row)) for row in witness[:4])
    visible, root = map(Q, witness[4:])
    y = a**6
    _, discount, group, mark, law = ENGINE["ordinary_heads"](data, circuit, rate, first, a)
    V29["adjust_heads"](group, mark, discount, data, circuit, rate, first, a, s, floor)
    raw, raw_mark = ([row.copy() for row in table] for table in (group, mark))
    v21, ball = ENGINE["v21"], ENGINE["BALL"]
    v21["pair_mark"](mark, ENGINE["ANCHOR"], ball, discount, law is not None)
    v21["joint_rows"](group, mark, B_PASS, THETA, y, s, u, h, physical, visible, root)
    costs = ENGINE["rounded_prices"](ENGINE["source_prices"](a, rate if law is not None else None, first))
    matrix = ENGINE["matrices"](a, s, floor)
    polynomial = ENGINE["polynomials"](costs, matrix, 13)
    actual, children = fields(witness), (*h, visible)
    dot = F5["dot"]
    for i in (0, 1):
        for key in group[i]:
            group[i][key] = min(group[i][key], polynomial[i][key]*discount.get(key, 1))
        forbidden = dot(v21["DIRECT"][i], children)
        for j in (0, 1):
            left, right = sorted((i, j))
            factors = {}
            for x, y_pos in v21["displacements"](left, right):
                pos = (left, right, x, y_pos)
                loss = forbidden + dot(v21["initial_overlap"](*pos), actual[:5])
                loss += sum(k*v21["link_overlap"](*pos, t, y, s) for t, k in enumerate(children))
                factor = v21["upper_factor"](loss)
                if ENGINE["edge_kind"](*pos) != 2:
                    factor *= pair[left, right]
                factors[x, y_pos] = factor
            key = (2-i-j, i+j)
            bound = pair_prices(i, j, factors, y, discount[key], raw[i][key],
                                raw_mark[i].get(key, 0), law is not None)
            group[i][key] = min(group[i][key], bound)
    dd, dg, _, _ = ENGINE["FAR"]["kernels"](ball, y, s)
    w = [F5["rounded"](m[i]*F5["exp_upper"](dot([k[0] for k in dd[i]], h) + dg[i][0]*visible))
         for i in (0, 1)]
    full = ENGINE["total"](w, costs, matrix)
    partial = ENGINE["partial"]
    tail = [full[i] - partial(polynomial[i], w, 1, 10) for i in (0, 1)]
    assert all(t >= partial(polynomial[i], m, 11, 13) >= 0 for i, t in enumerate(tail))
    return [F5["evaluate"](group[i], w)/h[i] for i in (0, 1)], tail


def run_rows(data, circuit):
    V12["check_local"]()
    data, circuit, pair = coverage(data, circuit)
    guard, other = ENGINE["field_coefficients"]()
    print(f"Medición {float(MEASUREMENT):g}, guard {float(guard):.12g}, other {float(other):.12g}.", flush=True)
    certified = set()
    for case in CASES:
        name, rate, first, _, _, witness = case
        floor = fields(witness)
        transcript = StringIO()
        try:
            with redirect_stdout(transcript):
                ENGINE["verify_model"](case, data, circuit, floor_fields=floor,
                                       disjoint_pair_factors=pair, adjust_heads=V29["adjust_heads"])
        except AssertionError as exc:
            # Únicamente se admite el diagnóstico final de los doce cocientes
            # de v22. Una excepción anterior, incluso en la cola, es un error.
            if len(exc.args) != 1 or not isinstance(exc.args[0], tuple):
                raise
            label, failed = exc.args[0]
            if label != name or not failed or not all(0 <= j < 12 and r >= 1 for j, r in failed):
                raise
            closed = False
            print(f"{name}, tasa {rate}, medición {float(MEASUREMENT):g}. "
                  f"Testigo fuera del presupuesto, mayor cociente {max(r for _, r in failed):.12f}.",
                  flush=True)
        else:
            closed = True
            certified.add((first, Q(rate)))
            # La salida histórica de v22 contiene la etiqueta literal de la esquina 1.
            print(transcript.getvalue().replace("medición=.001.", f"medición={float(MEASUREMENT):g}."),
                  end="", flush=True)
        head, tail = ordinary_contributions(case, data, circuit, pair, floor)
        h = tuple(map(Q, witness[2]))
        if closed:
            assert all(x + t/hi < 1 for x, t, hi in zip(head, tail, h))
        else:
            assert min(head) > 1
        print("  Grados 1–10 / presupuesto V y C:", *(f"{float(x):.10f}" for x in head),
              "| cola de 11 a infinito V y C:", *(f"{float(t):.10g}" for t in tail), flush=True)
    assert certified == {(None, Q("1e-6")), (4, Q("1.85e-7")), (7, Q("1.68e-7"))}


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para mantener las comparaciones racionales.")
    configure()
    run_rows(*geometry())

