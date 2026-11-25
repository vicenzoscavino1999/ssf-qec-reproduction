"""Esquina Bernoulli q=10^-5 en la familia general HGP (3,4). Python estándar."""
from contextlib import redirect_stdout
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from io import StringIO
from pathlib import Path


def corner():
    path = Path(__file__).with_name("verificar_H5A_3_4_esquina_2.py")
    spec = spec_from_file_location("esquina_2_bernoulli", path)
    engine = module_from_spec(spec)
    spec.loader.exec_module(engine)
    engine.MEASUREMENT = Q("1e-5")
    engine.THETA = engine.B_PASS = Q("1.05e-5")
    engine.PAIR_PROB = tuple(map(Q, (".001324", ".00146", ".00462")))
    engine.configure()
    return engine


ALPHA = Q(".002")
ALPHA_V22 = Q(".004")
CASES = (
    ("1.06e-6", ".568", ".670",
     (("1.112736783e-6", "1.121456913e-6"), ("1.10389837e-6", "1.123909725e-6"),
      ("4.842158624e-8", "1.391992152e-7"), ("6.738693813e-9", "2.383012885e-8"),
      "6.696689095e-7", "2.496451587e-6")),
    ("2e-6", ".570", ".664",
     (("1.254499757e-6", "1.261214422e-6"), ("1.245749105e-6", "1.258357196e-6"),
      ("4.136930521e-8", "1.310891565e-7"), ("5.97044344e-9", "2.315898158e-8"),
      "3.882530107e-7", "7.855996598e-7")),
)


def data_activity(suite, data, p, witness):
    """Cuatro cocientes que dependen de p y comprobación explícita de R < alfa."""
    f5 = suite.F5
    exp_upper, dot = f5["exp_upper"], f5["dot"]
    m, u, _, physical = (tuple(map(Q, row)) for row in witness[:4])
    guard, other = suite.ENGINE["field_coefficients"]()
    halo = data[0]
    g = 40*physical[0] + 27*physical[1] + Q(witness[5])
    R = max(dot(row, u) for row in suite.ENGINE["FAR"]["N"])
    q, tilt, eta = suite.MEASUREMENT, f5["TILT"], f5["ETA"]
    assert 0 < p < Q(1, 2) and R < ALPHA and max(u) < eta
    y = q*exp_upper(tilt)/(1-2*q)
    z = p*exp_upper(tilt+ALPHA)/(1-2*p)
    init = ((210, 252), (640, 784), (441, 513), (1344, 1596))
    ratios = []
    for i in (0, 1):
        exposure = dot((2*u[0], 2*u[1], 2*g), [r[i] for r in halo[1:]])
        mass = z*exp_upper(eta + 2*R + guard*(504, 628)[i] + other*halo[0][i] + exposure)
        cost = dot((suite.B_PASS, 2*u[0], 2*u[1], 2*g), [r[i] for r in init])
        simple = z*exp_upper(eta + 2*R + cost)*(1+y)**(3+i)
        ratios.extend((mass/m[i], (simple+physical[i])/u[i]))
    # v22 fija alfa=.004 en la rama Bernoulli. Este argumento entrega al motor
    # exactamente el z racional de la tasa física. La geometría Bernoulli es independiente de p.
    fixed_exp = exp_upper(tilt + ALPHA_V22)
    argument = z/(fixed_exp+2*z)
    assert argument*fixed_exp/(1-2*argument) == z
    return argument, ratios, R


def five_source_test(suite, data, a, s, floor, p):
    """Minorante para 1 VV + 4 CC a los núcleos y suelos indicados, aun con alfa >= 0."""
    with redirect_stdout(StringIO()):
        baseline = suite.ENGINE["obstruction"](
            data, a, s, floor_fields=floor, compositions=((1, 4),),
            adjust_heads=suite.V29["adjust_heads"],
        )
    # El único monomio tiene grado cinco en p. Quitar toda la inclinación
    # histórica por dato da un minorante válido para cualquier alfa no negativo.
    return baseline*(p/(Q("1e-6")*suite.F5["exp_upper"](ALPHA_V22)))**5


def verify():
    suite = corner()
    suite.V12["check_local"]()
    data, circuit, pair = suite.coverage(*suite.geometry())
    guard, other = suite.ENGINE["field_coefficients"]()
    print(f"q={float(suite.MEASUREMENT):g}, guard={float(guard):.12g}, other={float(other):.12g}.", flush=True)
    accepted = set()
    for p_text, a, s, witness in CASES:
        p = Q(p_text)
        argument, direct, R = data_activity(suite, data, p, witness)
        name = f"Bernoulli p={p_text}, q=1e-5"
        case = (name, argument, None, a, s, witness)
        floor = suite.fields(witness)
        transcript = StringIO()
        try:
            with redirect_stdout(transcript):
                suite.ENGINE["verify_model"](case, data, circuit, floor_fields=floor,
                                             disjoint_pair_factors=pair, adjust_heads=suite.V29["adjust_heads"])
        except AssertionError as exc:
            # v22 adjunta el nombre y los cocientes fallidos sólo al acabar la fila completa.
            if len(exc.args) != 1 or not isinstance(exc.args[0], tuple) or len(exc.args[0]) != 2:
                raise
            label, failed = exc.args[0]
            if label != name or not failed or not all(0 <= i < 12 and r >= 1 for i, r in failed):
                raise
            closed = False
            print(f"{name}. Testigo fuera del presupuesto, máximo {max(r for _, r in failed):.12f}.")
        else:
            closed = True
            accepted.add(p)
            assert all(0 < r < 1 for r in direct)
            # La cabecera original imprime la tasa auxiliar y la etiqueta fija ".001".
            _, separator, details = transcript.getvalue().partition(". Fila ")
            assert separator
            print(f"{name}. Fila {details}", end="")
        head, tail = suite.ordinary_contributions(case, data, circuit, pair, floor)
        h = tuple(map(Q, witness[2]))
        if closed:
            assert all(x + t/budget < 1 for x, t, budget in zip(head, tail, h))
        else:
            assert min(head) > 1
        print(f"  R={float(R):.12g} < {float(ALPHA):g}, mayor cociente directo de datos={float(max(direct)):.12f}.")
        print("  Grados 1–10 / presupuesto V y C:", *(f"{float(x):.10f}" for x in head),
              "| colas >=11 V y C:", *(f"{float(x):.10g}" for x in tail), flush=True)
        necessary = five_source_test(suite, data, Q(a), Q(s), floor, Q("2e-6"))
        assert necessary > 1
        print(f"  A p=2e-6, (a,s)=({a},{s}), minorante de e c L para 1 VV + 4 CC: {float(necessary):.7f}.",
              flush=True)
    assert accepted == {Q("1.06e-6")}


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para conservar las comparaciones racionales.")
    verify()
