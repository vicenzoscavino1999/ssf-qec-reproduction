"""Esquina circuital r_m=10^-5 de la familia general HGP (3,4). Python estándar."""
from contextlib import redirect_stdout
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from io import StringIO
from pathlib import Path


def corner():
    path = Path(__file__).with_name("verificar_H5A_3_4_esquina_2.py")
    spec = spec_from_file_location("esquina_2_para_circuito", path)
    suite = module_from_spec(spec)
    spec.loader.exec_module(suite)
    suite.MEASUREMENT = Q("1e-5")
    suite.THETA = suite.B_PASS = Q("1.15e-5")
    suite.PAIR_PROB = tuple(map(Q, (".00145", ".00160", ".00506")))
    suite.configure()
    suite.F5.update(ALPHA=Q(".0015"), ETA=Q(".000002"))
    return suite


def verify_row(suite, data, circuit, pair, tilt, case):
    """Comprueba la fila original y separa el tramo 1–10 de la cola infinita."""
    name, rate, _, _, _, witness = case
    suite.F5["TILT"] = Q(tilt)
    floor = suite.fields(witness)
    transcript = StringIO()
    try:
        with redirect_stdout(transcript):
            suite.ENGINE["verify_model"](
                case, data, circuit, floor_fields=floor, disjoint_pair_factors=pair,
                adjust_heads=suite.V29["adjust_heads"],
            )
    except AssertionError as exc:
        # v22 sólo adjunta el nombre de la fila a la última comparación.
        # Un fallo en las hipótesis o en la construcción de la cola se propaga.
        diagnostic = exc.args[0] if len(exc.args) == 1 else None
        if not isinstance(diagnostic, tuple) or len(diagnostic) != 2 or diagnostic[0] != name:
            raise
        failures = diagnostic[1]
        if not failures or not all(0 <= j < 12 and r >= 1 for j, r in failures):
            raise
        closed = False
        print(f"{name}, rc={rate}, rm={float(suite.MEASUREMENT):g}, λ={tilt}. "
              f"Testigo abierto, máximo {max(r for _, r in failures):.12f}.", flush=True)
    else:
        closed = True
        # v22 imprime una etiqueta fija de medición de su esquina original.
        print(transcript.getvalue().replace("medición=.001.", f"rm={float(suite.MEASUREMENT):g}, λ={tilt}."),
              end="", flush=True)

    head, tail = suite.ordinary_contributions(case, data, circuit, pair, floor)
    if closed:
        assert all(g + t/Q(h) < 1 for g, t, h in zip(head, tail, witness[2]))
    else:
        assert min(head) > 1  # En estos testigos, los diez primeros tamaños ya exceden ambas filas.
    print("  Grados 1–10 / h (VV, CC):", *(f"{float(g):.10f}" for g in head),
          "| cola absoluta ≥11:", *(f"{float(t):.10g}" for t in tail), flush=True)
    return closed


# λ, cierre esperado, r_c, cota de orden VV (CC usa 7), a, s y testigo de v22.
CASES = (
    (".015", True, "1.975e-7", 4, "0.558", "0.745",
     (("1.572938652e-6", "1.557254298e-6"), ("7.166414931e-7", "1.273103295e-6"),
      ("6.702895968e-8", "2.891927541e-7"), ("9.126381283e-9", "4.537049656e-8"),
      "9.099802684e-7", "1.804826158e-6")),
    (".015", True, "1.781e-7", 7, "0.5553", "0.761",
     (("1.420877124e-6", "1.407361776e-6"), ("1.01973753e-6", "1.157628947e-6"),
      ("1.064456029e-7", "3.047018754e-7"), ("1.267429447e-8", "4.925669331e-8"),
      "9.902643515e-7", "1.746589069e-6")),
    (".015", False, "2e-7", 4, "0.558", "0.745",
     (("1.576951538e-6", "1.561154514e-6"), ("7.183968169e-7", "1.27592756e-6"),
      ("6.685048701e-8", "2.893117791e-7"), ("9.109888147e-9", "4.543051088e-8"),
      "8.985674849e-7", "1.771066213e-6")),
    (".015", False, "2e-7", 7, "0.5562", "0.755",
     (("1.454765649e-6", "1.440387121e-6"), ("1.042981477e-6", "1.179342145e-6"),
      ("1.018621205e-7", "2.940175078e-7"), ("1.231036712e-8", "4.814385266e-8"),
      "9.004478329e-7", "1.522259211e-6")),
    (".0025", True, "2e-7", 4, "0.558", "0.745",
     (("1.57296384e-6", "1.557278799e-6"), ("7.166522521e-7", "1.273120517e-6"),
      ("6.702790554e-8", "2.891939147e-7"), ("9.126296679e-9", "4.537100856e-8"),
      "9.099108922e-7", "1.804617982e-6")),
    (".0025", True, "1.803e-7", 7, "0.5553", "0.761",
     (("1.420812423e-6", "1.407298773e-6"), ("1.019692416e-6", "1.157583846e-6"),
      ("1.064505175e-7", "3.046979088e-7"), ("1.267452183e-8", "4.925517312e-8"),
      "9.904861856e-7", "1.747171324e-6")),
    (".0025", False, "2e-7", 7, "0.556", "0.757",
     (("1.451026085e-6", "1.436717924e-6"), ("1.040434122e-6", "1.17701135e-6"),
      ("1.025926252e-7", "2.966284278e-7"), ("1.234819963e-8", "4.834940769e-8"),
      "9.098572413e-7", "1.532735321e-6")),
)


def main():
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para mantener las comprobaciones racionales.")
    suite = corner()
    suite.V12["check_local"]()
    data, circuit, pair = suite.coverage(*suite.geometry())
    for tilt, expected, rate, first, a, s, witness in CASES:
        name = "Circuito VV primero" if first == 4 else "Circuito en cualquier orden admitido"
        case = (name, rate, first, a, s, witness)
        closed = verify_row(suite, data, circuit, pair, tilt, case)
        assert closed == expected, (name, rate, tilt, closed)


if __name__ == "__main__":
    main()
