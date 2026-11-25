"""Fallo silencioso del [[225,9,6]] en los extremos uniformes de la tabla 5."""
import argparse
import json
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial
from pathlib import Path
from runpy import run_path
from time import perf_counter


HERE = Path(__file__).resolve().parent
circuit = run_path(str(HERE / "circuito_HGP_225.py"))
geometry = run_path(str(HERE / "comparar_HGP_225.py"))
pilot = circuit["pilot"]
N, M = circuit["N"], circuit["M"]
TRIALS = 120_000
EXPONENT = 7
# p, q=r_m, inclinación circuital, r_c VV primero, r_c en todo orden admitido.
ROWS = (
    ("5.40e-7", "1e-3", ".015",  "1.005e-7", "9.18e-8"),
    ("1.00e-6", "1e-4", ".015",  "1.85e-7",  "1.68e-7"),
    ("1.06e-6", "1e-5", ".015",  "1.975e-7", "1.781e-7"),
    ("1.06e-6", "1e-5", ".0025", "2.00e-7",  "1.803e-7"),
)


def single_readout_triples(H, checks, tags, run):
    """Cuenta posibles fallos de tres datos y una lectura sobre todos los empates."""
    maps = geometry["affine_maps"](checks[:12])
    winners = circuit["all_winners"](H, checks, tags)
    readouts = [1 << j for j in range(M)]
    assert all(not winners(d) for d in readouts)
    classes = circuit["terminal_classes"](winners)

    total_triples = upper = representatives = 0
    for nvv in range(4):
        axis = geometry["axis_orbits"](nvv, maps)
        ordered = possible = 0
        for (first, size_i), (second, size_j) in product(axis, repeat=2):
            sites = [12*i+j if k < nvv else 144+9*i+j
                     for k, (i, j) in enumerate(zip(first, second))]
            if len(set(sites)) < 3:
                continue
            representatives += 1
            multiplicity = size_i * size_j
            ordered += multiplicity
            support = pilot["mask"](sites)
            s0 = pilot["xor_at"](checks, support)
            data_tag = pilot["xor_at"](tags, support)
            for d in readouts:
                s = s0 ^ d
                paths = classes(s, d)
                terminal, actual = run(s)
                if terminal == d:
                    assert actual in paths
                if any(tag != data_tag for tag in paths):
                    possible += multiplicity
        permutations = factorial(nvv) * factorial(3-nvv)
        assert ordered % permutations == possible % permutations == 0
        total_triples += ordered // permutations
        upper += possible // permutations
    assert representatives == 292 and total_triples == comb(N, 3)
    assert 0 < upper <= M*total_triples-134_073_792  # exclusión de la fase 7
    return upper


def bernoulli(p, q, pairs, triples_zero, triples_one_upper):
    w = [p**a * (1-p)**(N-a) for a in range(4)]
    m = [q**b * (1-q)**(M-b) for b in range(3)]
    lower = pairs*w[2]*(m[1]+99*m[2]) + triples_zero*w[3]*m[0]
    pair_upper = pairs*w[2]*q*(1-q)**8
    triple_upper = w[3]*(triples_zero*m[0] + triples_one_upper*m[1]
                           + comb(N, 3)*(1-m[0]-M*m[1]))
    tail = 1-sum(comb(N, a)*w[a] for a in range(4))
    return circuit["certificate"](lower, pair_upper+triple_upper+tail)


def circuit_point(rc, rm, layers, checks, tags, run, pairs, triples, seed):
    laws, norm = circuit["joint_fibres"](layers, rc)
    circuit["confirm_laws"](layers, laws, norm, rc)
    h = circuit["ancilla_rate"](rc, rm)
    spectra = circuit["pair_spectra"](laws, pairs, norm)
    circuit["contrast_pair"](laws, norm, pairs[0], h)
    lower_two, upper_two = circuit["pair_bounds"](spectra, h)
    clean_three = circuit["clean_triple_mass"](laws, norm, triples) * (1-h)**M
    model = circuit["trial_inputs"](laws, norm)
    burden = sum(model["prob_bg"]) + M*h
    lower = lower_two + clean_three
    upper_base = upper_two + clean_three + model["tail5"]
    upper_total = upper_base + model["mass"][3]*min(1, burden) + model["mass"][4]
    samples = circuit["sampled_sectors"](model, h, checks, tags, run, TRIALS, seed,
                                          tail_exponent=EXPONENT)
    upper_sample = min(upper_total, upper_base+sum(part[0] for part in samples.values()))
    return {
        "determinista": circuit["certificate"](lower, upper_total),
        "intervalo_conjunto": circuit["certificate"](lower, upper_sample),
        "estratos": {k: part[1] for k, part in samples.items()},
    }


def main():
    start = perf_counter()
    H, checks, tags, _, _, run = pilot["build"]()
    geometry["physical_actions"](H, checks, tags)
    geometry["verify_exterior_bit"]()
    pairs = []
    assert pilot["row_sectors"](checks, tags, run, patterns=pairs) == 324
    assert all(d and not ((checks[u] | checks[v]) & d) for _, u, v, d in pairs)
    triples, census = circuit["physical_triples"](H, checks, tags, run)
    triples_zero = sum(d == 0 for _, d in triples)
    assert triples_zero == 36 and len(triples) == 1656
    one_upper = single_readout_triples(H, checks, tags, run)
    orders = {}
    contracts = {}
    for name in ("VV primero", "CC primero"):
        data, ancillas = circuit["schedule"](checks, cc_first=name == "CC primero")
        orders[name] = data
        contracts[name] = circuit["primitive_checks"](H, checks, data, ancillas)
    assert orders["VV primero"] != orders["CC primero"]
    preparation = perf_counter()-start
    # 4 filas, 2 calendarios y 2 límites unilaterales. La suma parcial da e^7>800.
    risk_bounds = len(ROWS)*len(orders)*2
    exp_lower = sum(Q(EXPONENT)**j/factorial(j) for j in range(20))
    assert Q(risk_bounds)/exp_lower < Q(".02")
    data_points, circuit_points = {}, []
    for p, m, tilt, rv, rc in ROWS:
        if (p, m) not in data_points:
            result = bernoulli(Q(p), Q(m), len(pairs), triples_zero, one_upper)
            data_points[p, m] = {"p": p, "q": m, **result}
        for name, rate in zip(orders, (rv, rc)):
            before = perf_counter()
            seed = 2026112400 + 100*len(circuit_points)
            result = circuit_point(Q(rate), Q(m), orders[name], checks, tags,
                                   run, pairs, triples, seed)
            circuit_points.append({"orden_ejecutado": name, "rc": rate, "rm": m,
                                   "lambda_del_certificado": tilt, **result,
                                   "evaluacion_s": perf_counter()-before})
    return {"codigo": "[[225,9,6]]", "evento": "H_X R=0 y clase X no trivial",
            "decodificador": "SSF/FGL global, beta=.01, mayor |F| y menor mascara fisica",
            "modelos": "Bernoulli independiente; circuito primitivo de una ronda Z, sector X",
            "calendarios": contracts, "patrones_de_pareja": len(pairs),
            "censo_trios_circuito": census, "trios_Bernoulli_cero_lecturas": triples_zero,
            "trios_Bernoulli_una_lectura_max": one_upper,
            "Bernoulli": list(data_points.values()), "circuito": circuit_points,
            "confianza_simultanea_circuito_min": "0.98",
            "riesgo_union": f"{risk_bounds}*exp(-{EXPONENT})",
            "ensayos_por_estrato": TRIALS, "preparacion_s": preparation,
            "total_pared_s": perf_counter()-start, "pico_MiB": pilot["peak_memory"]()/2**20}


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar con las aserciones activas.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--salida", type=Path, help="Guardar los intervalos en JSON.")
    args = parser.parse_args()
    output = json.dumps(main(), indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if args.salida:
        args.salida.write_text(output, encoding="utf-8")
    else:
        print(output, end="")
