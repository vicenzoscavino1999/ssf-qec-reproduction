"""Coste de la fila VV solicitada con el productor y las semillas de la tabla 6."""
import argparse
import json
import statistics
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
from runpy import run_path
from time import perf_counter

HERE = Path(__file__).resolve().parent
INPUT = HERE.parent

def one_point():
    start = perf_counter()
    producer = run_path(str(INPUT / "esquinas_HGP_225.py"))
    circuit = producer["circuit"]
    pilot = producer["pilot"]
    H, checks, tags, _, _, run = pilot["build"]()
    producer["geometry"]["physical_actions"](H, checks, tags)
    producer["geometry"]["verify_exterior_bit"]()
    pairs = []
    assert pilot["row_sectors"](checks, tags, run, patterns=pairs) == 324
    triples, _ = circuit["physical_triples"](H, checks, tags, run)
    data, anc = circuit["schedule"](checks, cc_first=False)
    contract = circuit["primitive_checks"](H, checks, data, anc)
    assert len(triples) == 1656 and contract["N_loc"] == 3240
    prepared = perf_counter()
    result = producer["circuit_point"](Fraction("2e-7"), Fraction("1e-5"),
                                      data, checks, tags, run, pairs, triples, 2026113000)
    evaluated = perf_counter()
    ref = json.loads((INPUT / "piloto_HGP_225_esquinas_20261124.json").read_bytes())
    row, = [v for v in ref["circuito"] if v["orden_ejecutado"] == "VV primero"
            and Fraction(v["rc"]) == Fraction("2e-7") and Fraction(v["rm"]) == Fraction("1e-5")]
    assert result["intervalo_conjunto"] == row["intervalo_conjunto"]
    assert result["determinista"] == row["determinista"]
    assert {str(k): v for k, v in result["estratos"].items()} == row["estratos"]
    return {**result, "preparacion_s": prepared-start, "evaluacion_s": evaluated-prepared,
            "pico_MiB": pilot["peak_memory"]()/2**20}

def comparison():
    runs = []
    for _ in range(3):
        start = perf_counter()
        proc = subprocess.run([sys.executable, "-I", "-B", "-X", "utf8", __file__, "--un-punto"],
                              check=True, capture_output=True)
        elapsed = perf_counter()-start
        trial = json.loads(proc.stdout)
        trial["pared_s"] = elapsed
        runs.append(trial)
        print(f"Certificado VV, proceso {len(runs)}/3: {elapsed:.3f} s", flush=True)
    return {"rc": "2e-7", "rm": "1e-5", "orden": "VV primero", "corridas": runs,
            "mediana_pared_s": statistics.median(r["pared_s"] for r in runs)}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--un-punto", action="store_true")
    args = parser.parse_args()
    if args.un_punto:
        print(json.dumps(one_point()))
    else:
        output = comparison()
        (HERE / "certificado_225_tiempos.json").write_text(json.dumps(output, indent=2), encoding="utf-8")
