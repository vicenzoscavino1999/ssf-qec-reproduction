"""Compare the 225-qubit schedules at identical rates, using the archived circuit counts."""
import argparse
import json
from fractions import Fraction
from pathlib import Path
from runpy import run_path


def compare(release):
    folder = release / "comparaciones"
    pilot = run_path(str(folder / "esquinas_HGP_225.py"))
    archived = json.loads((folder / "piloto_HGP_225_esquinas_20261124.json").read_text(encoding="utf-8"))
    circuit = pilot["circuit"]
    seed_matrix, checks, tags, _, _, decode = pilot["pilot"]["build"]()
    pairs = []
    assert pilot["pilot"]["row_sectors"](checks, tags, decode, patterns=pairs) == 324
    triples, _ = circuit["physical_triples"](seed_matrix, checks, tags, decode)
    assert len(triples) == 1656
    schedule, _ = circuit["schedule"](checks, cc_first=False)
    rows = []
    for i, (_, rm, _, _, rc) in enumerate(pilot["ROWS"]):
        cc = archived["circuito"][2*i+1]
        assert (cc["orden_ejecutado"], cc["rc"], cc["rm"]) == ("CC primero", rc, rm)
        seed = 2026112500 + 100*i
        vv = pilot["circuit_point"](Fraction(rc), Fraction(rm), schedule, checks, tags,
                                    decode, pairs, triples, seed)
        v, c = vv["intervalo_conjunto"], cc["intervalo_conjunto"]
        lower = Fraction(c["L"]) / Fraction(v["U"])
        upper = Fraction(c["U"]) / Fraction(v["L"])
        point = Fraction(c["estimacion"]) / Fraction(v["estimacion"])
        assert 1 < lower <= point <= upper
        fmt = circuit["decimal_bound"]
        rows.append({"rc": rc, "rm": rm, "first_VV": v, "first_CC": c,
                     "first_VV_deterministic": vv["determinista"],
                     "ratio_lower": fmt(lower, circuit["ROUND_FLOOR"]),
                     "ratio_upper": fmt(upper, circuit["ROUND_CEILING"]),
                     "ratio_point": fmt(point, circuit["ROUND_HALF_EVEN"]),
                     "first_VV_strata": vv["estratos"], "first_CC_strata": cc["estratos"]})
    return {"trials_per_stratum": pilot["TRIALS"], "exponent": pilot["EXPONENT"],
            "joint_noncoverage_bound": "16*exp(-7)", "rows": rows}


if __name__ == "__main__":
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("release", type=Path, help="Root of the reproduction package")
    print(json.dumps(compare(args.parse_args().release), indent=2, ensure_ascii=True))
