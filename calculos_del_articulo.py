"""Magnitudes de la versión limpia que utilizan varios de los suplementos."""
from collections import Counter
from fractions import Fraction as F
from math import comb
from pathlib import Path
import csv
import json
import runpy
import sys

ROOT = Path(__file__).resolve().parent
SECTORES = ROOT / "sectores"
COMPARACIONES = ROOT / "comparaciones"
WORK = ROOT / "tiras/workspace"

def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))

def exp_lower(x):
    # e^x = (e^(x/128))^128. Taylor con términos positivos y redondeo
    # hacia abajo antes de la potencia final mantienen la cota racional.
    term = total = F(1)
    for j in range(1, 81):
        term *= x / (128*j)
        total += term
    scale = 10**32
    return F(int(total*scale), scale)**128

def finite_codes():
    sys.path.insert(0, str(SECTORES))
    from intervalo import certificate, regional_interval
    p, q = F("1e-6"), F("1e-4")
    results = {}
    for n, m, c, cuts, lower, upper, static_limit in [
        (225, 108, 324, (6482,3028), "3.24009e-14", "3.27319e-14", "5.21e-12"),
        (900, 432, 864, (6575,3057), "8.63039e-14", "9.69644e-14", "8.42e-9")
    ]:
        data = read(SECTORES / "referencia" / f"resultado_{n}.json")
        stored = data["interval"]
        counts = {tuple(map(int,k.split(","))): v for k,v in stored["positive_census"].items()}
        assert counts[2,1] == c
        fresh = certificate(data["geometry"], counts, stored["first_step_exceptional_triples"])
        assert fresh == stored
        assert F(lower) <= F(fresh["lower"]) <= F(fresh["upper"]) <= F(upper)
        regional = regional_interval(n, p, q)
        t21 = c*p**2*q*(1-p)**(n-2)*(1-q)**(m-1)
        if n == 225:
            pilot = read(COMPARACIONES / "piloto_HGP_225_esquinas_20261124.json")
            cap = next(F(row["U"]) for row in pilot["Bernoulli"] if (F(row["p"]), F(row["q"])) == (p, q))
        else:
            cap = F(fresh["upper"])
        rho, share, target = (F(".0180"), F("98.23"), F(".90")) if n == 225 else (F(".1728"), F("85.26"), F(8))
        assert F(regional["lower"]) <= t21 <= cap <= (1+rho)*t21 <= F(regional["upper"])
        assert 100/(1+rho) >= share and F(regional["max_relative_percent"]) < target
        witness = 1/(1+((1-p)/p)**6)
        sector = c*p*p*q*(1-(n-2)*p-(m-1)*q)
        for order, mass in zip(cuts, [witness,sector]):
            assert exp_lower(F(".015")*(order+1))*mass > 5120*(n+m)
        assert comb(n,6)*p**6 / F(lower) < F(static_limit)
        results[str(n)] = {"P_sil": [fresh["lower"],fresh["upper"]], "C21": c, "regional": regional,
                          "cortes_suficientes": cuts, "fraccion_sin_flip_menor_que": static_limit}
    return results

def opening_figure():
    """Replay both Figure 1 inputs with the paper's complete FGL family."""
    from itertools import combinations

    audit = runpy.run_path(str(COMPARACIONES / "auditar_esquinas_HGP_225.py"))
    _, checks, generators = audit["incidence"]()
    tags, decode, _ = audit["decoder"](checks, generators)
    xor_sum, mask, sites = audit["bit_sum"], audit["mask"], audit["positions"]

    candidates = set()
    for generator in generators:
        for weight in range(1, len(generator) + 1):
            for qubits in combinations(generator, weight):
                footprint = xor_sum(checks[v] for v in qubits)
                if 2 * footprint.bit_count() >= 3 * weight:
                    candidates.add((footprint, weight, mask(qubits)))

    # Rows v0..v2 intersect columns v3..v5 at the checks drawn in panel (a).
    grid = [[list(sites(checks[row] & checks[col])) for col in (3, 4, 5)]
            for row in (0, 1, 2)]
    assert grid == [[[0], [6], [3]], [[4], [1], [7]], [[8], [5], [2]]]
    data = mask((3, 4))
    physical = checks[3] ^ checks[4]
    assert list(sites(physical)) == [0, 1, 4, 5, 6, 8]
    tex = next((ROOT / "manuscrito").glob("*.tex")).read_text(encoding="utf-8")
    end = tex.index(r"\label{fig:affine-trace}")
    figure = tex[tex.rindex(r"\begin{figure}", 0, end):end]
    results = {}
    for label, discrepancy, expected in (("correct", 0, [3, 4]), ("inverted", 1 << 2, [2, 0, 1, 5])):
        observed = physical ^ discrepancy
        residual, syndrome = data, observed
        flips = []
        while True:
            # Global maximum score, then largest action, then smallest physical mask.
            score, _, neg_action, footprint = max(
                (F(2 * (syndrome & phi).bit_count(), phi.bit_count()) - 1, weight, -action, phi)
                for phi, weight, action in candidates
            )
            if score < F(1, 100):
                break
            action = -neg_action
            flips.append(list(sites(action)))
            residual ^= action
            syndrome ^= footprint

        assert flips == [[v] for v in expected]
        assert "$" + r"\to ".join(f"v_{v}" for v in expected) + "$" in figure
        assert list(sites(residual)) == (list(range(6)) if discrepancy else [])
        assert syndrome == discrepancy and xor_sum(checks[v] for v in sites(residual)) == 0
        logical_tag = xor_sum(tags[v] for v in sites(residual))
        assert bool(logical_tag) == bool(discrepancy)
        assert decode(observed) == (syndrome, xor_sum(tags[v] for v in sites(data ^ residual)))
        results[label] = {"flips": flips, "residual": list(sites(residual)),
                          "recorded_checks": list(sites(syndrome)), "logical_X": bool(logical_tag)}
    return results

def affine_17():
    n, m = 177157, 88434
    p, q = F("5e-33"), F("2e-11")
    zp = p/(1-p)
    pair = zp**34/(1+zp**34)
    classes = 5202*p**17*(1-(n-17)*p-m*q)
    R = n*F("2e-16")+m*F("2.057e-11")
    assert pair > F("5.82e-1099") and classes > F("3.96e-546") and R < F("1.82e-6")
    for order, mass in [(253000,pair),(126000,classes)]:
        assert exp_lower(F(".01")*(order+1))*mass > 5120*R
    # D.3: mismo polinomio geométrico, tasas actuales del circuito afín.
    geom = runpy.run_path(str(ROOT / "certificados/sector_f17/verify_budget.py"))
    bounds, poly = geom["bounds"], geom["polynomial"]
    rc, rm = F("5e-33"), F("2e-11")
    pw, pa = 2*rc/3, 8*rc/15
    zw, zc = pw/(1-pw), pa/(1-pa)
    theta = (1-(1-2*pw)**19*((15-16*rc)/(15-8*rc))**35*(1-2*rm))/2
    tl, tu = bounds(theta)
    def rho(d):
        return bounds((d*zw*theta+zc/2)/(d*zw+zc))
    rl,ru = rho(17)
    sl,su = rho(2)
    w0 = 1-6471288*pw-3095190*pa
    ordinary = 1-(m-18)*tu
    assert 0 < w0 <= 1 and 0 < ordinary <= 1 and 0 < rl <= ru < 1 and 0 < sl <= su < 1
    lower = (306*w0*(17*zw+zc)**14*(2*zw+zc)*tl**3*ordinary
             *(1-ru)**14*(1-su)*poly(rl/(1-rl),sl/(1-sl)))
    assert lower > F("1.02e-499")
    return {"cortes_suficientes":[253000,126000], "L_fin_mayor_que":"1.02e-499"}

def comparisons():
    splitting = COMPARACIONES / "SPLITTING_HGP_225_VV_20261124"
    with (splitting/"cadenas.tsv").open(encoding="utf-8",newline="") as stream:
        final = [r for r in csv.DictReader(stream, delimiter="\t") if int(r["j"])==27]
    moves = Counter(int(r["moves"]) for r in final)
    assert len({(r["seed"],r["rep"]) for r in final}) == 36
    assert {int(r["n"]) for r in final} == {100000}
    assert moves == {0:26, 2:8, 4:1, 6:1}
    run = read(splitting/"resultado_splitting.json")
    rates = [tuple(map(F,level["tasas"])) for level in run["diagnosticos"]]
    assert len(rates)==28 and rates[0]==(F(".003"),F(".005")) and rates[-1]==(F("2e-7"),F("1e-5"))
    assert all(rm==min(F(".005"),50*rc) for rc,rm in rates)
    run61 = read(WORK / "T2_2_PILOTO_M6_EJECUCION_20261115/corrida/RESULTADO_T2_2.json")
    assert run61["escalera"]["niveles"]==25 and run61["escalera"]["indices_de_los_extremos"]["1/100000"]==18
    return {"225_niveles":28, "225_transiciones_finales":sum(k*v for k,v in moves.items()),
            "225_cadenas_moviles":36-moves[0], "225_cadenas_inmoviles":moves[0], "61_niveles":25}

def plot_coordinates():
    """Read the plotted numbers in the TeX and compare at their printed precision."""
    import re
    tex = next((ROOT/"manuscrito").glob("*.tex")).read_text(encoding="utf-8")
    checked = 0
    def number(printed, value):
        nonlocal checked
        unit = F(1,10**len(printed.partition(".")[2]))
        assert abs(F(printed)-value) < unit, (printed, float(value))
        checked += 1
    matched = read(COMPARACIONES/"postselection_v2/matched_schedules.json")["rows"]
    tables = re.findall(r"row\s+est\s+lo\s+hi\s*([\d.\s]+)\}", tex)
    vv, cc, sizes = [t.strip().splitlines() for t in tables]
    for group, schedule in ((vv,"first_VV"),(cc,"first_CC")):
        for idx,(row,record) in enumerate(zip(group,matched,strict=True)):
            tag, *values = row.split()
            assert int(tag)==4-idx
            limits = record[schedule]
            for printed,key in zip(values,("estimacion","L","U"),strict=True):
                number(printed,10**15*F(limits[key]))
    for row,n in zip(sizes,(225,900),strict=True):
        _, *values = row.split()
        data = read(SECTORES/"referencia"/f"resultado_{n}.json")["interval"]
        for printed,key in zip(values,("estimate","lower","upper"),strict=True):
            number(printed,10**14*F(data[key]))
    def figure(label):
        stop = tex.index(r"\label{"+label+"}")
        return tex[tex.rindex(r"\begin{figure",0,stop):stop]
    split = read(COMPARACIONES/"SPLITTING_HGP_225_VV_20261124/resultado_splitting.json")
    section = tex.split(r"\label{sec:splitting-comparison}", 1)[1].split(r"\section{", 1)[0]
    printed = re.search(r"splitting estimate of \$([\d.]+)\\cdot10\^\{-16\}", section)[1]
    number(printed, 10**16*F(str(split["estimacion"])))
    strips = read(ROOT/"tiras/tables/results/circuit.json")["instancias"]
    def cell(length, power):
        rate = f"1/{10**(-power)}"
        d = strips[f"TIRA_{length}x6"]["celdas"][rate+"|"+rate]["por_k"]["4"]
        lo,hi = map(F,d["intervalo_P_sil_exacto"])
        tail,mass = F(d["eps_ing"]),F(d["L_ing"])
        return lo,hi,100*(hi-lo)/(hi+lo),100*tail/(2*mass+tail)
    f5 = figure("fig:strip-comparison")
    for length,n in zip((6,12,24,48),(61,127,259,523)):
        global_value, component = re.findall(r"\("+str(n)+r",([\d.]+)\)",f5)
        lo,hi,comp,glob = cell(length,-5)
        number(global_value,glob)
        number(component,comp)
        if n==61:
            for value,expected in zip(re.search(r"axis cs:([\d.]+),\.36\)--\(axis cs:([\d.]+),\.36",f5).groups(),(lo,hi)):
                number(value,10**14*expected)
    run61 = read(WORK/"T2_2_PILOTO_M6_EJECUCION_20261115/corrida/RESULTADO_T2_2.json")
    number(re.search(r"coordinates \{\(([\d.]+),\.36\)\}",f5)[1],10**14*F(run61["estimacion"]["por_tasa"]["1/100000"]["x"]))
    f6 = figure("fig:strip-rate-target")
    grid = re.search(r"power comp61 comp127 comp259 comp523 global523\s*([^}]+)\}",f6)[1]
    for row,power in zip(grid.strip().splitlines(),(-6,-5,-4,-3),strict=True):
        tag,*values = row.split()
        assert int(tag)==power
        for printed,length in zip(values[:4],(6,12,24,48),strict=True):
            number(printed,cell(length,power)[2])
        number(values[4],cell(48,power)[3])
    return checked

if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Use Python con las aserciones activas.")
    sys.set_int_max_str_digits(100000)
    coordinates = plot_coordinates()
    data = {"Bernoulli":finite_codes(), "F17":affine_17(), "splitting":comparisons(),
            "opening_figure":opening_figure(), "coordenadas_cotejadas":coordinates,
            "coordenadas_cotejadas_figuras_4_5_6":coordinates}  # Zenodo v2 field name
    print(json.dumps(data, indent=2))

