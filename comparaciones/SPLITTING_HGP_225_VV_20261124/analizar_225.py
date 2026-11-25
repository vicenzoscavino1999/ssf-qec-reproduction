"""Recalcula el transporte, la precisión y el contraste a partir de la corrida guardada."""
import csv
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path
import statistics

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MATHS = run_path(str(ROOT/"T2_2_CODIGO_PILOTO_20261115/t2s_mat.py"))
CP = run_path(str(ROOT/"T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/t2_cp.py"))
GRID = 100_000_000
SIDES = ("left","right")
SEEDS, REPS = 12, 3

def half_values(row, side):
    prefix = SIDES[side]
    n = row["n"]//2
    values = []
    for h in (0,1):
        sx, ss = row[f"{prefix}{h}"], row[f"{prefix}{h}sq"]
        mean = sx/n
        centered = ss - sx*mean
        assert centered >= -2e-8*max(1,ss), (row["j"],centered)
        values.append((mean, max(0.0,centered/(n-1))))
    return values

def split_rhat(cohort,side):
    # Compare equal-length halves. Extended trajectories are grouped by length.
    n = cohort[0]["n"]//2
    vals = [v for row in cohort for v in half_values(row,side)]
    within = statistics.mean(v for _,v in vals)
    between_n = statistics.variance(m for m,_ in vals)
    if within < 1e-20:
        return "infinito" if between_n > 1e-18 else "indeterminado"
    return math.sqrt((n-1)/n+between_n/within)

def relative_range(x, lo, hi):
    maximum = max(abs(x/lo-1),abs(x/hi-1))
    minimum = 0 if lo<=x<=hi else min(abs(x/lo-1),abs(x/hi-1))
    return [float(100*minimum),float(100*maximum)]

def analyze():
    levels=json.loads((HERE/"modelo.json").read_text("utf-8"))["tasas"]
    last=len(levels)-1
    anchor_rows=(HERE/"ancla.txt").read_text("ascii").splitlines()
    fields=anchor_rows[0].split()
    pilot_n,pilot_hits,n,hits=map(int,fields[:4])
    assert n==min(200_000_000,max(1_000_000,(1600*(pilot_n+1)+pilot_hits)//(pilot_hits+1)))
    assert len(anchor_rows)==SEEDS+1 and all(len(row.split())==333 for row in anchor_rows[1:])
    iv=CP["cp_interval"](hits,n,alpha=".025")
    anchor={"l":F(iv["p_inf"]),"u":F(iv["p_sup"])}
    p0=F(hits,n)
    rows={}
    with (HERE/"cadenas.tsv").open() as src:
        for record in csv.DictReader(src,delimiter="\t"):
            row={k: (int(v) if k in ("seed","rep","j","n","moves","oracles") else float(v))
                 for k,v in record.items()}
            key=(row["seed"],row["rep"],row["j"])
            assert key not in rows and row["n"] in (100000,300000,700000,1500000,3100000)
            assert 0<=row["moves"]<=row["oracles"]<=row["n"]
            rows[key]=row
    assert set(rows)==set(itertools.product(range(SEEDS),range(REPS),range(last+1)))
    logratios={(l+1,r+1):[] for l in range(SEEDS) for r in range(REPS)}
    previous={}
    def unpack(line):
        ls=line.split()
        return tuple(map(int,ls[:4])), int(ls[4])/GRID, int(ls[5])
    def logistic(z):
        return 1/(1+math.exp(z))
    with (HERE/"histogramas.tsv").open() as hfile:
        next(hfile)
        for key,group in itertools.groupby(map(unpack,hfile),key=lambda v:v[0]):
            l,r,j,side=key
            items=[(z,count) for _,z,count in group]
            row=rows[l,r,j]
            nr=row["n"]
            assert all(v>0 for _,v in items) and sum(v for _,v in items)==nr
            sign=1 if side else -1
            mean=math.fsum(count*logistic(sign*z) for z,count in items)/nr
            original=statistics.mean(x for x,_ in half_values(row,side))
            assert abs(mean-original)<3e-9
            row[SIDES[side]+"_bins"]=len(items)
            if side:
                previous[l,r]=(items,nr,j)
            else:
                upper,nh,jprev=previous.pop((l,r))
                assert jprev==j-1
                estimate=MATHS["bennett"](upper,nh,items,nr)
                assert estimate["subdesbordamientos"]==0
                logratios[l+1,r+1].append(estimate["ln_q"])
    assert not previous and all(len(xs)==last for xs in logratios.values())
    t=MATHS["student_quantile"]()
    est=MATHS["estimate_rate"](logratios,last,anchor,p0,t)
    phat=est["x"]; blo,bhi=est["I"]
    nominal=MATHS["precision_measures"](phat,(blo,bhi))["E_nom"]
    transport=MATHS["precision_measures"](est["Qbar"],(max(0,est["Q_lo"]),est["Q_hi"]))["E_nom"]
    diagnostics=[]
    for j in range(last+1):
        group=[r for (l,k,lev),r in rows.items() if lev==j]
        rh={}
        for size in sorted({v["n"] for v in group}):
            cohort=[v for v in group if v["n"]==size]
            if len(cohort)<2:
                continue
            rh[str(size)]={"cadenas":len(cohort)}
            for side in range(2):
                if (side==0 and j==0) or (side==1 and j==last): continue
                rh[str(size)][SIDES[side]]=split_rhat(cohort,side)
        diagnostics.append({"nivel":j, "tasas":levels[j], "pasos":sum(r["n"] for r in group),
                            "transiciones":sum(r["moves"] for r in group),
                            "oraculos":sum(r["oracles"] for r in group),
                            "sin_transicion":sum(r["moves"]==0 for r in group),
                            "criterio_ff_incumplido":sum(r["rule"]>.05/math.sqrt(last+1) for r in group),
                            "rhat_por_longitud":rh})
    # Reference files enter only after the probability and its own interval have been computed.
    ref=json.loads((ROOT/"piloto_HGP_225_esquinas_20261124.json").read_bytes())
    row,=[r for r in ref["circuito"] if r["orden_ejecutado"]=="VV primero"
          and F(r["rc"])==F("2e-7") and F(r["rm"])==F("1e-5")]
    joint=row["intervalo_conjunto"]; det=row["determinista"]
    low,high=F(joint["L"]),F(joint["U"])
    lower_exact,upper_exact=F(det["L"]),F(det["U"])
    errs=relative_range(phat,low,high)
    result={"metodo":"Splitting Bravyi–Vargo/Bennett, 12 semillas y tres réplicas",
            "punto":{"rc":"2e-7","rm":"1e-5","orden":"VV primero","evento":"P_sil"},
            "ancla":{"calibracion_N":pilot_n,"calibracion_aciertos":pilot_hits,"N":n,
                     "aciertos":hits,"p":float(p0),"IC_97_5":[float(anchor["l"]),float(anchor["u"])],
                     "error_relativo_ancla_pct":relative_range(p0,anchor["l"],anchor["u"])[1]},
            "niveles":last+1,"semillas":SEEDS,"replicas_por_semilla":REPS,
            "pasos_totales":sum(v["pasos"] for v in diagnostics),
            "llamadas_evento_cadenas":sum(v["oraculos"] for v in diagnostics),
            "proporcion_transporte":float(est["Qbar"]),
            "transporte_por_semilla":[float(v) for v in est["Qbar_por_semilla"]],
            "intervalo_transporte_97_5_nominal":[float(est["Q_lo"]),float(est["Q_hi"])],
            "error_solo_transporte_pct":None if transport is None else float(100*transport),
            "estimacion":float(phat),"IC_nominal_95":[float(blo),float(bhi)],
            "error_intervalo_propio_pct":None if nominal is None else float(100*nominal),
            "referencia_conjunta":joint,"referencia_determinista":det,
            "error_vs_referencia_conjunta_pct":errs,
            "error_vs_referencia_determinista_pct":relative_range(phat,lower_exact,upper_exact),
            "intervalos_se_cortan":not (bhi<low or blo>high),"diagnosticos":diagnostics}
    (HERE/"resultado_splitting.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False),encoding="utf-8")
    print(f"Splitting P_sil={float(phat):.9g}, error frente a referencia {errs} %, "
          f"error nominal {result['error_intervalo_propio_pct']} %.",flush=True)
    return result

if __name__=="__main__":
    analyze()
