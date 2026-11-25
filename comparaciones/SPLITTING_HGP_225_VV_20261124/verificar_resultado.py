"""Revisa leyes físicas, semillas y cocientes de Bennett desde los registros."""
from fractions import Fraction as F
from itertools import groupby
import json
import math
from pathlib import Path
from runpy import run_path
from statistics import mean, stdev
import struct

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

def check_physics():
    independent = run_path(str(ROOT/"auditar_esquinas_HGP_225.py"))
    _, checks, gens = independent["incidence"]()
    tags, decode, _ = independent["decoder"](checks,gens)
    calendar, layers = independent["greedy_calendar"](checks,False)
    info = json.loads((HERE/"modelo.json").read_bytes())
    assert layers == info["contrato"]["capas"]
    reserved=[]
    def word(columns):
        return columns[0] | columns[1]<<64
    with (HERE/"modelo.bin").open("rb") as src:
        def read(fmt):
            return struct.unpack(fmt,src.read(struct.calcsize(fmt)))
        assert src.read(8)==b"SSF22501"
        actions,_=read("<II")
        src.seek(108*16+108*4096*4+actions*20,1)
        tests,=read("<I")
        src.seek(tests*36,1)
        levels,=read("<I")
        for rc,rm in info["tasas"]:
            assert read("<dd")==tuple(map(float,(F(rc),F(rm))))
        for s in range(333):
            count,=read("<I")
            keys=[]
            for _ in range(count):
                a,b,c,d,tag=read("<QQQQI")
                observed,offset=word((a,b)),word((c,d))
                data_syndrome=observed^offset
                assert data_syndrome in ((0,checks[s]) if s<225 else (0,))
                e=int(data_syndrome!=0)
                if not e: assert tag==0
                keys.append((e,offset))
            assert keys[0]==(0,0) and len(keys)==len(set(keys))
            reserved.append(keys)
        checked=0
        for rc,rm in info["tasas"]:
            laws,norm,h=independent["circuit_laws"](calendar,F(rc),F(rm))
            for s,keys in enumerate(reserved):
                if s<225:
                    expected={ (e,d):F(w,norm) for e in (0,1) for d,w in laws[s][e].items() if w }
                else:
                    expected={(0,0):1-h,(0,1<<(s-225)):h}
                assert set(expected)==set(keys)
                assert read("<"+str(len(keys))+"d")==tuple(float(expected[k]) for k in keys)
                checked+=len(keys)
        assert not src.read(1) and levels==len(info["tasas"])
    seeds=(HERE/"ancla.txt").read_text().splitlines()[1:]
    for line in seeds:
        syndrome=offset=tag=0
        for s,k in enumerate(map(int,line.split())):
            e,d=reserved[s][k]
            if e:
                syndrome^=checks[s]
                tag^=tags[s]
            offset^=d
        end,correction=decode(syndrome^offset)
        assert end==offset and correction!=tag
    assert len(seeds)==12
    return checked, len(seeds)

def check_transport():
    # Independent implementation on ordinary probabilities (the producer used log-sum-exp).
    def obs(v):
        seed,rep,j,side,z,m=map(int,v.split())
        return (seed,rep,j,side), math.exp(z/1e8), m
    paths={}
    last_high=None
    with (HERE/"histogramas.tsv").open() as src:
        next(src)
        for (seed,rep,j,side),items in groupby(map(obs,src),lambda el:el[0]):
            counts=[(x,n) for _,x,n in items]
            total=sum(n for _,n in counts)
            freqs=[(x,n/total) for x,n in counts]
            if side:
                last_high=(seed,rep,j,freqs)
                continue
            a,b,old,high=last_high
            assert (a,b,old)==(seed,rep,j-1)
            last_high=None
            c=1.0
            for _ in range(3):
                numerator=math.fsum(m/(1+c*x) for x,m in high)
                denominator=math.fsum(m*c*x/(1+c*x) for x,m in freqs)
                c*=numerator/denominator
            paths.setdefault((seed,rep),[]).append(c)
    assert last_high is None and len(paths)==36 and {len(v) for v in paths.values()}=={27}
    ratios=[mean(math.prod(paths[l,r]) for r in range(3)) for l in range(12)]
    output=json.loads((HERE/"resultado_splitting.json").read_bytes())
    for new,old in zip(ratios,output["transporte_por_semilla"],strict=True):
        assert math.isclose(new,old,rel_tol=5e-12)
    scale=F(output["ancla"]["aciertos"],output["ancla"]["N"])
    x=float(scale)*mean(ratios)
    assert math.isclose(x,output["estimacion"],rel_tol=5e-12)
    ref=json.loads((ROOT/"piloto_HGP_225_esquinas_20261124.json").read_bytes())["circuito"]
    target,=[r for r in ref if (F(r["rc"]),F(r["rm"]),r["orden_ejecutado"])==(F("2e-7"),F("1e-5"),"VV primero")]
    bounds=[float(target["intervalo_conjunto"][e]) for e in ("L","U")]
    actual=[(x/bounds[1]-1)*100,(x/bounds[0]-1)*100]
    assert min(actual)>10 and output["IC_nominal_95"][0]>max(bounds)
    assert all(abs(a-b)<1e-9 for a,b in zip(actual,output["error_vs_referencia_conjunta_pct"]))
    return len(paths)*len(paths[0,0]),x

if __name__=="__main__":
    checked,seeds=check_physics()
    ratios,estimate=check_transport()
    print(f"{checked} probabilidades de etiqueta reconstruidas en los 28 niveles. {seeds}/12 semillas silenciosas por reducción estabilizadora.")
    print(f"{ratios} razones de Bennett recalculadas con aritmética directa y la misma estimación, {estimate:.12g}.")
