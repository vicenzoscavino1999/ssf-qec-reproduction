"""Check the strip cutoff and rare-event comparisons against the recorded data."""
import json
import sys
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from math import ceil, comb
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent

def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))

def binomial_cdf(x, n, p):
    with localcontext() as ctx:
        ctx.prec = 65
        term = (1-p)**n
        total = term
        for k in range(1, x+1):
            term = term*(n-k+1)*p/(k*(1-p))
            total += term
        return +total

def main():
    sys.set_int_max_str_digits(100_000)  # Extremos racionales del fichero sellado de tiras.
    strips = ROOT/'INTEGRACION_TIRAS_Y_CERTIFICADO_CIRCUITAL_20261116'/'STRIP_COMPONENT_CERTIFICATES_20261116'/'results'
    circuit = read(strips/'circuit.json')['instancias']
    cell = circuit['TIRA_48x6']['celdas']['1/100000|1/100000']['por_k']['4']
    lo, hi = map(F,cell['intervalo_P_sil_exacto'])
    comp = (hi-lo)/(hi+lo)
    tail, mass = F(cell['eps_ing']), F(cell['L_ing'])
    glob = tail/(2*mass+tail)
    assert abs(glob-F(cell['e_ing'])) < F('1e-22')
    assert [ceil(100000*e) for e in (glob,comp)] == [99513,1221]

    strip_runs = ROOT/'T2_2_PILOTO_M6_EJECUCION_20261115'/'corrida'
    old = read(strip_runs/'RESULTADO_T2_2.json')
    assert round(read(strip_runs/'SUPERVISION_T2_2.json')['segundos']/3600,2) == 2.49
    for rate, low_pct, high_pct in [('1/100000','34.46','36.07'),('1/1000000','36.10','36.13')]:
        row = old['estimacion']['por_tasa'][rate]
        x = F(row['x'])
        lo,hi = map(F,row['J'])
        assert [lo,hi] == list(map(F,old['referencias']['J'][rate]))
        ref = circuit['TIRA_6x6']['celdas'][rate+'|'+rate]['por_k']['4']
        assert [lo,hi] == list(map(F,ref['intervalo_P_sil_exacto']))
        assert F(low_pct) <= 100*(x/hi-1) <= 100*(x/lo-1) <= F(high_pct)
        if rate == '1/100000':
            assert F('8.1358e-14') <= lo <= hi <= F('8.2332e-14')
            assert format(float(x),'.5e') == '1.10704e-13'

    split = ROOT/'SPLITTING_HGP_225_VV_20261124'
    row = read(split/'resultado_splitting.json')
    budget = read(split/'tiempos_splitting.json')
    timing = read(split/'certificado_225_tiempos.json')
    assert round(sum(budget.values()),2) == 161.91
    assert round(median(run['pared_s'] for run in timing['corridas']),2) == 8.31
    ref = row['referencia_conjunta']
    lo,hi = F(ref['L']),F(ref['U'])
    x = F(str(row['estimacion']))
    assert F('3.33938e-16') <= lo <= hi <= F('3.34508e-16')
    assert format(float(x),'.5e') == '4.91234e-16'
    assert F('46.85') <= 100*(x/hi-1) <= 100*(x/lo-1) <= F('47.11')
    assert [F(ref[k]) for k in ('L','U')] == [F(timing['corridas'][0]['intervalo_conjunto'][k]) for k in ('L','U')]

    out = read(HERE/'resultado.json')
    sample, census = out['sample'], out['census']
    S = 225*224//2*108
    assert sample['calls'] == 10**6 and sample['hits'] == 119
    assert census['calls'] == S and census['hits'] == 324
    assert round(sample['process_s'],3) == .345 and round(census['process_s'],3) == .795
    assert round(read(HERE/'preparacion.json')['build_s'],3) == 1.535
    f = F(324,S)
    assert f == F(1,8400) and (1-f)/f == 8399
    assert F(out['var_f_exact']) == f*(1-f)/sample['calls']
    assert 1-F(sample['hits'],sample['calls'])/f == F(1,2500)
    low,high = map(D,out['cp95'])
    # Verificación separada de las dos colas por suma de masas binomiales.
    assert abs(1-binomial_cdf(118,1000000,low)-D('0.025')) < D('1e-11')
    assert abs(binomial_cdf(119,1000000,high)-D('0.025')) < D('1e-11')
    assert D('.98582e-4') <= low < high <= D('1.42400e-4')
    p,q = F('1e-6'),F('1e-4')
    contribution = comb(225,2)*108*p*p*q*(1-p)**223*(1-q)**107*f
    pilot = read(ROOT/'piloto_HGP_225_esquinas_20261124.json')
    reference = next(r for r in pilot['Bernoulli'] if F(r['p']) == p and F(r['q']) == q)
    assert 324 == pilot['patrones_de_pareja']
    assert F('.9823') <= contribution/F(reference['U']) <= contribution/F(reference['L']) <= F('.9892')
    assert out['normal_planning_calls'] == {'1/10':3226442,'1/100':322644127}
    print('Strip cutoff, splitting estimates and subset-sampling results verified, including both binomial confidence limits.')

if __name__ == '__main__':
    main()
