"""Muestreo y censo en procesos nuevos; genera los dos resultados con el mismo SSF."""
import argparse
import hashlib
import json
import os
from fractions import Fraction as F
from math import ceil, comb, sqrt
from pathlib import Path
from runpy import run_path
from statistics import NormalDist
from subprocess import run
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SHA_MODEL = '8581eadc0e728ef7c6d1b97f82a7e8b1807f1859e5dd7beb049c2194897b12fe'

def process(mode):
    program = HERE / ('estrato_225.exe' if os.name == 'nt' else 'estrato_225')
    start = perf_counter()
    result = run([str(program), mode], cwd=HERE, text=True, capture_output=True,
                 check=True)
    wall = perf_counter() - start
    row = json.loads(result.stdout)
    assert row['mode'] == mode
    row['process_s'] = wall
    return row

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE / 'repeticion.json')
    args = parser.parse_args()
    with (HERE / 'modelo.bin').open('rb') as source:
        assert hashlib.file_digest(source, 'sha256').hexdigest() == SHA_MODEL
    # Neither process reads the census or previous pilot counts.
    sample = process('sample')
    census = process('census')
    n, h = sample['calls'], sample['hits']
    size = comb(225, 2)*108
    assert n == 1_000_000 and sample['seed'] == 202611242151
    assert census['calls'] == size and census['seed'] is None
    reference = json.loads((ROOT / 'piloto_HGP_225_esquinas_20261124.json').read_text(encoding='utf-8'))
    assert census['hits'] == reference['patrones_de_pareja'] == 324
    pilot = next(row for row in reference['Bernoulli']
                 if F(row['p']) == F('1e-6') and F(row['q']) == F('1e-4'))
    p, q = F(pilot['p']), F(pilot['q'])
    weight = size*p**2*(1-p)**223*q*(1-q)**107
    frequency = F(census['hits'], size)
    rel_factor = (1-frequency)/frequency
    assert frequency == F(1,8400) and rel_factor == 8399
    cp_interval = run_path(str(ROOT/'T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114'/'t2_cp.py'))['cp_interval']
    interval = cp_interval(h, n, alpha='0.05', target='1e-16')
    bounds = [F(interval[key]) for key in ('p_inf','p_sup')]
    point = F(h, n)
    z = NormalDist().inv_cdf(.975)
    weighted = weight*frequency
    output = {
        'stratum': [2,1], 'data_readout_sites': [225,108], 'rates': ['1e-6','1e-4'],
        'model_sha256': SHA_MODEL, 'sample': sample, 'census': census,
        'f_exact': str(frequency), 'weight': float(weight), 'true_contribution': float(weighted),
        'full_reference_pilot': [pilot['L'], pilot['U']],
        'share_of_full_reference': [float(weighted/F(pilot['U'])),float(weighted/F(pilot['L']))],
        'sample_f': str(point), 'cp95': [str(x) for x in (interval['p_inf'],interval['p_sup'])],
        'sample_relative_error': float(abs(point-frequency)/frequency),
        'sample_relative_interval_from_point': [float((point-bounds[0])/point),float((bounds[1]-point)/point)] if h else None,
        'weighted_sample': float(weight*point), 'weighted_cp95': [float(weight*x) for x in bounds],
        'var_f_exact': str(frequency*(1-frequency)/n),
        'relative_variance_exact': str(rel_factor/n), 'relative_sd': sqrt(float(rel_factor/n)),
        'normal_planning_calls': {str(e):ceil(z*z*float(rel_factor)/float(e)**2) for e in (F(1,10),F(1,100))},
    }
    args.output.write_text(json.dumps(output, indent=2)+'\n',encoding='utf-8')
    print(json.dumps(output, indent=2))

if __name__ == '__main__':
    main()
