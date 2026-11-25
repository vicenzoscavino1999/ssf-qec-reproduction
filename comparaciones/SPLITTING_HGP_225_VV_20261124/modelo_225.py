"""Prepara las leyes físicas VV y las tablas del SSF del piloto, sin leer su resultado."""
from array import array
from fractions import Fraction as F
from inspect import getclosurevars
import json
from pathlib import Path
from random import Random
from runpy import run_path
import struct
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OLD = ROOT / "T2_2_CODIGO_PILOTO_20261115"
MASK64 = (1 << 64) - 1

def halves(value):
    return value & MASK64, value >> 64

def build(folder):
    root = Path(folder)
    circuit = run_path(str(ROOT / "circuito_HGP_225.py"))
    maths = run_path(str(OLD / "t2s_mat.py"))
    pilot = circuit["pilot"]
    H, checks, tags, _, _, run = pilot["build"]()
    data, anc = circuit["schedule"](checks, cc_first=False)
    contract = circuit["primitive_checks"](H, checks, data, anc)
    assert contract["N_loc"] == 3240
    def readout(rc):
        return min(F(".005"), 50 * rc)
    # El número medio de fibras ocupadas no supera el de fuentes primitivas.
    rates = maths["build_ladder"](lambda rc: 3240*rc + 108*readout(rc),
                                  ends=(F(".003"), F(".0001"), F("2e-7")))
    assert readout(rates[-1]) == F(".00001")
    positions = pilot["bits"]
    local = getclosurevars(run.__wrapped__).nonlocals
    rects = [[9*i+c for i in positions(row) for c in positions(col)]
             for row in H for col in checks[:12]]
    assert len(rects)==108 and all(z==sorted(z) and len(z)==12 for z in rects)
    levels, labels = [], None
    for rc in rates:
        laws, norm = circuit["joint_fibres"](data, rc)
        if rc in (rates[0], rates[-1]):
            circuit["confirm_laws"](data, laws, norm, rc)
        h = circuit["ancilla_rate"](rc, readout(rc))
        keys, vals = [], []
        for site, law in enumerate(laws):
            keys_i = sorted((e, d) for e in (0, 1) for d, count in law[e].items() if count)
            assert keys_i[0]==(0,0) and sum(law[e][d] for e,d in keys_i)==norm
            keys.append(keys_i)
            vals.append([float(F(law[e][d],norm)) for e,d in keys_i])
        for j in range(108):
            keys.append([(0,0),(0,1<<j)])
            vals.append([float(1-h),float(h)])
        assert labels is None or labels==keys
        labels = keys
        levels.append(vals)
    with (root/'modelo.bin').open('wb') as file:
        def put(fmt,*args):
            file.write(struct.pack('<'+fmt,*args))
        file.write(b'SSF22501')
        put('II',len(local["actions"]),local["low_mask"])
        for sites in rects:
            put('QQ', *halves(sum(1 << i for i in sites)))
        for row in local["lookups"]:
            a=array('I',row)
            assert a.itemsize==4
            if sys.byteorder!='little':
                a.byteswap()
            file.write(a.tobytes())
        for phi,tag in zip(local["action_phi"],local["action_tag"],strict=True):
            put('QQI',*halves(phi),tag)
        rng=Random(596142003)
        tests=[0]+[1<<j for j in range(108)]+checks
        tests.extend(sum(1<<j for j in rng.sample(range(108),rng.randrange(2,31))) for _ in range(1250))
        put('I',len(tests))
        for syn in tests:
            end,tag=run(syn)
            put('QQQQI',*halves(syn),*halves(end),tag)
        put('I',len(rates))
        for rc in rates:
            put('dd',float(rc),float(readout(rc)))
        for site,keys in enumerate(labels):
            put('I',len(keys))
            for e,d in keys:
                obs=d ^ (checks[site] if e else 0)
                put('QQQQI',*halves(obs),*halves(d),tags[site] if e else 0)
        for lev in levels:
            for row in lev:
                assert all(x>0 for x in row)
                put(str(len(row))+'d',*row)
    meta = {"tasas": [[str(x), str(readout(x))] for x in rates],
            "fibras": len(labels), "validaciones_decodificador": len(tests),
            "instancia": "[[225,9,6]]", "orden":"VV primero", "contrato":contract,
            "punto": ["2e-7","1e-5"], "codigo_heredado": ["piloto_HGP_225.py",
            "circuito_HGP_225.py","T2_2_CODIGO_PILOTO_20261115/t2s_mat.py"]}
    (root/'modelo.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    return meta

if __name__=='__main__':
    build(HERE)
