"""Check deterministic totals, full class histograms, independent censuses and the prior pilot."""
import json
import sys
from decimal import Decimal as D
from collections import Counter
from itertools import combinations
from pathlib import Path
from modelo import bits,seed
from mayorante import cell_choices


def sectors(file):
    out={}
    for row in file.read_text().splitlines():
        col=row.split()
        if row.startswith(('2 ','3 ')):
            out[int(col[0]),int(col[1])]=(int(col[4]),col[9] if len(col)>9 else '')
        if row.startswith('triples '):
            out[3,0]=(int(col[2]),col[9] if len(col)>9 else '')
    return out


def check(candidate,archive):
    # Compare the local parity recurrence to labelled subset enumeration.
    for d in (1,2,3):
        for a in range(5):
            for b in range(min(3,d)+1):
                hist=Counter()
                for data in combinations(range(5*d),a):
                    parity=0
                    for site in data:
                        parity^=1<<(site//5)
                    for readings in combinations(range(d),b):
                        hist[parity^sum(1<<j for j in readings)]+=1
                for missing in range(d+1):
                    assert cell_choices(d,missing,a,b)==hist[(1<<missing)-1]
    for n in (225,900):
        a=json.loads((candidate/f'resultado_{n}.json').read_text())
        b=json.loads((archive/f'resultado_{n}.json').read_text())
        assert a['interval']==b['interval']
        t=a['geometry']['sheets']
        cols=seed(t)
        A,B=len(cols),9*t
        row=[[j for j,h in enumerate(cols) if h>>r&1] for r in range(B)]
        for r,js in enumerate(row):
            for j,h in enumerate(cols):
                inside={A*i+j for i in js}|{A*A+B*r+c for c in bits(h)}
                outer=[]
                for i in js:
                    for c in bits(h):
                        touched={A*i+k for k in row[c]}|{A*A+B*k+c for k in bits(cols[i])}
                        assert len(touched)==7 and len(touched&inside)==2
                        outer.extend(touched-inside)
                assert len(outer)==len(set(outer))==60
        for key in ('source_configurations_or_orbit_queries','winner_queries','local_table_probes'):
            assert a[key]==b[key]
        result=sectors(candidate/f'{n}_pairs.txt')|sectors(candidate/f'{n}_triples.txt')
        direct=sectors(archive/f'directo_{n}.txt')|sectors(archive/f'pares_sin_filtro_{n}.txt')
        assert result==direct and set(result)=={(2,0),(2,1),(2,2),(3,0)}
        assert int(a['interval']['positive_census']['3,0'])==result[3,0][0]
        ob=(candidate/f'{n}_one_readout.txt').read_text().split('orbits ')[-1].split()
        ref=(archive/f'{n}_one_readout.txt').read_text().split('orbits ')[-1].split()
        assert ob[:5]==ref[:5]
        assert D(a['interval']['max_relative_percent'])<10
        if n==225:
            prior=json.loads((archive/'piloto_previo_225.json').read_text())
            assert D(a['interval']['lower'])<=D(prior['L'])<D(prior['U'])<=D(a['interval']['upper'])
            assert D(prior['L'])<=D(a['interval']['estimate'])<=D(prior['U'])
    print('OK: two complete certificates, logical histograms, direct censuses and the 225-qubit pilot.')


if __name__=='__main__':
    if not __debug__:
        raise RuntimeError('Keep Python assertions enabled.')
    check(Path(sys.argv[1]),Path(sys.argv[2]))

