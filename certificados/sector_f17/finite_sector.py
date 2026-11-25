"""Produce or independently verify the fixed 1214-background F17 certificate.

Python 3.10+, standard library only. No global circuit is simulated.
Usage: python finite_sector.py produce certificate.jsonl
       python finite_sector.py verify certificate.jsonl
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import json, sys, time

def domain():
    result=[]
    for k in range(6):
        for xs in combinations(range(2,17),k):
            b=set(xs)
            k1={13,16}<=b
            k2={4,13}<=b and 16 not in b
            u={13,15}<=b and 4 not in b and 16 not in b
            if k==0 or k1 or k2 or (k<=3 and not u): result.append(xs)
    assert len(result)==1214
    return result

def masks():
    lines=[sum(1<<(17*x+(a*x+c)%17) for x in range(17))
           for a in range(17) for c in range(17)]
    lines += [sum(1<<(17*x+y) for y in range(17)) for x in range(17)]
    return lines

def initial(lines,xs):
    s=0
    for i in range(291,306): s ^= lines[i]
    return s ^ (1<<2) ^ (1<<5) ^ (1<<22) ^ sum(1<<(17*x) for x in xs)

def classification(residual, physical):
    if physical: return 'detectable'
    return 'silent' if residual else 'success'

def produce(path):
    start=time.monotonic(); lines=masks(); count=Counter()
    with path.open('x',encoding='utf-8') as out:
        for xs in domain():
            assert time.monotonic()-start<60, '60-second reproduction bound'
            s=initial(lines,xs); correction=0; trace=[]
            while True:
                scores=[(s&line).bit_count() for line in lines]
                h=max(scores); winner=scores.index(h) if h>=9 else None
                trace.append({'counts':scores,'winner':winner})
                if winner is None: break
                s ^= lines[winner]; correction ^= 1<<winner
                assert len(trace)<=259
            residual=correction ^ sum(1<<i for i in range(291,306))
            physical=0
            for i,line in enumerate(lines):
                if residual>>i&1: physical ^= line
            label=classification(residual,physical); count[label]+=1
            out.write(json.dumps({'B':xs,'trace':trace,'residual_lines':[i for i in range(306) if residual>>i&1],
                                 'label':label},separators=(',',':'))+'\n')
    assert count=={'silent':810,'detectable':404},count
    print(json.dumps({'operation':'produce','cases':sum(count.values()),'classes':count,'seconds':time.monotonic()-start}))

def verify(path):
    # Independent geometry: build lines from each point's incident indices.
    start=time.monotonic()
    incidence=[[17*a+(y-a*x)%17 for a in range(17)]+[289+x]
               for x in range(17) for y in range(17)]
    pointsets=[set() for _ in range(306)]
    for p,ids in enumerate(incidence):
        for i in ids: pointsets[i].add(p)
    all_backgrounds=[xs for k in range(6) for xs in combinations(range(2,17),k)]
    # A second domain predicate: omitted U and otherwise unselected weights 4/5.
    expected=[]; omitted=Counter()
    for xs in all_backgrounds:
        b=set(xs); known=(13 in b and (16 in b or 4 in b))
        omit=(13 in b and 15 in b and 4 not in b and 16 not in b) or (len(b)>=4 and not known)
        if omit: omitted[(len(b-{16}),int(16 in b))]+=1
        else: expected.append(xs)
    assert len(expected)==1214 and sum(omitted.values())==3730
    assert omitted=={(2,0):1,(3,0):11,(4,0):935,(5,0):1782,(3,1):286,(4,1):715}
    classes=Counter(); polynomial=Counter(); cases=0; states=0
    with path.open(encoding='utf-8') as f:
        for line in f:
            assert time.monotonic()-start<60, '60-second verification bound'
            r=json.loads(line); xs=tuple(r['B']); assert xs==expected[cases]; cases+=1
            s={17*x+y for x in range(2,17) for y in range(17)} ^ {2,5,22} ^ {17*x for x in xs}
            correction=set(); initial_weight=len(s)
            for k,step in enumerate(r['trace']):
                counts=[0]*306
                for p in s:
                    for i in incidence[p]: counts[i]+=1
                assert counts==step['counts']; states+=1
                h=max(counts); winner=counts.index(h) if h>=9 else None
                assert winner==step['winner']
                if winner is None:
                    assert k==len(r['trace'])-1
                else:
                    old=len(s);s.symmetric_difference_update(pointsets[winner]);correction.symmetric_difference_update({winner})
                    assert old-len(s)==2*h-17>0
                    assert k<len(r['trace'])-1
            assert r['trace'][-1]['winner'] is None and len(r['trace'])<=initial_weight+1
            residual=correction ^ set(range(291,306)); physical=set()
            for i in residual: physical.symmetric_difference_update(pointsets[i])
            assert physical==s ^ {2,5,22} ^ {17*x for x in xs}
            label='detectable' if physical else ('silent' if residual else 'success')
            assert label==r['label'] and sorted(residual)==r['residual_lines']
            classes[label]+=1
            if label=='silent':
                assert residual==set(range(17))|set(range(289,306))
                polynomial[(len(set(xs)-{16}),int(16 in xs))]+=1
    assert cases==len(expected) and classes=={'silent':810,'detectable':404}
    assert polynomial=={(0,0):1,(1,0):2,(2,0):21,(3,0):100,(4,0):66,(5,0):220,
                       (1,1):3,(2,1):33,(3,1):78,(4,1):286}
    print(json.dumps({'operation':'verify','cases':cases,'states':states,'classes':classes,
                     'omitted':sum(omitted.values()),'seconds':time.monotonic()-start,'status':'PASS'}))

if __name__=='__main__':
    mode,path=sys.argv[1:]; {'produce':produce,'verify':verify}[mode](Path(path))
