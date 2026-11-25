"""Independent seed symmetries on each HGP axis; only score-and-size ties are relaxed."""
from itertools import permutations,product
from math import prod
from pathlib import Path
import struct
from modelo import seed,bits,mask


def group(t):
    base = seed(1)
    affine = []
    for a,b,c,d,dx,dy in product(range(3),repeat=6):
        if (a*d-b*c)%3 == 0:
            continue
        p = tuple(3*((a*x+b*y+dx)%3)+(c*x+d*y+dy)%3 for x,y in product(range(3),repeat=2))
        line = tuple(base.index(mask(p[i] for i in bits(h))) for h in base)
        affine.append((line,p))
    if t==1:
        result=affine
    else:
        stab=[z for z in affine if z[1][0]==0 and z[0][6]==6]
        result=[]
        for s0,s1,swap in product(stab,stab,range(2)):
            l=tuple(12*(i^swap)+j for i,u in enumerate((s0,s1)) for j in u[0])
            p=tuple(9*(i^swap)+j for i,u in enumerate((s0,s1)) for j in u[1])
            result.append((l,p))
    assert len(set(result))==len(result)==(432 if t==1 else 288)
    h=seed(t)
    for l,p in result:
        assert all(mask(p[v] for v in bits(col))==h[l[j]] for j,col in enumerate(h))
    return result


def prepare(t,path,triples_only=False):
    groupmaps=group(t)
    with Path(path).open('wb') as out:
        def pack(fmt,*v):out.write(struct.pack('<'+fmt,*v))
        out.write(b'SSFZERO1' if triples_only else b'SSFORBT1')
        for nvv in range(4):
            perms=[a+b+(3,) for a in permutations(range(nvv))
                   for b in permutations(range(nvv,3))]
            pack('I',len(perms))
            for last in (0,1):
                kinds=(0,)*nvv+(1,)*(3-nvv)+(() if triples_only else (last,))
                sizes=tuple((12*t,9*t)[j] for j in kinds)
                index,rows={},[]
                for z in product(*(range(k) for k in sizes)):
                    if z in index:
                        continue
                    orbit={tuple(g[k][u] for k,u in zip(kinds,z)) for g in groupmaps}
                    row=len(rows)
                    for zz in orbit:index[zz]=row
                    rows.append((z,len(orbit),tuple(sorted(orbit)) if triples_only else ()))
                assert len(index)==prod(sizes)==sum(v for _,v,_ in rows)
                pack('I',len(rows))
                for z,sz,members in rows:
                    pack('4BI',*(z+(0,) if triples_only else z),sz)
                    pack('I'*len(perms),*(index[tuple(z[i] for i in pi[:len(z)])] for pi in perms))
                    for coords in members:
                        pack('3B',*coords)

