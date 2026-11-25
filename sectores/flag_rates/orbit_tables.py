"""Orbit tables for fixed numbers of VV, CC, and readout positions, using the released seed automorphisms."""
import argparse
from itertools import permutations, product
from math import prod
from pathlib import Path
import struct
import sys

def axis_rows(group,kinds,perms):
    sizes=[len(group[0][kind]) for kind in kinds]
    indices={}
    records=[]
    for z in product(*(range(size) for size in sizes)):
        if z in indices:
            continue
        orbit={tuple(g[k][v] for k,v in zip(kinds,z)) for g in group}
        i=len(records)
        indices.update(dict.fromkeys(orbit,i))
        records.append((z,len(orbit)))
    assert sum(size for _,size in records)==len(indices)==prod(sizes)
    return [(z,size,[indices[tuple(z[i] for i in pi)] for pi in perms]) for z,size in records]

def prepare(repo,path,t,a,b):
    sys.path.insert(0,str(repo/'sectores'))
    from orbitas import group
    sym=group(t)
    with Path(path).open('wb') as out:
        def pack(fmt,*values):out.write(struct.pack('<'+fmt,*values))
        out.write(b'ABORBIT1')
        pack('3I',t,a,b)
        for v in range(a+1):
            perms=[x+y+z for x in permutations(range(v)) for y in permutations(range(v,a))
                   for z in permutations(range(a,a+b))]
            pack('I',len(perms))
            first=(0,)*v+(1,)*(a-v)
            previous=None
            for side in (0,1):
                kinds=first+(side,)*b
                rows=previous if previous is not None and b==0 else axis_rows(sym,kinds,perms)
                previous=rows
                pack('I',len(rows))
                for coords,size,perm_rows in rows:
                    pack('B'*len(coords),*coords)
                    pack('I'*(1+len(perms)),size,*perm_rows)
if __name__=='__main__':
    repo,output,t,a,b=sys.argv[1:]
    prepare(Path(repo),Path(output),int(t),int(a),int(b))
