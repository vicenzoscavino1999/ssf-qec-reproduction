"""Generate the two SSF instances and their exact global competition tables."""
from array import array
from fractions import Fraction
from itertools import combinations
from functools import reduce
from operator import xor
from pathlib import Path
from random import Random
import struct
import sys


def bits(x):
    while x:
        bit = x & -x
        yield bit.bit_length() - 1
        x ^= bit


def mask(indices):
    return sum(1 << i for i in indices)


def tagsum(table, support):
    return reduce(xor, (table[i] for i in bits(support)), 0)


def kernel(columns):
    pivots, result = {}, []
    for j, h in enumerate(columns):
        word = 1 << j
        while h:
            index = h.bit_length() - 1
            if index not in pivots:
                pivots[index] = (h, word)
                break
            a, b = pivots[index]
            h, word = h ^ a, word ^ b
        if h == 0:
            result.append(word)
    return result


def seed(t):
    assert t in (1, 2)
    base = [[3*x+(s*x+b) % 3 for x in range(3)] for s in range(3) for b in range(3)]
    base += [[3*x+y for y in range(3)] for x in range(3)]
    cols = [mask(9*((i+int(j == 6 and point == 0)) % t)+point for point in points)
            for i in range(t) for j, points in enumerate(base)]
    return cols


def build(t, path):
    cols = seed(t)
    A, B = len(cols), 9*t
    at = [[j for j, h in enumerate(cols) if h >> i & 1] for i in range(B)]
    rows = [mask(js) for js in at]
    assert all(len(js) == 4 for js in at) and all(h.bit_count() == 3 for h in cols)
    assert max((h & k).bit_count() for h, k in combinations(cols, 2)) == 1
    reached = cols[0]
    while True:
        more = reduce(int.__or__, (h for h in cols if h & reached), reached)
        if more == reached:
            break
        reached = more
    assert reached.bit_count() == B
    cb = kernel(cols)
    assert len(cb) == A - B
    retained=[col for j,col in enumerate(seed(1)) if j!=6]
    assert len(retained)-len(kernel(retained))==9
    assert tagsum(cols,(1<<6)-1)==0
    code = [tagsum(cb, i) for i in range(1 << len(cb))]
    words = [w for w in code if w.bit_count() == 6]
    assert min(w.bit_count() for w in code if w) == 6
    info = [v.bit_length()-1 for v in cb]  # distinct leading bits, triangular
    assert len(set(info)) == len(cb)
    check = [col << (B*i) for i in range(A) for col in cols]
    check += [mask(B*i+j for i in js) for js in at for j in range(B)]
    n, nr = len(check), A*B
    tags = [mask(a*len(info)+k for a, w in enumerate(cb) if w >> i & 1
                 for k, d in enumerate(info) if j == d) for i in range(A) for j in range(A)]
    tags += [0]*(B*B)
    generators = [[A*i+j for i in js]+[A*A+B*r+p for p in bits(col)]
                  for r, js in enumerate(at) for j, col in enumerate(cols)]
    assert all(tagsum(check, mask(g)) == tagsum(tags, mask(g)) == 0 for g in generators)
    assert any(tagsum(tags, w << (A*i)) for w in words for i in range(A))
    rects = [[B*i+p for i in js for p in bits(col)] for js in at for col in cols]
    shapes = {}
    for sub in range(1, 128):
        h = mask(3*i+j for i in range(4) for j in range(3)
                 if (sub >> i ^ sub >> (4+j)) & 1)
        if h and 2*h.bit_count() >= 3*sub.bit_count():
            shapes.setdefault(h, []).append(sub)
    best_rep, global_shapes = {}, []
    for sites, g in zip(rects, generators):
        glob = {h: mask(sites[i] for i in bits(h)) for h in shapes}
        global_shapes.append(glob)
        for h, subs in shapes.items():
            phi = glob[h]
            for sub in subs:
                f = mask(g[i] for i in bits(sub))
                prev = best_rep.get(phi)
                if prev is not None:
                    assert tagsum(tags, prev) == tagsum(tags, f)
                if prev is None or (-f.bit_count(), f) < (-prev.bit_count(), prev):
                    best_rep[phi] = f
    acts = sorted(best_rep, key=lambda h: (-best_rep[h].bit_count(), best_rep[h]))
    assert {h for h in acts if h.bit_count() == 3} == set(check[:A*A])
    assert all(best_rep[h].bit_count() == 1 for h in acts if h.bit_count() == 3)
    priorities = {h: len(acts)-j for j, h in enumerate(acts)}
    nbits = len(acts).bit_length()
    sizes = {h.bit_count() for h in shapes}
    scores = sorted({Fraction(k, d) for d in sizes for k in range(d+1) if 100*(2*k-d) >= d})
    ranks = {(k, d): scores.index(Fraction(k, d))+1 for d in sizes for k in range(d+1)
             if 100*(2*k-d) >= d}
    choices = []
    for syn in range(4096):
        rank, tied = 0, []
        for h in shapes:
            r = ranks.get(((h & syn).bit_count(), h.bit_count()), 0)
            if r > rank:
                rank, tied = r, [h]
            elif r and r == rank:
                tied.append(h)
        choices.append((rank << nbits, tied))
    colpairs = {(u, v) for w in words for u, v in combinations(bits(w), 2)}
    pairs = [(A*i+u, A*i+v) for i in range(A) for u, v in sorted(colpairs)]
    rng = Random(721663+t)
    # Independent global scan with exact rational comparison and physical tie order.
    alist = [(h,h.bit_count(),tagsum(tags,best_rep[h])) for h in acts]
    def direct(s):
        mark = 0
        while True:
            bn, bd, action, toggle = 0, 1, 0, 0
            for h, d, tag in alist:
                k = (h & s).bit_count()
                if 100*(2*k-d) >= d and k*bd > bn*d:
                    bn, bd, action, toggle = k, d, h, tag
            if not action:
                return s, mark
            s ^= action
            mark ^= toggle
    checksums = []
    for i in range(36):
        e = mask(rng.sample(range(n),i%6))
        d = mask(rng.sample(range(nr),i%4))
        syn = tagsum(check,e)^d
        terminal, corr = direct(syn)
        checksums.append((syn,terminal,corr))
    with Path(path).open('wb') as out:
        def pack(fmt,*v):
            out.write(struct.pack('<'+fmt,*v))
        def small(i):
            pos = list(bits(i))
            pack('B',len(pos))
            pack('H'*len(pos),*pos)
        out.write(b'SSFMARK1')
        pack('6I',n,nr,len(rects),len(acts),nbits,len(pairs))
        for g in rects:
            pack('12H',*g)
        for glob in global_shapes:
            pri = {h:priorities[v] for h,v in glob.items()}
            tbl = array('I',(r+max((pri[h] for h in hh),default=0) for r,hh in choices))
            assert tbl.itemsize == 4
            if sys.byteorder != 'little':
                tbl.byteswap()
            out.write(tbl.tobytes())
        for h, _, tag in alist:
            small(h)
            pack('Q',tag)
        for h, tag in zip(check,tags):
            small(h)
            pack('Q',tag)
        for i,j in pairs:
            pack('HH',i,j)
        pack('I',len(checksums))
        for syn,end,tag in checksums:
            small(syn)
            small(end)
            pack('Q',tag)
        unique = sorted(shapes)
        pack('I',len(unique))
        for _, tied in choices:
            pack('B',len(tied))
            pack('B'*len(tied),*(unique.index(h) for h in tied))
        for glob in global_shapes:
            pack('I'*len(unique),*(len(acts)-priorities[glob[h]] for h in unique))
        out.write(bytes(best_rep[h].bit_count() for h in acts))
    return dict(sheets=t,qubits=n,checks=nr,logical_qubits=len(cb)**2,witness_weight=6,
                min_classical_distance=6,classical_weight_six=len(words),marked_pairs=len(pairs),
                actions=len(acts),global_oracle_tests=len(checksums),
                checked_seed_rank=B,seed_columns=cols,information_columns=info)

