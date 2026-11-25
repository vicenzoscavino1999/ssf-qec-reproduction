"""Deterministic upper counts for a spurious perfect first flip or blocked originals."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import comb
from modelo import bits,mask


def choose(n,k):
    return comb(n,k) if 0 <= k <= n else 0


@lru_cache(None)
def cell_choices(d,missing,x,y):
    # x external data, y readouts inside the footprint. Even/odd at each cell.
    if min(x,y) < 0:
        return 0
    if d == 0:
        return int(x == y == missing == 0)
    if missing > x+y or (missing-x-y)%2:
        return 0
    parity = int(missing>0)
    return sum(choose(5,j)*cell_choices(d-1,missing-parity,x-j,y-k)
               for j in range(min(5,x)+1) for k in range(min(1,y)+1) if (j+k)%2 == parity)


def false_perfect(n,m,a,b):
    def phi(sub):
        return mask(3*i+j for i in range(4) for j in range(3) if (sub>>i ^ sub>>(4+j))&1)
    shapes = {}
    for sub in range(1,128):
        h=phi(sub)
        if h and 2*h.bit_count()>=3*sub.bit_count():
            shapes[h]=sub
    @lru_cache(None)
    def completions(d,missing,left):
        return sum(cell_choices(d,missing,x,y)*choose(n-7-5*d,left-x)*choose(m-d,b-y)
                   for x in range(left+1) for y in range(b+1))
    groups={}
    for h,rep in shapes.items():
        cnt=0
        for sub in range(128):
            if sub.bit_count()>a or sub&rep==rep or sub&(127^rep)==(127^rep):
                continue
            d=h.bit_count()
            cnt+=completions(d,(h&~phi(sub)).bit_count(),a-sub.bit_count())
        # Row-only, single-row footprints occur in three first-axis generators;
        # single-column footprints occur in four second-axis generators.
        # Either line is uniquely identified by 2 cells. Other shapes fit once
        # unless all lit rows are in just one row or likewise for columns.
        rowids={j//3 for j in bits(h)}
        colids={j%3 for j in bits(h)}
        div=3 if len(rowids)==1 else 4 if len(colids)==1 else 1
        groups[h]=Q(m*cnt,div)
    return sum(groups.values())


def no_perfect(n,m,a,b):
    # Each original overlaps at most 24 other sites and touches at most 4 readouts.
    if a==3:
        return min(Q(1), 3*Q(24,n-1)*(Q(48,n-2)+Q(4*b,m))
                   +Q(choose(7*b,3),comb(n,3)))
    return min(Q(1),Q(24*(a-1),n-1)+Q(4*b,m))


if __name__=='__main__':
    for n,m in ((225,108),(900,432)):
        for a,b in ((3,0),(3,1),(3,2),(3,3),(4,0),(4,1)):
            f=false_perfect(n,m,a,b)
            print(n,a,b,'false_perfect',float(f/(comb(n,a)*comb(m,b))),
                  'blocked',float(no_perfect(n,m,a,b)))

