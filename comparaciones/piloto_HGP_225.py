"""Piloto del estimador directo en la incidencia afín de orden tres. Python 3.10+."""
from argparse import ArgumentParser
from array import array
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from math import comb
from platform import system
import sys
from time import perf_counter


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def mask(items):
    return sum(1 << i for i in items)


def basis(vectors):
    rows = {}
    for v in vectors:
        while v:
            j = v.bit_length()-1
            if j not in rows:
                rows[j] = v
                break
            v ^= rows[j]
    return rows


def xor_at(table, support):
    value = 0
    for i in bits(support):
        value ^= table[i]
    return value


def build():
    # Mismo orden de líneas, puntos y desempate que check_shared_background de v2.
    lines = [tuple(3*x+(s*x+b)%3 for x in range(3)) for s in range(3) for b in range(3)]
    lines += [tuple(3*x+y for y in range(3)) for x in range(3)]
    at = [[j for j, line in enumerate(lines) if c in line] for c in range(9)]
    H = [mask(js) for js in at]
    cols = [mask(line) for line in lines]
    assert all(len(x)==4 for x in at) and all(len(x)==3 for x in lines)
    assert all(len(set(x) & set(y))==1 for x,y in combinations(at,2))
    assert len(basis(H))==9
    classical = [c for c in range(1 << 12) if xor_at(cols,c)==0]
    assert Counter(c.bit_count() for c in classical)=={0:1,6:6,12:1}
    check = [mask(9*i+c for c in lines[j]) for i in range(12) for j in range(12)]
    check += [mask(9*i+c for i in at[r]) for r in range(9) for c in range(9)]
    generators = [[12*i+j for i in at[r]] + [144+9*r+c for c in lines[j]]
                  for r in range(9) for j in range(12)]
    hz = [mask(qubits) for qubits in generators]
    hx = [mask(v for v,psi in enumerate(check) if (psi >> l) & 1) for l in range(108)]
    assert len(basis(hx))==len(basis(hz))==108 and all(xor_at(check, row)==0 for row in hz)
    cb = list(basis(classical).values())
    info = next(J for J in combinations(range(12),3)
                if len(basis(mask(i for i,j in enumerate(J) if (c >> j) & 1) for c in cb))==3)
    zlogical = [mask(12*i+j for i in bits(c)) for c in cb for j in info]
    assert all((z & row).bit_count()%2==0 for z in zlogical for row in hz)
    assert len(basis(hx+zlogical))==117
    tag = [mask(a for a,z in enumerate(zlogical) if (z >> i) & 1) for i in range(225)]
    witness = next(c for c in classical if c.bit_count()==6)
    assert xor_at(check,witness)==0 and xor_at(tag,witness)!=0

    rect_checks, best_rep = [], {}
    for g, qubits in enumerate(generators):
        r,j = divmod(g,12)
        rect_checks.append([9*i+c for i in at[r] for c in lines[j]])
        for sub in range(1,128):
            action = mask(qubits[l] for l in bits(sub))
            phi = xor_at(check,action)
            if phi and 2*phi.bit_count() >= 3*action.bit_count():
                old = best_rep.get(phi)
                if old is None or (-action.bit_count(),action)<(-old.bit_count(),old):
                    best_rep[phi] = action
    actions = sorted(best_rep.items(),key=lambda x:(-x[1].bit_count(),x[1]))
    assert min(phi.bit_count() for phi,_ in actions)==3
    assert {(phi,f) for phi,f in actions if phi.bit_count()==3}=={(check[i],1 << i) for i in range(144)}
    idx = {phi:k for k,(phi,_) in enumerate(actions)}
    scores = sorted({Q(h,n) for phi,_ in actions for n in [phi.bit_count()]
                     for h in range(n+1) if 100*(2*h-n)>=n})
    score_rank = {v:i+1 for i,v in enumerate(scores)}
    score_bits = (len(actions)+1).bit_length()
    low_mask = (1 << score_bits)-1
    # A igualdad de score el índice físico menor tiene mayor prioridad.
    action_phi = [phi for phi,_ in actions]
    action_tag = [xor_at(tag,f) for _,f in actions]
    local_shapes = {}
    for sub in range(1,128):
        local = 0
        for l in bits(sub):
            local ^= (7 << (3*l)) if l<4 else sum(1 << (3*i+l-4) for i in range(4))
        if local and 2*local.bit_count()>=3*sub.bit_count():
            local_shapes[local] = local.bit_count()
    for phi in local_shapes:
        rows = [(phi >> (3*i)) & 7 for i in range(4)]
        assert 0 not in rows or all(row in (0,7) for row in rows)
    local_winners = []
    for s in range(4096):
        best, choices = 0,[]
        for phi,n in local_shapes.items():
            h = (phi & s).bit_count()
            if 100*(2*h-n)<n:
                continue
            rank = score_rank[Q(h,n)]
            if rank>best:
                best,choices=rank,[phi]
            elif rank==best:
                choices.append(phi)
        local_winners.append((best,choices))
    lookups, readers = [],[[] for _ in range(108)]
    for g, sites in enumerate(rect_checks):
        for l,c in enumerate(sites):
            readers[c].append((g,1 << l))
        priority={phi:len(actions)-idx[mask(sites[l] for l in bits(phi))] for phi in local_shapes}
        lookups.append(array('I',((rank << score_bits)+max((priority[phi] for phi in choices),default=0)
                                  for rank,choices in local_winners)))
    assert all(len(row)==12 for row in readers)

    @lru_cache(maxsize=16384)
    def run(s):
        states = {}
        for c in bits(s):
            for g,bit in readers[c]:
                states[g] = states.get(g,0)|bit
        chosen=max((lookups[g][v] for g,v in states.items()),default=0)
        if not chosen:
            return s,0
        k=len(actions)-(chosen & low_mask)
        terminal,logical=run(s ^ action_phi[k])
        return terminal, logical ^ action_tag[k]

    def direct(s):
        logical=0
        while True:
            best,h0,n0=-1,0,1
            for i,(phi,_) in enumerate(actions):
                h,n=(s & phi).bit_count(),phi.bit_count()
                if 100*(2*h-n)>=n and h*n0>h0*n:
                    best,h0,n0=i,h,n
            if best<0:
                return s,logical
            s ^= action_phi[best]
            logical ^= action_tag[best]

    for s in [0,1,mask((0,1,4,8)),mask(range(9))]+[check[i]^check[224-i]^(1 << (i%108)) for i in range(0,110,7)]:
        assert run(s)==direct(s)
    return H,check,tag,witness,len(actions),run


def row_sectors(check, tag, run, *, patterns=None):
    """Suma positiva, después de la reducción de dos errores a su fila VV."""
    by_row = []
    for i in range(12):
        hits = [0, 0]
        for u,v in combinations(range(12*i, 12*i+12), 2):
            for d in [0]+[1 << (9*i+c) for c in range(9)]:
                end, correction = run(check[u] ^ check[v] ^ d)
                if end == d and correction != tag[u] ^ tag[v]:
                    hits[int(d != 0)] += 1
                    if patterns is not None:
                        patterns.append((i,u,v,d))
        by_row.append(tuple(hits))
    assert by_row == [(0,27)]*12
    return sum(b for _,b in by_row)


def bounds(count):
    p, q = Q(1,10**12), Q(1,1000)
    pair = p*p*(1-p)**223
    low = count*pair*q*(1-q)**107
    rest = comb(225,3)*p**3
    high = count*pair*q*(1-q)**8 + rest
    estimate = 2*low*high/(low+high)
    error = (high-low)/(high+low)
    p0 = (1-p)**225*((1-q)**108 + 108*q*(1-q)**107)
    assert 0 < low < estimate < high < p0 <= 1
    assert (estimate-low)/low == (high-estimate)/high == error < Q(1,20)
    # Redondeo hacia fuera de los intervalos publicados en la ficha.
    assert Q("2.9110e-25") < low and high < Q("3.2142e-25")
    assert Q("0.9946146") < p0 and high/p0 < Q("3.23160e-25")
    displayed = Q("3.05513e-25")
    assert max(displayed/low-1, 1-displayed/high) < Q(".0495")
    return low, high, estimate, error, p0, high/p0, rest


def sector_counts(a, b, check, tag, run):
    accepted = silent = count = 0
    assert (a,b) in ((2,0),(2,1),(3,0))
    for measurement in combinations(range(108),b):
        d = mask(measurement)
        if a == 2:
            entries = ((check[i]^check[j],tag[i]^tag[j]) for i,j in combinations(range(225),2))
        else:
            entries = ((check[i]^check[j]^check[k],tag[i]^tag[j]^tag[k])
                       for i,j,k in combinations(range(225),3))
        for s,e in entries:
            end,correction = run(s ^ d)
            count += 1
            if end == d:
                accepted += 1
                silent += correction != e
    assert count == comb(225,a)*comb(108,b)
    return count,accepted,silent


def upper_corner(H, check, tag, run, patterns, triples):
    """Intervalo en el extremo uniforme con dos medidas o tres datos desglosados."""
    p, q = Q(375,10**9), Q(1,1000)
    both = 0
    for i,u,v,d in patterns:
        for bit in range(108):
            if bit//9 == i:
                continue
            measurement = d | (1 << bit)
            terminal,correction = run(check[u] ^ check[v] ^ measurement)
            if terminal == measurement and correction != tag[u] ^ tag[v]:
                both += 1
    assert len(patterns)==324 and both==324*99
    assert triples==(comb(225,3),1790256,36)
    # Tres filas con datos y la cuarta con la medida no forman las cuatro
    # rectas de ningún punto. Sólo pueden ganar filas completas hasta limpiarse.
    pencils = set(H)
    safe_row_choices = sum(
        1 for three in combinations(range(12),3) for fourth in range(12)
        if fourth not in three and (mask(three) | (1 << fourth)) not in pencils
    )
    assert safe_row_choices == comb(12,3)*9 - 9*4
    safe = safe_row_choices * 12**3 * 9
    assert safe == 30233088
    pair = p**2*(1-p)**223
    triple = p**3*(1-p)**222
    m0,m1,m2 = (q**k*(1-q)**(108-k) for k in range(3))
    low = pair*(len(patterns)*m1 + both*m2) + triple*triples[2]*m0
    top_two = pair*len(patterns)*q*(1-q)**8
    top_three = triple*(triples[2]*m0 + triples[0]*(1-m0) - safe*m1)
    rest = comb(225,4)*p**4
    high = top_two + top_three + rest
    estimate = 2*low*high/(low+high)
    error = (high-low)/(high+low)
    p0 = (1-p)**225*(m0+108*m1)
    assert Q('4.49916e-14') < low < estimate < high < Q('5.38835e-14')
    assert p0 > Q('.99453') and high/p0 < Q('5.41798e-14')
    displayed = Q('4.90377e-14')
    assert max(displayed/low-1,1-displayed/high) < Q('.09')
    assert (estimate-low)/low == (high-estimate)/high == error < Q('.09')
    return low,high,estimate,error,p0,high/p0,rest


def peak_memory():
    if system() == 'Windows':
        import ctypes
        from ctypes import wintypes
        class MemoryCounters(ctypes.Structure):
            _fields_ = [('cb',wintypes.DWORD), ('PageFaultCount',wintypes.DWORD)] + [
                (n,ctypes.c_size_t) for n in ('PeakWorkingSetSize','WorkingSetSize',
                'QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage',
                'QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage')]
        counters = MemoryCounters()
        counters.cb = ctypes.sizeof(counters)
        get_process = ctypes.windll.kernel32.GetCurrentProcess
        get_process.restype = wintypes.HANDLE
        fn = ctypes.windll.psapi.GetProcessMemoryInfo
        fn.argtypes = [wintypes.HANDLE, ctypes.POINTER(MemoryCounters), wintypes.DWORD]
        if not fn(get_process(),ctypes.byref(counters),counters.cb):
            raise ctypes.WinError()
        return counters.PeakWorkingSetSize
    import resource
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == 'darwin' else 1024*value


def exhaustive_checks(check, tag, run, triples=None):
    """Contraste opcional por fuerza bruta. Cada sector usa el host entero."""
    # En una sola fila iluminada, su recta singleton domina toda huella
    # que también salga de ella. La prioridad entre singletons es el índice.
    lines = check[:12]
    assert all(0 < line < 512 for line in lines)
    row_reduction = []
    for initial in range(512):
        state, changed = initial, 0
        while True:
            lit, neg_index = max(((state & line).bit_count(),-j) for j,line in enumerate(lines))
            if lit < 2:
                break
            j = -neg_index
            state ^= lines[j]
            changed ^= 1 << j
        row_reduction.append((state, changed))
    for i in range(12):
        for s,(end,changed) in enumerate(row_reduction):
            assert run(s << (9*i)) == (end << (9*i), xor_at(tag,changed << (12*i)))
    print("Contraste de todos los 512 síndromes de una fila en sus 12 posiciones.",flush=True)
    for a,b,expected in [(2,0,(24876,0)), (2,1,(2679129,324)), (3,0,(1790256,36))]:
        count,accepted,silent = triples if a==3 and triples is not None else sector_counts(a,b,check,tag,run)
        assert count == comb(225,a)*comb(108,b) and (accepted,silent)==expected
        print(f"Contraste |E|={a}, |D|={b}: {count} entradas, {accepted} aceptadas, {silent} silenciosas.",flush=True)


if __name__=='__main__':
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O.")
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--contraste',action='store_true',help='Evaluar también tres sectores completos por separado.')
    parser.add_argument('--extremo',action='store_true',help='Evaluar el punto p=3.75e-7, q=1e-3.')
    args = parser.parse_args()
    start=perf_counter()
    H,check,tag,witness,action_count,run=build()
    print("Filas de H:", ','.join(f'{x:03x}' for x in H), flush=True)
    print("Testigo VV:",list(bits(witness)),"huellas físicas:",action_count,flush=True)
    prepared=perf_counter()
    patterns = [] if args.extremo else None
    count = row_sectors(check,tag,run,patterns=patterns)
    triples = None
    if args.extremo:
        triples = sector_counts(3,0,check,tag,run)
        low,high,estimate,error,p0,qhigh,rest=upper_corner(H,check,tag,run,patterns,triples)
        print(f"p=3.75e-7, q=1e-3. 32076 fallos de pareja con un segundo bit exterior. "
              f"Cota de cuatro o más datos {float(rest):.12g}.")
    else:
        low,high,estimate,error,p0,qhigh,rest=bounds(count)
        print(f"p=1e-12, q=1e-3. 7920 entradas de sector, 324 con fallo. "
              f"Cota de tres o más datos {float(rest):.12g}.")
    done=perf_counter()
    print(f"P_sil entre {float(low):.14g} y {float(high):.14g} (decimales aproximados).")
    print(f"Estimación {float(estimate):.14g}, error relativo garantizado {100*float(error):.8f} %.")
    print(f"P0 >= {float(p0):.12g}, extremo superior de Q_sil {float(qhigh):.12g} (aprox.).")
    print(f"Preparación {prepared-start:.6f} s. Evaluación {done-prepared:.6f} s. "
          f"Total {done-start:.6f} s. Memoria pico {peak_memory()/2**20:.3f} MiB.",flush=True)
    if args.contraste:
        exhaustive_checks(check,tag,run,triples)
