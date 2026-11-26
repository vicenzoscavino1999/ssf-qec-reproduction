"""Certificados de una ronda de extracción Z en el HGP afín (3,4)."""
from argparse import ArgumentParser
from bisect import bisect_right
from collections import Counter, defaultdict
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
import json
from math import comb, prod, isqrt, factorial
from pathlib import Path
import platform
from random import Random
from runpy import run_path
import sys
from time import perf_counter


HERE = Path(__file__).resolve().parent
pilot = run_path(str(HERE / "piloto_HGP_225.py"))
bits, mask = pilot["bits"], pilot["mask"]
N, M, L = 225, 108, 10


def schedule(checks, *, cc_first=False):
    """Primer color libre en cada incidencia. El orden histórico procesa primero los VV."""
    data_layers = [dict() for _ in checks]
    check_layers = [dict() for _ in range(M)]
    wire_order = list(range(N))
    if cc_first:
        wire_order = wire_order[144:] + wire_order[:144]
    for c,used in enumerate(check_layers):
        for v in wire_order:
            if checks[v] >> c & 1:
                layer = next(t for t in range(1,L+1) if t not in used and t not in data_layers[v])
                used[layer] = v
                data_layers[v][layer] = c
    assert all(len(row)==7 for row in check_layers)
    assert [len(row) for row in data_layers]==[3]*144+[4]*81
    return data_layers,check_layers


def joint_fibres(data_layers, rc):
    """Leyes enteras de (paridad X final, D efectivo local) de los 12 sitios del dato."""
    denominator, rate = 15*rc.denominator,rc.numerator
    wait_x, pair_x = 10*rate,4*rate
    normalizer = denominator**(L+2)

    @lru_cache(maxsize=None)
    def distribution(times):
        degree = len(times)
        law = {(0,0):1}
        past = 0
        for t in range(L+2):
            if t in times:
                current = 1 << times.index(t)
                past |= current
                # Un X de control omite los checks hasta t en D=s xor H_X E.
                options = [(0,0,denominator-3*pair_x),(1,past,pair_x),
                           (0,current,pair_x),(1,past ^ current,pair_x)]
            else:
                options = [(0,0,denominator-wait_x),(1,past,wait_x)]
            next_law = defaultdict(int)
            for (e,d),n in law.items():
                for a,b,mass in options:
                    next_law[e ^ a,d ^ b] += n*mass
            law = next_law
        assert sum(law.values())==normalizer and len(law)<=2**(degree+1)
        return tuple(law.get((e,d),0) for e in range(2) for d in range(1 << degree))

    laws = []
    for row in data_layers:
        times = tuple(sorted(row))
        positions = [row[t] for t in times]
        packed = distribution(times)
        stride = 1 << len(times)
        laws.append(tuple({mask(positions[i] for i in bits(d)):packed[e*stride+d]
                           for d in range(stride)} for e in range(2)))
    return laws,normalizer


def ancilla_rate(rc,rm):
    """Tres esperas, preparación y Pauli de lectura, más inversión clásica."""
    return (1-(1-4*rc/3)**5*(1-2*rm))/2


GRID = 1 << 128


def product_grid(values):
    result = GRID
    for v in values:
        result = (result*v)//GRID
    return result


def pair_spectra(laws,patterns,normalizer, *, restricted_patterns=None):
    """Walsh por fila. En la cota inferior cada cable tiene D exterior vacío."""
    sums = [[0]*10 for _ in range(2)]
    rows = [defaultdict(list) for _ in range(12)]
    for i,u,v,d in patterns:
        rows[i][u%12,v%12].append(d >> (9*i))
    for i,choices in enumerate(rows):
        if not choices:
            continue
        area = 511 << (9*i)
        for exterior_free,total in enumerate(sums):
            # Fuera de la fila, un CC toca a lo sumo una celda de la fila.
            fixed = []
            sign_factors = [[[],[]] for _ in range(9)]
            for v,law in enumerate(laws):
                if 12*i<=v<12*i+12:
                    continue
                local = defaultdict(int)
                for d,n in law[0].items():
                    if exterior_free or d & ~area==0:
                        local[(d & area) >> (9*i)] += n
                nonzero = [j for j in local if j]
                if nonzero:
                    bit, = nonzero
                    assert bit.bit_count()==1
                    factors = sign_factors[bit.bit_length()-1]
                    for parity in (0,1):
                        factors[parity].append((local[0]+(-1)**parity*local[bit])*GRID//normalizer)
                else:
                    fixed.append(local[0]*GRID//normalizer)
            constant = product_grid(fixed)
            signs = [[product_grid(row) for row in options] for options in sign_factors]
            terms = [[[(d >> (9*i),n) for d,n in law[e].items()] for e in (0,1)]
                     for law in laws[12*i:12*i+12]]
            for t in range(512):
                background = product_grid([constant]+[p[t>>j & 1] for j,p in enumerate(signs)])
                transformed = [tuple(sum(n if (t & d).bit_count()%2==0 else -n for d,n in row)
                                     *GRID//normalizer for row in pair) for pair in terms]
                for (u,v),allowed in choices.items():
                    sign = sum(1 if (t & d).bit_count()%2==0 else -1 for d in allowed)
                    if sign or (restricted_patterns is not None and not exterior_free):
                        coefficient = product_grid([background]+[
                            pair[int(j==u or j==v)] for j,pair in enumerate(transformed)])
                        total[t.bit_count()] += sign*coefficient
                        if restricted_patterns is not None and not exterior_free:
                            for d in allowed:
                                key = (i, 12*i+u, 12*i+v, d << (9*i))
                                moments = restricted_patterns.setdefault(key, [0]*10)
                                moments[t.bit_count()] += (-1)**((t & d).bit_count()%2)*coefficient
    return sums


def pair_bounds(spectra,h):
    inside = [sum(Q(x,512*GRID)*(1-2*h)**j for j,x in enumerate(row)) for row in spectra]
    # Cada producto de 225 momentos pasa por <500 redondeos de error <=1/GRID.
    # La inversión suma 324 máscaras con factores de módulo <=1.
    slack = Q(324*512,GRID)
    exterior = (1-h)**99 + 99*h*(1-h)**98
    lower,upper = (inside[0]-slack)*exterior,inside[1]+slack
    assert 0<lower<=upper
    return lower,upper


def contrast_pair(laws,norm,pattern,h):
    """Convolución positiva de una entrada de fila, independiente de la inversión Walsh."""
    row,u,v,target = pattern
    shift,area = 9*row,511 << (9*row)
    spectra = pair_spectra(laws,[pattern],norm)
    for free in (0,1):
        state = {0:GRID}
        for site,law in enumerate(laws):
            weights = defaultdict(int)
            for d,n in law[int(site==u or site==v)].items():
                if free or d & ~area==0:
                    weights[(d & area) >> shift] += n
            updated = defaultdict(int)
            for a,left in state.items():
                for b,right in weights.items():
                    updated[a^b] += left*right
            state = {a:value//norm for a,value in updated.items()}
        independent = sum(Q(n,GRID)*h**k*(1-h)**(9-k) for d,n in state.items()
                          for k in [(d^(target >> shift)).bit_count()])
        fourier = sum(Q(x,512*GRID)*(1-2*h)**j for j,x in enumerate(spectra[free]))
        # Suelo por celda en cada convolución positiva, con masa total <=1.
        assert abs(independent-fourier)<Q(512*(N+1),GRID)


def all_winners(H,checks,tags):
    """Todas las huellas físicas de máxima puntuación, con su clase."""
    rects = [[9*i+c for i in bits(row) for c in bits(line)] for row in H for line in checks[:12]]
    stripes = [7 << (3*i) for i in range(4)] + [mask(3*i+j for i in range(4)) for j in range(3)]
    shapes = {}
    for sub in range(1,128):
        phi = pilot["xor_at"](stripes,sub)
        if phi and 2*phi.bit_count()>=3*sub.bit_count():
            shapes[phi] = sub
    assert len(shapes)==63
    ranks = sorted({Q(h,n) for phi in shapes for n in [phi.bit_count()] for h in range(n+1)
                    if 100*(2*h-n)>=n})
    rank = {(h,n):ranks.index(Q(h,n))+1 if 100*(2*h-n)>=n else 0
            for n in {s.bit_count() for s in shapes} for h in range(n+1)}
    table = []
    for s in range(4096):
        scores = [(rank[(s & phi).bit_count(),phi.bit_count()],phi) for phi in shapes]
        best = max(k for k,_ in scores)
        table.append((best,[phi for k,phi in scores if k==best] if best else []))
    readers = [[] for _ in range(M)]
    actions = []
    for g,sites in enumerate(rects):
        r,j = divmod(g,12)
        data = [12*i+j for i in bits(H[r])] + [144+9*r+c for c in bits(checks[j])]
        acts = {}
        for phi,subset in shapes.items():
            qubits = mask(data[a] for a in bits(subset))
            physical = mask(sites[a] for a in bits(phi))
            assert physical==pilot["xor_at"](checks,qubits)
            acts[phi] = (physical,pilot["xor_at"](tags,qubits))
        actions.append(acts)
        for a,c in enumerate(sites):
            readers[c].append((g,1 << a))

    @lru_cache(maxsize=65536)
    def winners(s):
        states = {}
        for c in bits(s):
            for g,bit in readers[c]:
                states[g] = states.get(g,0) | bit
        best = max((table[z][0] for z in states.values()),default=0)
        return frozenset(actions[g][phi] for g,z in states.items() if table[z][0]==best
                         for phi in table[z][1])
    return winners


def terminal_classes(winners):
    """Clases de las correcciones que pueden terminar en d al seguir cualquier empate."""
    @lru_cache(maxsize=65536)
    def classes(s,d):
        if s.bit_count()<=d.bit_count():
            return frozenset([0]) if s==d and not winners(d) else frozenset()
        return frozenset(t^a for phi,t in winners(s) for a in classes(s^phi,d))
    return classes


def interior_triples(H,checks,tags):
    """Pares de órbitas con al menos una lectura interior que admite fallo."""
    geometry = run_path(str(HERE / "comparar_HGP_225.py"))
    maps = geometry["affine_maps"](checks[:12])
    classes = terminal_classes(all_winners(H,checks,tags))

    studied = hits = 0
    promising = []
    for nvv in range(4):
        reps = geometry["axis_orbits"](nvv,maps)
        for (left,nl),(right,nr) in product(reps,repeat=2):
            sites = tuple(12*r+c if i<nvv else 144+9*r+c for i,(r,c) in enumerate(zip(left,right)))
            if len(set(sites))!=3:
                continue
            support = mask(sites)
            syn,tag = pilot["xor_at"](checks,support),pilot["xor_at"](tags,support)
            union = checks[sites[0]] | checks[sites[1]] | checks[sites[2]]
            d = union
            allowed = []
            while True:
                if (syn.bit_count()-2*(syn & d).bit_count()>=6-nvv and
                    all(((d >> j) & 511).bit_count()<=1 for j in range(0,108,9))):
                    studied += 1
                    if any(t!=tag for t in classes(syn ^ d,d)):
                        hits += 1
                        allowed.append(d)
                if d==0:
                    break
                d = (d-1) & union
            if allowed:
                promising.append((nvv,left,right,nl*nr,tuple(allowed)))
    return promising,studied,hits


def physical_triples(H,checks,tags,run):
    geometry = run_path(str(HERE / "comparar_HGP_225.py"))
    maps = geometry["affine_maps"](checks[:12])
    rows,studied,hits = interior_triples(H,checks,tags)
    candidates = set()
    for nvv,left,right,_,indicators in rows:
        for d in indicators:
            cells = [divmod(j,9) for j in bits(d)]
            def transforms(labels,axis):
                return {(tuple(m[int(i>=nvv)][x] for i,x in enumerate(labels)),
                         tuple(m[axis][c[axis]] for c in cells)) for m in maps}
            for (rows_a,cr),(cols_a,cc) in product(transforms(left,0),transforms(right,1)):
                sites = tuple(sorted(12*r+c if i<nvv else 144+9*r+c
                                     for i,(r,c) in enumerate(zip(rows_a,cols_a))))
                candidates.add((sites,mask(9*r+c for r,c in zip(cr,cc))))
    confirmed = []
    for sites,d in sorted(candidates):
        syn = checks[sites[0]] ^ checks[sites[1]] ^ checks[sites[2]]
        tagsum = tags[sites[0]] ^ tags[sites[1]] ^ tags[sites[2]]
        end,corr = run(syn ^ d)
        if end==d and corr!=tagsum:
            confirmed.append((sites,d))
    return confirmed,{"representantes":292,"representantes_posibles":len(rows),
                      "mascaras_examinadas":studied,"mascaras_posibles":hits,
                      "entradas_fisicas_posibles":len(candidates),"entradas_fisicas_silenciosas":len(confirmed)}


def clean_triple_mass(laws,normalizer,confirmed):
    zero = [law[0][0] for law in laws]
    ratios = defaultdict(int)
    for sites,d in confirmed:
        u,v,w = sites
        ratios[zero[u]*zero[v]*zero[w]] += sum(
            a*b*laws[w][1].get(d ^ x ^ y,0) for x,a in laws[u][1].items() for y,b in laws[v][1].items())
    clean = Q(prod(zero),normalizer**N)
    return clean*sum(Q(n,z) for z,n in ratios.items())


class Choice:
    """Muestreo categórico exacto con 64 bits y desempate racional en las fronteras."""
    def __init__(self,items):
        pairs = [(v,w) for v,w in items if w]
        self.values = [v for v,_ in pairs]
        weights = [w for _,w in pairs]
        if not all(isinstance(x,int) for x in weights):
            from math import lcm
            common = lcm(*(Q(w).denominator for w in weights))
            weights = [int(w*common) for w in weights]
        self.total = sum(weights)
        prefix = 0
        self.exact = []
        for w in weights[:-1]:
            prefix += w
            self.exact.append(prefix << 64)
        self.grid = [n//self.total for n in self.exact]
        self.boundaries = set(self.grid)

    def draw(self,rng):
        if len(self.values)==1:
            return self.values[0]
        x = rng.getrandbits(64)
        if x in self.boundaries:
            index = bisect_right(self.exact,x*self.total+rng.randrange(self.total))
        else:
            index = bisect_right(self.grid,x)
        return self.values[index]


class Background:
    """Subconjunto de Bernoulli independientes. Muestreo disperso y exacto."""
    def __init__(self,probs,first=0):
        self.start = first
        if len(probs)==1:
            self.empty = 1-probs[0]
            self.split = None
            return
        m = len(probs)//2
        self.left,self.right = Background(probs[:m],first),Background(probs[m:],first+m)
        a,b = self.left.empty,self.right.empty
        self.empty = a*b
        self.split = Choice([(1,(1-a)*b),(2,a*(1-b)),(3,(1-a)*(1-b))])

    def occupied(self,rng,out):
        if self.split is None:
            out.append(self.start)
        else:
            side = self.split.draw(rng)
            if side & 1:
                self.left.occupied(rng,out)
            if side & 2:
                self.right.occupied(rng,out)


def trial_inputs(laws,normalizer):
    par0 = [sum(law[0].values()) for law in laws]
    assert len(set(par0[:144]))==len(set(par0[144:]))==1
    p = [Q(normalizer-par0[v],normalizer) for v in [0,224]]
    odds = [x/(1-x) for x in p]
    pzero = (1-p[0])**144*(1-p[1])**81
    type_weights = [[(i,comb(144,i)*comb(81,k-i)*odds[0]**i*odds[1]**(k-i))
                     for i in range(k+1)] for k in range(5)]
    masses = [pzero*sum(w for _,w in row) for row in type_weights]
    threshold = 1-sum(masses)
    assert 0 < threshold < masses[4] < masses[3]
    prob_bg = [1-Q(law[0][0],z) for law,z in zip(laws,par0)]
    tree = Background(prob_bg)
    return {"mass":masses,"tail5":threshold,"types":{k:Choice(type_weights[k]) for k in (3,4)},
            "positive":[Choice(law[1].items()) for law in laws],
            "other":[Choice([(d,n) for d,n in law[0].items() if d]) for law in laws],
            "prob_bg":prob_bg,"tree":tree,"occupied":Choice([(0,tree.empty),(1,1-tree.empty)])}


def sampled_sectors(model,h,checks,tags,run,trials,seed, *, tail_exponent=6):
    """Marca canónica para tres datos y fondo ocupado. Cuatro datos por condicionamiento."""
    if not isinstance(tail_exponent, int) or tail_exponent < 1:
        raise ValueError("El exponente de Chernoff debe ser un entero positivo.")
    read_size = Choice([(j,comb(M,j)*h**j*(1-h)**(M-j)) for j in range(M+1)])
    marker = Choice(enumerate(model["prob_bg"]+[h]*M))
    total_load = sum(model["prob_bg"])+M*h
    out = {}
    for k in (3,4):
        rng = Random(seed+k)
        hits = 0
        eligible = 0
        for _ in range(trials):
            nv = model["types"][k].draw(rng)
            data = rng.sample(range(144),nv)+rng.sample(range(144,225),k-nv)
            marked = marker.draw(rng) if k==3 else -1
            if marked in data:
                continue
            others = []
            if model["occupied"].draw(rng):
                model["tree"].occupied(rng,others)
            dirty = [v for v in others if v not in data and v!=marked]
            j = read_size.draw(rng)
            dirty.extend(N+c for c in rng.sample(range(M),j) if N+c!=marked)
            if marked>=0:
                if any(v<marked for v in dirty):
                    continue
                dirty.append(marked)
            d = syn = logical = 0
            for v in data:
                d ^= model["positive"][v].draw(rng)
                syn ^= checks[v]
                logical ^= tags[v]
            for v in dirty:
                d ^= model["other"][v].draw(rng) if v<N else 1 << (v-N)
            # Una palabra lógica no trivial necesita al menos seis columnas VV.
            if syn.bit_count()-2*(syn & d).bit_count()<6-nv:
                continue
            eligible += 1
            end,corr = run(syn ^ d)
            hits += end==d and corr!=logical
        # Chernoff inferior de un binomial de media mu: P(X<=k)<=exp(-(mu-k)^2/(2mu)).
        # Para cero se usa directamente (1-p)^n. El riesgo por estrato es <= exp(-t).
        t = tail_exponent
        limit = t if hits==0 else hits+t+isqrt(t*t+2*t*hits)+1
        upper = min(Q(1),Q(limit,trials))
        mult = model["mass"][k]*(total_load if k==3 else 1)
        out[k] = (mult*upper,{"ensayos":trials,"semilla":seed+k,"decodificados":eligible,
                              "silenciosos_con_marca_si_corresponde":hits,"limite_exitos":limit})
    return out


def primitive_checks(H,checks,data,anc):
    """Contraste directo de Clifford en paralelo sobre bits de casos distintos."""
    events = [[] for _ in range(L+2)]
    for t in range(L+2):
        for v,row in enumerate(data):
            events[t].append(("CNOT",v,row[t]) if t in row else ("espera_dato",v,-1))
        for c,row in enumerate(anc):
            if t not in row:
                kind = "preparacion" if t==0 else "lectura" if t==L+1 else "espera_ancilla"
                events[t].append((kind,-1,c))
    flat = [event for layer in events for event in layer]
    counts = Counter(event[0] for event in flat)
    assert counts==dict(CNOT=756,espera_dato=1944,espera_ancilla=324,preparacion=108,lectura=108)
    width = [16 if k=="CNOT" else 8 if k=="lectura" else 4 for k,_,_ in flat]
    extra = [defaultdict(int) for _ in flat]
    singles = sum(width)
    rng = Random(3261123)
    for k in range(68):
        sites = range(len(flat)) if k>=64 else rng.sample(range(len(flat)),2+k)
        for i in sites:
            extra[i][rng.randrange(1,width[i])] |= 1 << (singles+k)
    x,z = [0]*(N+M),[0]*(N+M)
    reported = [0]*M
    expected_x,expected_s = [0]*N,[0]*M
    expected_z = [0]*(N+M)
    case = pos = 0
    for t,layer in enumerate(events):
        for kind,v,c in layer:
            if kind=="CNOT":
                x[N+c] ^= x[v]
                z[v] ^= z[N+c]
        for kind,v,c in layer:
            future = [j for time,j in data[v].items() if time>t] if v>=0 else []
            for category in range(width[pos]):
                test = (1 << case) | extra[pos].get(category,0)
                case += 1
                if kind=="CNOT":
                    pc,pa = divmod(category,4)
                    flip = False
                elif v>=0:
                    pc,pa,flip = category,0,False
                elif kind=="lectura":
                    pa,flip = divmod(category,2)
                    pc = 0
                else:
                    pc,pa,flip = 0,category,False
                # Orden de Paulis I,X,Y,Z. Trayectoria física por conjugación.
                if v>=0:
                    if pc in (1,2):
                        x[v] ^= test
                    if pc in (2,3):
                        z[v] ^= test
                if c>=0:
                    if pa in (1,2):
                        x[N+c] ^= test
                    if pa in (2,3):
                        z[N+c] ^= test
                    if flip:
                        reported[c] ^= test
                # Fórmula cerrada: sólo las CNOT posteriores leen el X del control.
                if pc in (1,2):
                    expected_x[v] ^= test
                    for j in future:
                        expected_s[j] ^= test
                if (pa in (1,2))!=bool(flip):
                    expected_s[c] ^= test
                if pc in (2,3):
                    expected_z[v] ^= test
                if pa in (2,3):
                    expected_z[N+c] ^= test
                    for time,wire in anc[c].items():
                        if time>t:
                            expected_z[wire] ^= test
            pos += 1
    assert x[:N]==expected_x and z==expected_z and [x[N+c]^reported[c] for c in range(M)]==expected_s
    assert singles==22464 and pos==3240
    fibre_sizes = [13]*N+[5]*M
    neighbourhood = [{i} for i in range(N+M)]
    for row in H:
        for line in checks[:12]:
            rect = mask(9*i+c for i in bits(row) for c in bits(line))
            readers = {v for v,s in enumerate(checks) if s & rect} | {N+c for c in bits(rect)}
            for v in readers:
                neighbourhood[v].update(readers)
    degrees = [sum(fibre_sizes[v] for v in row)-1 for row in neighbourhood]
    doubled = sum(k*d for k,d in zip(fibre_sizes,degrees))
    assert sum(fibre_sizes)==pos+N==3465 and doubled%2==0
    return {"capas":[len([v for v in data if t in v]) for t in range(1,L+1)],
            "localizaciones":dict(counts),"N_loc":pos,"vertices_host":pos+N,
            "aristas_host":doubled//2,"grado_host":max(degrees),"categorias_unitarias":singles,
            "asignaciones_simultaneas":68}


def confirm_laws(data,laws,norm,rc):
    """Caracteres de (E,s) desde las inserciones Pauli, frente a la convolución de (E,D)."""
    den,k = 15*rc.denominator,rc.numerator
    for v,calendar in enumerate(data):
        full = mask(calendar.values())
        sites = sorted(calendar.values())
        for tmask in range(1 << len(sites)):
            t = mask(sites[j] for j in bits(tmask))
            for u in (0,1):
                def char(e,s):
                    return -1 if ((u*e+(t & s).bit_count())%2) else 1
                expected = 1
                for layer in range(L+2):
                    future = mask(c for time,c in calendar.items() if time>layer)
                    if layer in calendar:
                        here = 1 << calendar[layer]
                        expected *= den-12*k+4*k*(char(1,future)+char(0,here)+char(1,future^here))
                    else:
                        expected *= den-10*k+10*k*char(1,future)
                actual = sum(n*char(e,d ^ (full if e else 0))
                             for e in (0,1) for d,n in laws[v][e].items())
                assert actual==expected and abs(actual)<=norm


def decimal_bound(q,rounding):
    with localcontext() as ctx:
        ctx.prec,ctx.rounding = 15,rounding
        return format(Decimal(q.numerator)/Decimal(q.denominator),".14E")


def certificate(l,u):
    assert 0<l<=u<=1
    point = decimal_bound(2*l*u/(l+u),ROUND_HALF_EVEN)
    eps = max(Q(point)/l-1,1-Q(point)/u)
    return {"L":decimal_bound(l,ROUND_FLOOR),"U":decimal_bound(u,ROUND_CEILING),
            "estimacion":point,"error_relativo_pct_sup":decimal_bound(100*eps,ROUND_CEILING)}


def uniform_rectangles():
    base = run_path(str(HERE/"verificar_H5A_grados_3_4.py"))
    mat,down,causes = base["typed_counts"]()
    winners = base["local_signatures"]()
    alpha = tuple(sum(r) for r in base["first_extinction"](winners,False)[1].values())
    sigma = tuple(1+sum(r) for r in base["first_extinction"](winners,True)[1].values())
    assert alpha==(Q(7,3),3) and sigma==(Q(17,3),Q(103,15))
    witnesses = (
        ("2e-13","1e-6",".006064","4.766674047e-7",".00001516720759",".000007395863118",".000001039067636",".0007047078695"),
        ("1.1e-13","1e-4",".00521","4.582013048e-7",".00001683579891",".000007635823356",".0001060205702",".001384281556"),
        ("4.5e-16","1e-3",".001735","1.966675993e-7",".00002330732619",".000007438300174",".00132367094",".003980000608"),
    )
    results = []
    dot = lambda a,b: sum(x*y for x,y in zip(a,b))
    ex = base["exp_upper"]
    for values in witnesses:
        rc,rm,y,a,b,g,m,v = map(Q,values)
        odds_data = (1-rc)**(-12)-1
        odds_anc = (1-rc)**(-5)/(1-rm)-1
        R = (a+3*g,b+4*g,0,m)
        bv = max(dot(row,R) for row in mat[:2])
        bd = dot(mat[3],R)
        hv,hc,hg,hm = [dot(row,(a,b,g,m))+k*v for row,k in zip(mat,down)]
        tau = Q(1,100)
        ratios = ((odds_data*ex(sigma[0]*bv+tau+hv)/a)**3/y**7,
                  odds_data*ex(sigma[1]*bv+tau+hc)/(b*y**3),
                  y*ex(hg)*dot(causes,(a,b,g))/g,odds_anc*ex(bd+tau+hm)/m,
                  y*(ex(hg)-1)/v)
        assert max(R)<1 and all(q<Q(".995") for q in ratios)
        results.append({"rc_max":values[0],"rm_max":values[1],"y":values[2],
                        "u_VV_CC_P_C_v":values[3:],"mayor_cociente_sup":decimal_bound(max(ratios),ROUND_CEILING)})
    return results


def labelled_transfer():
    """Certificado uniforme de núcleos con etiquetas de circuito en las fuentes."""
    def previous(name):
        return run_path(str(HERE/("verificar_H5A_3_4_"+name+".py")))
    passive = previous("fondo_pasivo")
    first,adjacent = previous("medidas_y_un_dato"),previous("pares_adyacentes")
    two,corners = previous("dos_componentes"),previous("esquinas")
    profiles = adjacent["stripe_profiles"]()
    ranks,costs = passive["step_costs"](profiles,11)
    coverage = two["local_checks"](first["budgets"])
    # Un grupo de A contiene las componentes de E que alojan sus datos reales.
    # En un grupo grande, también los costes de pares y de un solo dato caben
    # en la tarifa lineal. Una fuente de paridad física cero sólo añade un sitio.
    for (i,j),(ordinary,overlap) in zip(((0,0),(0,1),(1,1)),costs):
        assert 1<=min(ordinary,overlap) and max(ordinary,overlap)<=ranks[i]+ranks[j]
    passive["certificate"](first["exp_upper"],corners["endpoint_counts"],coverage,ranks,costs)
    rc,rm = Q("2.9e-10"),Q(".001")
    data_odds = (1-Q(2,3)*rc)**(-8)*(1-Q(4,5)*rc)**(-4)-1
    check_odds = (1-Q(2,3)*rc)**(-5)/(1-rm)-1
    p0,q0 = Q("2.5e-9"),Q(".001")
    assert data_odds < p0/(1-2*p0) and check_odds < q0/(1-2*q0)
    return {"rc_max":"2.9e-10","rm_max":"1e-3","inclinacion":".01",
            "odds_dato_sup":decimal_bound(data_odds,ROUND_CEILING),
            "odds_check_sup":decimal_bound(check_odds,ROUND_CEILING),
            "cobertura":[int(v) for v in coverage],"rangos":[str(x) for x in ranks],
            "costes_vv_vc_cc":[[str(x) for x in pair] for pair in costs]}


def contrast_labelled_cores(H,checks,tags):
    """Proyección de trayectorias filtradas conservando las marcas de paridad cero."""
    choices = all_winners(H,checks,tags)
    reps = {}
    for r,row in enumerate(H):
        for j,line in enumerate(checks[:12]):
            qubits = [12*i+j for i in bits(row)]+[144+9*r+c for c in bits(line)]
            for f in range(1,128):
                action = mask(qubits[i] for i in bits(f))
                phi = pilot["xor_at"](checks,action)
                if phi and 2*phi.bit_count()>=3*f.bit_count():
                    reps[phi] = min(reps.get(phi,(0,0)),(-f.bit_count(),action))

    @lru_cache(maxsize=8192)
    def trace(s):
        wins = choices(s)
        if not wins:
            return ()
        phi,_ = min(wins,key=lambda x:reps[x[0]])
        return ((phi,reps[phi][1]),)+trace(s^phi)

    rects = [mask(9*i+c for i in bits(row) for c in bits(line)) for row in H for line in checks[:12]]
    readers = [mask(g for g,p in enumerate(rects) if p & col) for col in checks]
    wires = (0,1,12,144)
    options = [[None,(0,col & -col),(1,0),(1,col^(col & -col))] for col in map(checks.__getitem__,wires)]
    counts = Counter()
    for labels in product(*options):
        active = {v:a for v,a in zip(wires,labels) if a is not None}
        source = mask(active)
        data = mask(v for v,(e,_) in active.items() if e)
        intrinsic = 0
        original = 0
        for v,(_,d) in active.items():
            intrinsic ^= d
            original |= checks[v]
        for bmask in range(16):
            b = mask((0,9,20,30)[i] for i in bits(bmask))
            d = intrinsic ^ b
            s = pilot["xor_at"](checks,data)^d
            path = trace(s)
            counts["entradas"] += 1
            terminal = s
            for phi,_ in path:
                terminal ^= phi
            if terminal!=d:
                continue
            counts["filtradas"] += 1
            touched,coverage = source,original
            for phi,action in path:
                coverage |= phi
                touched |= action
            passive = b & ~coverage
            assert trace(s^passive)==path
            remaining = set(bits(touched))
            parts = 0
            while remaining:
                part = {min(remaining)}
                while True:
                    reach = 0
                    for v in part:
                        reach |= readers[v]
                    extended = {v for v in remaining if readers[v] & reach}
                    if extended==part:
                        break
                    part = extended
                members,ports,d0 = mask(part),0,0
                assert members & source
                for v in part:
                    ports |= checks[v]
                    if v in active:
                        d0 ^= active[v][1]
                d0 ^= b & coverage & ports
                inner = pilot["xor_at"](checks,data & members)^d0
                local = tuple(step for step in path if step[1] & members)
                assert trace(inner)==local and trace(inner^passive)==local
                remaining -= part
                parts += 1
            if parts>=2:
                counts["varios_nucleos"] += 1
            if any(e==0 for e,_ in active.values()):
                counts["con_fuente_de_paridad_cero"] += 1
    assert counts["entradas"]==4096 and min(counts.values())>0 and len(counts)==4
    return dict(counts)


def main(output):
    start = perf_counter()
    H,checks,tags,_,_,run = pilot["build"]()
    data,anc = schedule(checks)
    assert sum(Q(6)**j/factorial(j) for j in range(20))>400
    contract = primitive_checks(H,checks,data,anc)
    rectangles = uniform_rectangles()
    transfer = labelled_transfer()
    contrast = contrast_labelled_cores(H,checks,tags)
    patterns = []
    pilot["row_sectors"](checks,tags,run,patterns=patterns)
    # Para una pareja, la medida interior indispensable está fuera de ambas franjas.
    assert len(patterns)==324 and all(d and not (checks[u] | checks[v]) & d for _,u,v,d in patterns)
    previous = run_path(str(HERE/"comparar_HGP_225.py"))
    previous["verify_exterior_bit"]()
    previous["physical_actions"](H,checks,tags)
    confirmed,census = physical_triples(H,checks,tags,run)
    assert (census["representantes_posibles"],len(confirmed))==(5,1656)
    assert sum(d==0 for _,d in confirmed)==36
    prepared = perf_counter()
    results = []
    fixed_trials = 120000
    points = [("1e-6",["0","5e-7","1e-6"]),("3.75e-7",["1e-3"]),("1.5e-9",["1e-3"])]
    for rc_text,read_rates in points:
        rc = Q(rc_text)
        laws,normalizer = joint_fibres(data,rc)
        confirm_laws(data,laws,normalizer,rc)
        spectra = pair_spectra(laws,patterns,normalizer)
        triples = clean_triple_mass(laws,normalizer,confirmed)
        model = trial_inputs(laws,normalizer)
        contrast_pair(laws,normalizer,patterns[0],ancilla_rate(rc,Q(read_rates[0])))
        for rm_text in read_rates:
            h = ancilla_rate(rc,Q(rm_text))
            l2,u2 = pair_bounds(spectra,h)
            l3 = triples*(1-h)**M
            mass = model["mass"]
            burden = sum(model["prob_bg"])+M*h
            l = l2+l3
            base_u = u2+l3+model["tail5"]
            deterministic = certificate(l,base_u+mass[4]+mass[3]*min(1,burden))
            row = {"rc":rc_text,"rm":rm_text,"par_L":decimal_bound(l2,ROUND_FLOOR),
                   "par_U":decimal_bound(u2,ROUND_CEILING),"trios_fondo_nulo":decimal_bound(l3,ROUND_FLOOR),
                   "P_exactamente_3_sup":decimal_bound(mass[3],ROUND_CEILING),
                   "P_exactamente_4_sup":decimal_bound(mass[4],ROUND_CEILING),
                   "P_al_menos_5_sup":decimal_bound(model["tail5"],ROUND_CEILING),
                   "carga_marca_sup":decimal_bound(burden,ROUND_CEILING),"determinista":deterministic}
            if rc_text!="1.5e-9":
                strata = sampled_sectors(model,h,checks,tags,run,fixed_trials,2026112300+100*len(results))
                upper = min(base_u+sum(x[0] for x in strata.values()),
                            base_u+mass[4]+mass[3]*min(1,burden))
                row["confianza_conjunta_98_pct"] = certificate(l,upper)
                row["muestreo"] = {str(k):v[1] for k,v in strata.items()}
            results.append(row)
            interval = row.get("confianza_conjunta_98_pct",deterministic)
            print(f"rc={rc_text}, rm={rm_text}, error (%) {interval['error_relativo_pct_sup']}",flush=True)
    data = {"fecha":"2026-11-23","python":sys.version,"sistema":platform.platform(),
            "contrato":contract,"minimo_fallos_primitivos_piloto":3,"uniforme_nucleos":transfer,
            "uniforme_arbol_sin_compresion":rectangles,"contraste_nucleos":contrast,"censo_trios":census,
            "puntos":results,"preparacion_s":prepared-start,"total_interno_s":perf_counter()-start,
            "pico_MiB":pilot["peak_memory"]()/2**20}
    if output:
        output.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf8")
    else:
        print(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False))


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar con las comprobaciones de Python activas.")
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--salida",type=Path,help="Guardar las cotas y las mediciones en JSON.")
    arguments = parser.parse_args()
    main(arguments.salida)
