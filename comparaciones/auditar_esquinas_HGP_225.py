"""Contraste independiente de los once intervalos del HGP afín [[225,9,6]]."""
import argparse
import json
from collections import defaultdict
from decimal import Decimal
from fractions import Fraction as F
from functools import lru_cache, reduce
from itertools import combinations, permutations, product
from math import comb, factorial, isqrt, prod
from operator import xor
from pathlib import Path
from runpy import run_path

N, M = 225, 108


def positions(word):
    while word:
        bit = word & -word
        yield bit.bit_length() - 1
        word -= bit


def mask(sites):
    return sum(1 << v for v in sites)


def bit_sum(words):
    return reduce(xor, words, 0)


def incidence():
    points = list(product(range(3), repeat=2))
    lines = [tuple(i for i, (x, y) in enumerate(points) if (y-m*x) % 3 == b)
             for m, b in product(range(3), repeat=2)]
    lines += [tuple(i for i, (x, y) in enumerate(points) if x == b) for b in range(3)]
    through = [tuple(j for j, line in enumerate(lines) if a in line) for a in range(9)]
    checks = [mask(9*i+a for a in lines[j]) for i, j in product(range(12), repeat=2)]
    checks += [mask(9*i+b for i in through[a]) for a, b in product(range(9), repeat=2)]
    generators = [tuple([12*i+j for i in through[a]] + [144+9*a+b for b in lines[j]])
                  for a, j in product(range(9), range(12))]
    kernel = [word.bit_count() for word in range(1 << 12)
              if not bit_sum(mask(lines[i]) for i in positions(word))]
    assert sorted(kernel) == [0]+[6]*6+[12]
    return lines, checks, generators


def decoder(checks, generators):
    # Cociente directo por las 108 filas estabilizadoras. En síndrome físico cero,
    # un resto no nulo equivale a una clase lógica X no trivial.
    pivots = {}
    for sites in generators:
        v = mask(sites)
        while v and v.bit_length() in pivots:
            v ^= pivots[v.bit_length()]
        if v:
            pivots[v.bit_length()] = v
    assert len(pivots) == M
    tags = []
    for site in range(N):
        v = 1 << site
        for n in sorted(pivots, reverse=True):
            if v & (1 << (n-1)):
                v ^= pivots[n]
        tags.append(v)
    best = {}
    for sites in generators:
        assert bit_sum(checks[v] for v in sites) == bit_sum(tags[v] for v in sites) == 0
        for k in range(1, len(sites)):
            for subset in combinations(sites, k):
                phi = bit_sum(checks[v] for v in subset)
                if not phi or 2*phi.bit_count() < 3*k:
                    continue
                action = (-k, mask(subset), bit_sum(tags[v] for v in subset))
                assert phi not in best or best[phi][2] == action[2]
                best[phi] = min(best.get(phi, action), action)
    entries = sorted(best.items(), key=lambda item: item[1])
    assert len(entries) == 6273
    footprints = [p for p, _ in entries]
    classes = [a[2] for _, a in entries]
    touches = [0] * M
    by_size = defaultdict(int)
    for j, phi in enumerate(footprints):
        by_size[phi.bit_count()] |= 1 << j
        for c in positions(phi):
            touches[c] |= 1 << j
    scores = defaultdict(list)
    for n, group in by_size.items():
        for h in range(n+1):
            if 100*(2*h-n) >= n:
                scores[F(h, n)].append((h, group))
    priorities = [scores[r] for r in sorted(scores, reverse=True)]
    universe = (1 << len(entries))-1

    @lru_cache(maxsize=8192)
    def best_moves(s):
        # Los planos binarios guardan |s intersección huella| de cada
        # acción física. Se comparan globalmente sin tablas de rectángulos.
        planes = []
        for c in positions(s):
            carry = touches[c]
            for i, plane in enumerate(planes):
                planes[i] ^= carry
                carry &= plane
                if not carry:
                    break
            if carry:
                planes.append(carry)
        eq = {}
        for group in priorities:
            chosen = 0
            for hits, candidates in group:
                if hits not in eq:
                    a = universe if hits < 1 << len(planes) else 0
                    for i, plane in enumerate(planes):
                        a &= plane if hits & (1 << i) else ~plane
                    eq[hits] = a
                chosen |= eq[hits] & candidates
            if chosen:
                return chosen
        return 0

    @lru_cache(maxsize=24576)
    def run(s):
        winners = best_moves(s)
        if not winners:
            return s, 0
        i = (winners & -winners).bit_length()-1
        end, tag = run(s ^ footprints[i])
        return end, tag ^ classes[i]

    @lru_cache(maxsize=32768)
    def possible(s, target, tag):
        if s.bit_count() <= target.bit_count():
            return s == target and tag != 0 and not best_moves(s)
        return any(possible(s ^ footprints[i], target, tag ^ classes[i])
                   for i in positions(best_moves(s)))

    return tags, run, possible


def affine_frames(lines):
    # Las imágenes de tres puntos no colineales determinan la aplicación.
    lookup = {frozenset(line): i for i, line in enumerate(lines)}
    actions = []
    for o, ex, ey in permutations(range(9), 3):
        if any({o, ex, ey}.issubset(line) for line in lines):
            continue
        perm = []
        for x, y in product(range(3), repeat=2):
            coord = [((1-x-y)*f(o)+x*f(ex)+y*f(ey)) % 3
                     for f in (lambda z: z//3, lambda z: z%3)]
            perm.append(3*coord[0]+coord[1])
        on_lines = tuple(lookup[frozenset(perm[p] for p in line)] for line in lines)
        actions.append((on_lines, tuple(perm)))
    assert len(set(actions)) == 432
    return actions


def axis_orbits(types, group):
    todo = set(product(*(range(12 if t == 0 else 9) for t in types)))
    while todo:
        rep = min(todo)
        images = {tuple(g[t][a] for t, a in zip(types, rep)) for g in group}
        assert images <= todo
        todo -= images
        yield rep, len(images)


def physical(types, rows, cols):
    return tuple(12*i+j if t == 0 else 144+9*i+j for t, i, j in zip(types, rows, cols))


def censes(lines, checks, tags, run, possible):
    pairs = []
    for i in range(12):
        for u, v in combinations(range(12*i, 12*i+12), 2):
            for d in [0]+[1 << (9*i+a) for a in range(9)]:
                end, result = run(checks[u] ^ checks[v] ^ d)
                if end == d and result != tags[u] ^ tags[v]:
                    assert d and not d & (checks[u] | checks[v])
                    pairs.append((i, u, v, d))
                    for c in range(M):
                        if c//9 != i:
                            extra = d ^ (1 << c)
                            stop, corr = run(checks[u] ^ checks[v] ^ extra)
                            assert stop == extra and corr != tags[u] ^ tags[v]

    group = affine_frames(lines)
    total_sites = one_readout = reps = 0
    candidates = set()
    for nvv in range(4):
        types = (0,)*nvv + (1,)*(3-nvv)
        axes = list(axis_orbits(types, group))
        unordered = factorial(nvv)*factorial(3-nvv)
        count = one = 0
        for (rs, nr), (cs, nc) in product(axes, repeat=2):
            sites = physical(types, rs, cs)
            if len(set(sites)) < 3:
                continue
            reps += 1
            count += nr*nc
            a, b, c = (checks[v] for v in sites)
            syn = a ^ b ^ c
            tag = bit_sum(tags[v] for v in sites)
            one += nr*nc*sum(possible(syn ^ (1 << j), 1 << j, tag) for j in range(M))
            # D interior estable: en cada fila se elige cero o un punto.
            by_row = defaultdict(list)
            for j in positions(a | b | c):
                by_row[j//9].append(1 << j)
            for selection in product(*([0]+row for row in by_row.values())):
                d = sum(selection)
                if not possible(syn ^ d, d, tag):
                    continue
                cells = [divmod(j, 9) for j in positions(d)]
                def images(axis, labels):
                    return {(tuple(g[t][x] for t, x in zip(types, labels)),
                             tuple(g[axis][cell[axis]] for cell in cells)) for g in group}
                for (r1, rr), (c1, cc) in product(images(0, rs), images(1, cs)):
                    site_key = tuple(sorted(physical(types, r1, c1)))
                    candidates.add((site_key, mask(9*r+c for r, c in zip(rr, cc))))
        assert count % unordered == one % unordered == 0
        total_sites += count//unordered
        one_readout += one//unordered
    triples = []
    for sites, d in sorted(candidates):
        s = bit_sum(checks[v] for v in sites)
        end, tag = run(s ^ d)
        if end == d and tag != bit_sum(tags[v] for v in sites):
            triples.append((sites, d))
    assert total_sites == comb(N, 3) and reps == 292
    return pairs, triples, one_readout, len(candidates)



def binomial_counts(n, prob):
    a, b = prob.numerator, prob.denominator
    return [comb(n, i)*a**i*(b-a)**(n-i) for i in range(n+1)], b**n


def bernoulli_point(p, q, pairs, triple0, triple1):
    ad, den_d = binomial_counts(N, p)
    bs, den_m = binomial_counts(M, q)
    # Condicionado al tamaño, cada soporte de datos y de lecturas es equiprobable.
    pairs_p = F(ad[2]*pairs, den_d*comb(N, 2))
    trios_p = F(ad[3], den_d*comb(N, 3))
    lower = pairs_p*(F(bs[1], M)+F(99*bs[2], comb(M, 2)))/den_m
    lower += trios_p*triple0*F(bs[0], den_m)
    any_exterior = sum(F(bs[b]*comb(99, b-1), comb(M, b)) for b in range(1, 101))
    upper = pairs_p*any_exterior/den_m
    upper += trios_p*(triple0*bs[0]+F(triple1*bs[1], M)+comb(N, 3)*sum(bs[2:]))/den_m
    return lower, upper+F(sum(ad[4:]), den_d)


def greedy_calendar(checks, cc_first):
    wires = sorted(range(N), key=lambda v: (int((v < 144) == cc_first), v))
    busy = [0]*N
    gates = [[] for _ in wires]
    layers = [0]*10
    for j in range(M):
        occupied = 0
        for v in wires:
            if not (checks[v] & (1 << j)):
                continue
            free = 1023 & ~(busy[v] | occupied)
            slot = free & -free
            assert slot
            t = slot.bit_length()
            gates[v].append((t, j))
            layers[t-1] += 1
            occupied |= slot
            busy[v] |= slot
        assert occupied.bit_count() == 7
    return [sorted(row) for row in gates], layers


def circuit_laws(calendar, rc, rm):
    # Estado cronológico (X del dato, bits brutos de ancilla). El CNOT ideal
    # lee el control antes de las Paulis posteriores. D se forma al final.
    unit = 15*rc.denominator
    k = rc.numerator
    wait = defaultdict(int, {(0, 0): unit-15*k})
    gate = defaultdict(int, {(0, 0): unit-15*k})
    for pauli in range(1, 4):
        wait[int(pauli in (1, 2)), 0] += 5*k
    for a, b in product(range(4), repeat=2):
        if a or b:
            gate[int(a in (1, 2)), int(b in (1, 2))] += k

    @lru_cache(None)
    def local(times):
        states = {(0, 0): 1}
        for t in range(12):
            bit = 1 << times.index(t) if t in times else 0
            updated = defaultdict(int)
            for (e, raw), prob in states.items():
                recorded = raw ^ (bit if e else 0)
                for (flip, flag), w in (gate if bit else wait).items():
                    updated[e ^ flip, recorded ^ (bit if flag else 0)] += prob*w
            states = updated
        full = (1 << len(times))-1
        return [[states[e, d ^ (full if e else 0)] for d in range(full+1)] for e in (0, 1)]

    laws = []
    norm = unit**12
    for row in calendar:
        sites = [j for _, j in row]
        weights = local(tuple(t for t, _ in row))
        physical_masks = [mask(sites[i] for i in positions(d)) for d in range(1 << len(row))]
        laws.append([dict(zip(physical_masks, sub)) for sub in weights])
        assert sum(map(sum, weights)) == norm
        marginal = (1-(1-16*rc/15)**len(row)*(1-4*rc/3)**(12-len(row)))/2
        assert F(sum(weights[1]), norm) == marginal
    h = rm
    for _ in range(5):
        h += (1-2*h)*2*rc/3
    return laws, norm, h


PREC = 208
GRID = 1 << PREC


def convolve(dist, weights, den):
    # Convolución XOR positiva de hasta nueve bits. El suelo de cada peso
    # y de cada celda pierde a lo sumo (16+512)/GRID en norma l1 por paso.
    scaled = [(j, (w*GRID)//den) for j, w in weights.items() if w]
    new = defaultdict(int)
    for d, mass in dist.items():
        for shift, weight in scaled:
            new[d ^ shift] += mass*weight
    return {d: m >> PREC for d, m in new.items() if m >= GRID}


def pair_mass(pairs, laws, norm, h):
    rows = [defaultdict(list) for _ in range(12)]
    for i, u, v, d in pairs:
        rows[i][u, v].append(d >> (9*i))
    totals = [0, 0]
    for i, targets in enumerate(rows):
        shift = 9*i
        area = 511 << shift
        inside = range(12*i, 12*i+12)
        backgrounds = []
        for free in (False, True):
            dist = {0: GRID}
            for v, law in enumerate(laws):
                if v in inside:
                    continue
                projected = defaultdict(int)
                for d, w in law[0].items():
                    if free or d & ~area == 0:
                        projected[(d & area) >> shift] += w
                dist = convolve(dist, projected, norm)
            for bit in range(9):
                dist = convolve(dist, {0: h.denominator-h.numerator, 1 << bit: h.numerator},
                                h.denominator)
            backgrounds.append(dist)
        for pair, accepted in targets.items():
            dist = {0: GRID}
            for v in inside:
                projected = {d >> shift: w for d, w in laws[v][int(v in pair)].items()}
                dist = convolve(dist, projected, norm)
            for mode, other in enumerate(backgrounds):
                totals[mode] += sum(w*other.get(d ^ target, 0)
                                    for d, w in dist.items() for target in accepted)
    # 225 cables + nueve lecturas, más el producto de dos distribuciones
    # subnormalizadas (exacto en GRID^2). La pérdida se suma sobre los patrones.
    loss = F(len(pairs)*234*(512+16), GRID)
    exterior = (1-h)**99+99*h*(1-h)**98
    return F(totals[0], GRID**2)*exterior, F(totals[1], GRID**2)+loss


def controlled_triples(triples, laws, norm, h):
    all_zero = prod(law[0][0] for law in laws)

    @lru_cache(None)
    def complement(removed):
        return all_zero//removed

    numerator = 0
    for (u, v, w), target in triples:
        a, b, c = (laws[j] for j in (u, v, w))
        inner = sum(i*j*c[1].get(x ^ y ^ target, 0)
                    for x, i in a[1].items() for y, j in b[1].items())
        numerator += complement(a[0][0]*b[0][0]*c[0][0])*inner
    return F(numerator, norm**N)*(1-h)**M


def parity_counts(laws, norm):
    odds = [F(sum(law[1].values()), norm) for law in laws]
    assert len(set(odds[:144])) == len(set(odds[144:])) == 1
    a, da = binomial_counts(144, odds[0])
    b, db = binomial_counts(81, odds[-1])
    masses = [F(sum(a[i]*b[k-i] for i in range(k+1)), da*db) for k in range(5)]
    tail = F(sum(w*sum(b[max(0, 5-i):]) for i, w in enumerate(a)), da*db)
    assert sum(masses)+tail == 1
    return masses, tail


def printed_interval(lo, hi, record):
    def ulp(s):
        d = Decimal(s)
        return F(10)**d.as_tuple().exponent
    low, high, point = (F(record[k]) for k in ('L', 'U', 'estimacion'))
    # El cálculo Walsh original añade 324*512/2^128 en cada borde del par.
    original_slack = F(324*512, 1 << 128)
    assert 0 <= lo-low <= ulp(record['L'])+2*original_slack
    assert 0 <= high-hi <= ulp(record['U'])+2*original_slack
    assert abs(point - 2*lo*hi/(lo+hi)) <= ulp(record['estimacion'])+4*original_slack
    assert 100*max(point/lo-1, 1-point/hi) <= F(record['error_relativo_pct_sup'])


def verify_confidence(trials):
    # El límite es monótono en todos los conteos. Al invertir la cola de
    # Chernoff en el caso x>0 el riesgo es <= e^-7; para x=0 basta (1-p)^n.
    previous = 7
    for x in range(1, trials+1):
        upper = x+7+isqrt(49+14*x)+1
        assert upper >= previous and (upper-x)**2 >= 14*upper
        previous = upper
    assert 16/sum(F(7)**j/factorial(j) for j in range(20)) < F('.0146')


def check_proposal(model, laws, norm, masses):
    def compare(choice, weights):
        # U = 2^64 * total equiprobable microceldas. El caso de frontera
        # se desempata con el entero independiente uniforme del muestreador.
        assert all(x % (1 << 64) == 0 for x in choice.exact)
        assert choice.grid == [x//choice.total for x in choice.exact]
        assert choice.boundaries == set(choice.grid)
        cumulative = [0]+[x >> 64 for x in choice.exact]+[choice.total]
        counts = {v: b-a for v, a, b in zip(choice.values, cumulative, cumulative[1:])}
        weights = {v: w for v, w in weights.items() if w}
        total = sum(weights.values())
        assert counts.keys() == weights.keys()
        assert all(n*total == weights[v]*choice.total for v, n in counts.items())

    marginals = []
    for j, (zero, positive) in enumerate(laws):
        p = F(sum(positive.values()), norm)
        marginals.append(p)
        assert model['prob_bg'][j] == 1-F(zero[0], sum(zero.values()))
        compare(model['positive'][j], positive)
        compare(model['other'][j], {d: w for d, w in zero.items() if d})
    a, da = binomial_counts(144, marginals[0])
    b, db = binomial_counts(81, marginals[-1])
    for k in (3, 4):
        weights = {i: a[i]*b[k-i] for i in range(k+1)}
        assert F(sum(weights.values()), da*db) == masses[k]
        compare(model['types'][k], weights)

    def check_tree(node):
        if node.split is None:
            assert node.empty == 1-model['prob_bg'][node.start]
        else:
            left, right = node.left.empty, node.right.empty
            assert node.empty == left*right
            compare(node.split, {1:(1-left)*right, 2:left*(1-right), 3:(1-left)*(1-right)})
            check_tree(node.left)
            check_tree(node.right)
    check_tree(model['tree'])
    empty = model['tree'].empty
    compare(model['occupied'], {0:empty, 1:1-empty})


def replay_point(record, laws, norm, h, checks, tags, run, pairs, triples, sampler, index):
    l2, u2 = pair_mass(pairs, laws, norm, h)
    clean = controlled_triples(triples, laws, norm, h)
    masses, tail = parity_counts(laws, norm)
    load = M*h + sum(1-F(law[0][0], sum(law[0].values())) for law in laws)
    lo = l2+clean
    base = u2+clean+tail
    det = base+masses[3]*min(1, load)+masses[4]
    printed_interval(lo, det, record['determinista'])
    # Para reproducir la MISMA muestra se conserva el generador de propuestas
    # y sus semillas. Las leyes de entrada y el decodificador son los de esta auditoría.
    model = sampler['trial_inputs'](laws, norm)
    assert model['mass'] == masses and model['tail5'] == tail
    check_proposal(model, laws, norm, masses)
    seed = 2026112400+100*index
    n = 120_000
    samples = sampler['sampled_sectors'](model, h, checks, tags, run, n, seed, tail_exponent=7)
    rest = 0
    for k, multiplier in ((3, load), (4, 1)):
        observed = record['estratos'][str(k)]
        assert samples[k][1] == observed
        x = observed['silenciosos_con_marca_si_corresponde']
        u = 7 if x == 0 else x+7+isqrt(49+14*x)+1
        assert observed['limite_exitos'] == u and observed['ensayos'] == n
        rest += masses[k]*multiplier*min(1, F(u, n))
    printed_interval(lo, min(det, base+rest), record['intervalo_conjunto'])


def main(folder):
    result = json.loads((folder/'piloto_HGP_225_esquinas_20261124.json').read_text(encoding='utf8'))
    lines, checks, gen = incidence()
    tags, run, possible = decoder(checks, gen)
    pairs, triples, one, possible_clean = censes(lines, checks, tags, run, possible)
    triple0 = sum(d == 0 for _, d in triples)
    assert (len(pairs), triple0, one, possible_clean, len(triples)) == (
        result['patrones_de_pareja'], result['trios_Bernoulli_cero_lecturas'],
        result['trios_Bernoulli_una_lectura_max'], result['censo_trios_circuito']['entradas_fisicas_posibles'],
        result['censo_trios_circuito']['entradas_fisicas_silenciosas'])
    print(f'Censo: {len(pairs)} patrones, {len(pairs)*99} extensiones, {triple0} tríos sin lecturas, '
          f'{one} con una lectura como máximo, {len(triples)} configuraciones limpias.', flush=True)
    phenom = [('5.40e-7','1e-3'),('1e-6','1e-4'),('1.06e-6','1e-5')]
    assert len(result['Bernoulli']) == len(phenom)
    for (p, q), row in zip(phenom, result['Bernoulli']):
        assert (F(p), F(q)) == (F(row['p']), F(row['q']))
        lo, hi = bernoulli_point(F(p), F(q), len(pairs), triple0, one)
        printed_interval(lo, hi, row)
    print('3/3 intervalos Bernoulli por masas binomiales positivas.', flush=True)
    verify_confidence(120_000)
    sampler = run_path(str(folder/'circuito_HGP_225.py'))
    circuit = [('1e-3','.015','1.005e-7','9.18e-8'),('1e-4','.015','1.85e-7','1.68e-7'),
               ('1e-5','.015','1.975e-7','1.781e-7'),('1e-5','.0025','2e-7','1.803e-7')]
    assert len(result['circuito']) == 2*len(circuit)
    for idx, z in enumerate(result['circuito']):
        m, tilt, vv, cc = circuit[idx//2]
        rc = F(cc if idx % 2 else vv)
        order = 'CC primero' if idx % 2 else 'VV primero'
        assert (F(z['rc']), F(z['rm']), F(z['lambda_del_certificado']), z['orden_ejecutado']) == (
            rc, F(m), F(tilt), order)
        gates, counts = greedy_calendar(checks, cc_first=bool(idx % 2))
        assert counts == result['calendarios'][order]['capas'] and sum(counts) == 756
        laws, norm, h = circuit_laws(gates, rc, F(m))
        replay_point(z, laws, norm, h, checks, tags, run, pairs, triples, sampler, idx)
        print(f"Circuito {idx+1}/8, {order}, r_c={z['rc']}, r_m={m}: ambos intervalos y muestra conformes.", flush=True)
    print('11/11. Los 16 estratos reproducen sus conteos. Riesgo conjunto < 0.0146.')


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Este verificador requiere las aserciones activas.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--carpeta', type=Path, default=Path(__file__).resolve().parent)
    main(parser.parse_args().carpeta)


