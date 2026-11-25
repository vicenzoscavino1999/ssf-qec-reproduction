"""Encuentros de grupos y consultas de medida. Certificado H5-A, Python 3.10+."""

from fractions import Fraction as F
from heapq import heapify, heappop, heappush
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verificar_H5A_primer_apagado_por_grupos as last
import verificar_H5A_propagacion_entre_fuentes as spatial


def data_from_check(radius):
    vv, cc, vc, cv = spatial.ball_sizes(radius)
    return vv * cv + cc * vc


def rows():
    p, q = F('1e-30'), F('2e-11')
    odds, exp = last.prior.odds, last.prior.exp
    sp, sq = (odds(r) * exp(F('.015')) for r in (p, q))
    light_sites = 18 + 100 * 18
    _, _, measure_neighbours, _, rect = spatial.previous.geometry()
    assert (measure_neighbours, rect) == (1420146, 83810)
    vv, cc, vc, cv = spatial.ball_sizes(5)
    vv_next, cc_next, _, _ = spatial.ball_sizes(6)
    same = max(vv * vv_next + vc**2, cv**2 + cc * cc_next)
    contact = light_sites * rect * data_from_check(4)
    outside = same + light_sites * rect * data_from_check(3)
    reads = 3 * light_sites * rect + 18
    assert reads == 457099758 and reads * sq < F('.00929')

    data, passive = F('1.4e-30'), F('2.031e-11')
    bd, bm = F('.00587'), F('.000118')
    assert outside * data + 135 * measure_neighbours * passive < bd
    assert data_from_check(5) * data + rect * passive < bm
    assert sq * exp(2 * bm) < passive
    # De g_{r+1}=272g_r+1 se deduce, para las cuatro bolas y r>=1,
    # b_{r+1}=272 b_r+offset. Los offsets son <= b_1.
    one = spatial.ball_sizes(1)
    offsets = tuple(b - 272*a for a, b in zip(one, spatial.ball_sizes(2), strict=True))
    assert offsets == (18, 17, 17, 18)
    assert all(0 <= off <= b for off, b in zip(offsets, one, strict=True))
    growth = 273**2
    r_group = F(growth, 32**5)
    assert 0 < r_group < 1
    extra = r_group / (32 * (1 - r_group))
    # Dos inclinaciones: una por el núcleo de la medida y otra por una unión
    # de grupos que se simula junta. Su producto nunca excede 2.
    tilted = F('5.8e-30')
    theta = F(7, 5)
    r_union = growth / theta**40
    assert theta**2 < 2 and 0 < r_union < 1
    union_extra = r_union / (theta * (1 - r_union))
    union_contacts = light_sites * rect * data_from_check(3) * tilted * (1 + union_extra)
    assert union_contacts < F('.000016')
    for factor, root, group in ((1, data, F('5e-29')), (2, tilted, F('2.8e-28'))):
        overhead = contact * (root + extra * group) + union_contacts + 2 * bd + reads * sq
        assert factor * sp * exp(same * root + overhead) < root
        assert 32 * factor * sp * exp(same * group + overhead) < group
    active = sq * data_from_check(6) * tilted / (theta * (1 - r_union))
    assert passive + active < F('8.6e-10')
    print('Consultas frescas por dato <=', reads, 'en la reconstrucción del núcleo.')
    print('Fila espacial <= 1.4e-30 por dato y < 8.6e-10 por medida.')
    return data, passive + active


def long_reach():
    # Si una luz supera el radio admitido R, los nueve checks consumidos
    # en ese paso implican (101/99)^R <= 2n en una ejecución sin medidas.
    # La inducción avanza de cinco en cinco para incluir el suelo entero.
    r = F(101, 99)
    for base, step, start in ((4, 5, 2100), (3, 40, 22000)):
        assert r > F(start + step, start)
        for n in range(start, start + step):
            assert r**(base + (n-1)//step) > 2*n
    print('Sin medidas, los radios cubren el recorrido desde 2100 y 22000 datos, '
          'respectivamente antes y después de agrupar.')


def layered_graph(depth):
    """Grafo 17-regular. Cada arista corresponde a un check distinto."""
    sizes = [2 * 8**j * 9**(depth-j) for j in range(depth+1)]
    sizes = sizes + [sizes[-1]] + sizes[::-1]
    layers = []
    first = 0
    for size in sizes:
        layers.append(range(first, first + size))
        first += size
    neighbours = [[] for _ in range(first)]

    def edge(v, w):
        assert v != w
        neighbours[v].append(w)
        neighbours[w].append(v)

    for layer in (layers[0], layers[-1]):
        size = len(layer)
        for v in layer:
            for shift in range(1, 5):
                edge(v, layer[(v-layer.start+shift) % size])
        for j in range(size//2):
            edge(layer[j], layer[j+size//2])
    centre = layers[depth+1]
    for j in range(len(centre)//2):
        edge(centre[j], centre[j+len(centre)//2])
    for arm in (layers[:depth+2], layers[::-1][:depth+2]):
        for j, (source, dest) in enumerate(zip(arm, arm[1:])):
            incoming = 8 if j == depth else 9
            assert len(source)*8 == len(dest)*incoming
            for k in range(len(source)*8):
                edge(source[k//8], dest[k % len(dest)])
    assert all(len(row) == 17 for row in neighbours)
    return neighbours, layers


def run(neighbours, seeds):
    occupied = [False] * len(neighbours)
    for v in seeds:
        occupied[v] = True
    hits = [sum(occupied[w] != occupied[v] for w in row)
            for v, row in enumerate(neighbours)]
    pending = [(-h, v) for v, h in enumerate(hits) if h >= 9]
    heapify(pending)
    trace = []
    while pending:
        neg_h, v = heappop(pending)
        if -neg_h != hits[v]:
            continue
        assert not occupied[v]
        trace.append((v, hits[v]))
        occupied[v] = True
        hits[v] = 17 - hits[v]
        for w in neighbours[v]:
            hits[w] += -1 if occupied[w] else 1
            if hits[w] >= 9:
                heappush(pending, (-hits[w], w))
    assert max(hits) < 9
    return trace, sum(hits)//2


def contrast():
    depth = 3
    graph, layers = layered_graph(depth)
    left, right, centre = layers[0], layers[-1], layers[depth+1]
    solo = [run(graph, group) for group in (left, right)]
    joint, terminal = run(graph, [*left, *right])
    assert terminal == 0 and all(boundary == len(centre)*8 for _, boundary in solo)
    indices = [0, 0]
    reached = [set(left), set(right)]
    for step, (v, h) in enumerate(joint, 1):
        options = []
        for k, (trace, _) in enumerate(solo):
            if indices[k] < len(trace):
                w, lights = trace[indices[k]]
                options.append((lights, -w, k))
        winner = max(options, default=None)
        if winner is None or -winner[1] != v:
            assert v in centre
            contributors = [sum(w in side for w in graph[v]) for side in reached]
            assert max(contributors) <= 8 and sum(contributors) == h >= 9
            # La fila central no compite en ninguna copia aislada.
            assert all(all(w not in centre for w, _ in trace) for trace, _ in solo)
            print(f'Contraste físico m={depth}: primera diferencia en el paso {step}, '
                  f'con {contributors[0]}+{contributors[1]} luces. '
                  f'Síndromes finales {solo[0][1]}, {solo[1][1]} y {terminal}.')
            return
        assert winner[0] == h
        side = winner[2]
        reached[side].add(v)
        indices[side] += 1
    raise AssertionError('Falta el primer cambio de corrección.')


def third_group():
    """Cierre monótono por capas de igual vecindad. Las entradas son grados."""
    depth, trunk = 3, 3
    degrees = {}

    def link(first, second, forward, backward):
        degrees.setdefault(first, {})[second] = forward
        degrees.setdefault(second, {})[first] = backward

    arms = [[(side, j) for j in range(depth+1)] for side in 'ABC']
    for arm in arms:
        link(arm[0], arm[0], 9, 9)
        for before, after in zip(arm, arm[1:]):
            link(before, after, 8, 9)
    one, two = 'unión AB', 'unión con C'
    link(arms[0][-1], one, 8, 5)
    link(arms[1][-1], one, 8, 4)
    path = [one, *(('tramo', j) for j in range(trunk))]
    for before, after in zip(path, path[1:]):
        link(before, after, 8, 9)
    link(path[-1], two, 8, 8)
    link(arms[2][-1], two, 8, 8)
    link(two, two, 1, 1)
    assert all(sum(row.values()) == 17 for row in degrees.values())
    # Proporciones de las poblaciones. Un múltiplo par común realiza los
    # grados intercapas y permite completar los grafos internos.
    size = {one: F(1)}
    todo = [one]
    while todo:
        v = todo.pop()
        for w, deg in degrees[v].items():
            required = size[v] * deg / degrees[w][v]
            if w in size:
                assert size[w] == required
            else:
                size[w] = required
                todo.append(w)
    assert len(size) == len(degrees) and min(size.values()) > 0

    def closure(seeds):
        seen = set(seeds)
        while True:
            new = {v for v, row in degrees.items()
                   if v not in seen and sum(d for w, d in row.items() if w in seen) >= 9}
            if not new:
                return seen
            seen.update(new)

    def ends_of_reached_checks(active):
        return active | {w for v in active for w in degrees[v]}

    a, b, c = [arm[0] for arm in arms]
    alone = [closure([v]) for v in (a, b, c)]
    joined = closure([a, b])
    assert one not in alone[0] | alone[1] and two not in alone[2]
    assert all(ends_of_reached_checks(local).isdisjoint(ends_of_reached_checks(alone[2]))
               for local in alone[:2])
    assert path[-1] in joined and two not in joined
    assert {one, two}.isdisjoint(closure([a, c]) | closure([b, c]))
    marks = {}
    for mask in range(8):
        occupied = closure([v for j, v in enumerate((a, b, c)) if mask & (1<<j)])
        boundary = any(w not in occupied for v in occupied for w in degrees[v])
        marks[mask] = int(not boundary)
    assert [marks[m] for m in range(8)] == [1, 0, 0, 0, 0, 0, 0, 1]
    assert closure([a, b, c]) == set(degrees)
    print('Tres grupos: sólo la corrida conjunta de A y B llega al contacto con C.')


def verify():
    sector = rows()
    long_reach()
    contrast()
    third_group()
    # Cotas de las filas anteriores. Recalcular las dos mayores, cuyo redondeo
    # para publicación no conserva el margen que usamos al incorporar el sector.
    previous = spatial.previous
    stored_upper = (
        (F('3.68e-11'), F('2.04e-11')),  # propagación, árboles direccionales
        (F('1.1e-29'), F('2.04e-11')),   # geometría anterior
        (F('2e-16'), F('2.057e-11')),    # ocurrencias con enlace directo
        (F('3.53e-11'), F('2.05e-11')),  # primeros usos
    )
    terms = (previous.two_data(), previous.activity_counts(), *stored_upper, sector)
    ref = F('1.907e-6') * last.prior.odds(F('1e-30')) / last.prior.odds(F('1e-6'))
    totals = [sum(row) + ref for row in zip(*terms, strict=True)]
    assert totals[0] < F('1.94e-5') and totals[1] < F('.00070')
    assert max(totals) + F('.004298') < F('.004998')
    print('Mayorantes acumulados en la meta:', *(f'{float(t):.10g}' for t in totals))
    print('PASS parcial en p=1e-30, q=2e-11. Resto conjunto abierto.')


if __name__ == '__main__':
    verify()
