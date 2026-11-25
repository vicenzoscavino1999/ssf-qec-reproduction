"""HGP (3,4). Posición del cuarto dato y cierre por expansión local. Biblioteca estándar."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
from runpy import run_path


HERE = Path(__file__).resolve().parent
engine = run_path(str(HERE / "verificar_H5A_3_4_cola_por_lectores.py"))
cross = run_path(str(HERE / "verificar_H5A_3_4_reuniones_cruzadas.py"))
blocks = engine["v19"]
M = 32
FOUR = ((1, 3), (2, 2))
FOUR_AND_THREE = (*FOUR, (0, 3), (1, 2))


def closure_bound(v, c):
    """Coordenadas y centros incidentes por eje con la condición de expansión E32."""
    n = v + c
    coords, readers = [], []
    for stars, singles, degree, opposite in ((c, v, 4, 3), (v, c, 3, 4)):
        assert stars+n+1 <= M
        coords.append(degree*stars+singles+(degree-2)*n)
        # Estrellas repetidas reducen el extremo de este mayorante.
        readers.append(opposite*singles+(opposite*degree-degree+1)*stars
                       +(opposite*(degree-2)-degree+1)*n)
    return tuple(coords), tuple(readers)


def multiple_limits():
    # Un eje de radio cuatro contiene a lo sumo un 4-ciclo que testifique multiplicidad.
    vv, cc = 9*((3-2)*(4-1)), 8*((4-2)*(3-1))
    # En ejes de longitud tres, intersección de dos subconjuntos de destinos en un producto.
    vc = 18**2-(18-1)*(18-1)
    cv = 24**2-(24-2)*(24-1)
    return ((vv, vc), (cv, cc))


def reader_geometry():
    """Máscaras del lector y vecindades cerradas de cada dato en el recubrimiento."""
    tree = engine["v21"]["SeedTree"]
    row, col = tree(1, 0), tree(0, 0)
    rows, cols = row.adjacent(0), col.adjacent(0)
    cells = tuple(product(rows, cols))
    sources = (tuple(product(rows, sorted(col.ball(0, 2)))),
               tuple(product(sorted(row.ball(0, 2)), cols)))
    assert tuple(map(len, sources)) == (40, 27)
    masks, reached = {}, {}
    for kind, sites in enumerate(sources):
        for x, y in sites:
            stripe = (set(product((x,), col.adjacent(y))) if kind == 0
                      else set(product(row.adjacent(x), (y,))))
            masks[kind, (x, y)] = sum(1 << j for j, cell in enumerate(cells) if cell in stripe)
    assert all(masks.values())
    for dst in (0, 1):
        index = {}
        for kind, sites in enumerate(sources):
            rx, ry = engine["v21"]["DISTANCES"][kind][dst]
            for x, y in sites:
                neighborhood = product(row.ball(x, rx), col.ball(y, ry))
                bits = 0
                for pos in neighborhood:
                    bits |= 1 << index.setdefault(pos, len(index))
                reached[dst, kind, (x, y)] = bits
    return sources, masks, reached


def triple_volume():
    """Cota por raíz para toda terna donante y su posible cuarto dato."""
    sources, masks, reached = reader_geometry()
    by_root = (Counter(), Counter())
    for v in range(4):
        c = 3 - v
        for xv, xc in product(combinations(sources[0], v), combinations(sources[1], c)):
            donors = tuple((0, x) for x in xv) + tuple((1, x) for x in xc)
            a, b, d = (masks[x] for x in donors)
            if min((a|b).bit_count(), (b|d).bit_count(), (a|d).bit_count()) < 2:
                continue
            if (a|b|d).bit_count() < 3:
                continue
            for i, count in enumerate((v, c)):
                by_root[i][v, c] += count
            for dst in (0, 1):
                union = 0
                for kind, site in donors:
                    union |= reached[dst, kind, site]
                outside = union.bit_count() - (v, c)[dst]
                assert outside >= 0
                key = v + (dst == 0), c + (dst == 1)
                for i, count in enumerate(key):
                    by_root[i][key] += count * outside
    for i, norm in enumerate((Q(3, 4), Q(4, 3))):
        for key in by_root[i]:
            by_root[i][key] *= norm
    return by_root


@lru_cache(None)
def odd_owners(v, c):
    """Vectores de celdas impares asignadas en un bloque de solapamiento conexo."""
    widths = (3,)*v + (4,)*c
    n = len(widths)
    assert 1 <= n <= 4
    if n == 1:
        return (widths,)
    labels = tuple(range(n))
    pairs = tuple(combinations(labels, 2))
    groups = [g for size in range(n, 1, -1) for g in combinations(labels, size)]
    inside = {g: tuple(j for j, (a, b) in enumerate(pairs) if a in g and b in g) for g in groups}
    pair_limits = tuple(2 if widths[a] == widths[b] else 1 for a, b in pairs)
    same_kind = [set(g) for g in groups if len({widths[i] for i in g}) == 1]
    out = set()

    def enumerate_cells(k, left, used_pairs, shared):
        if k == len(groups):
            reached = {0}
            for _ in labels:
                for group in shared:
                    if reached.intersection(group):
                        reached.update(group)
            if len(reached) < n:
                return
            for same in same_kind:
                collisions = sum(max(0, len(same.intersection(g))-1) for g in shared)
                if collisions > len(same):
                    return
            owned = list(left)
            for g in shared:
                if len(g) % 2:
                    owned[g[0]] += 1
            out.add(tuple(owned))
            return
        g, indices = groups[k], inside[groups[k]]
        limit = min([left[j] for j in g]+[pair_limits[j]-used_pairs[j] for j in indices])
        for count in range(limit+1):
            new_left = tuple(x-count if i in g else x for i, x in enumerate(left))
            new_pairs = tuple(x+count if j in indices else x for j, x in enumerate(used_pairs))
            enumerate_cells(k+1, new_left, new_pairs, shared+(g,)*count)

    enumerate_cells(0, widths, (0,)*len(pairs), ())
    return tuple(sorted(out))


@lru_cache(None)
def expanded_block(v, c):
    original = blocks["block_states"](v, c)
    widths = (3,)*v + (4,)*c
    out = {}
    for owned in odd_owners(v, c):
        ranks = [range(max(0, u-w//2), u+1) for u, w in zip(owned, widths)]
        for r in product(*ranks):
            key = blocks["clean_tuple"](r[:v]), blocks["clean_tuple"](r[v:])
            if key in original:
                balance = min(original[key], 2*sum(r)-sum(owned))
                out[key] = max(out.get(key, balance), balance)
    return out


def expanded_states(shape):
    return blocks["combined_states"](shape, expanded_block)


def improved_tables(data, circuit, volume, expansion):
    hist, universal, real = data[3], circuit[6], circuit[7]
    states = expanded_states if expansion else blocks["combined_states"]
    caps, circuit_caps = blocks["block_caps"](hist, volume, 4, universal, states=states)
    if expansion:
        cuts = blocks["separated_caps"](hist, universal, 4, states=states)
        for dst, extra in zip((caps, circuit_caps), cuts):
            for key, cap in extra.items():
                dst[key] = min(dst.get(key, cap), cap)
    apply = engine["v18"]["apply_caps"]
    head = apply(hist, caps)
    all_real_caps = {}
    for i in (0, 1):
        for key, row in head[i].items():
            if 3 <= sum(key) <= 4:
                for k in range(1, max(real[i][key])+1):
                    all_real_caps[(i, *key), k] = sum(weight for r, weight in row.items() if r >= k)
    new_universal = apply(universal, circuit_caps)
    new_real = apply(apply(real, circuit_caps), all_real_caps)
    return (*data[:3], head), (*circuit[:6], new_universal, new_real)


CASES = (
    ("General, Bernoulli", "5.30e-7", None, ".567", ".685",
     (("9.997404931638269e-7", "1.16609506010893e-6"),
      ("7.014209285135773e-7", "7.534346312947361e-7"),
      ("3.882446482721063e-8", "1.416450718206368e-7"),
      ("4.920597729886739e-9", "2.06751689910055e-8"),
      "5.531768194698058e-7", "1.888311839184505e-6")),
    ("General, VV primero", "9.92e-8", 4, ".558", ".750",
     (("1.418817501511518e-6", "1.625864155460227e-6"),
      ("4.578329336245391e-7", "8.586713081195475e-7"),
      ("5.110492363940325e-8", "2.659048003189979e-7"),
      ("6.68674451626325e-9", "3.717250911461154e-8"),
      "7.050632238707232e-7", "1.441412884898932e-6")),
    ("General, cualquier orden admitido", "9.05e-8", 7, ".555", ".780",
     (("1.295293254418869e-6", "1.484553790989989e-6"),
      ("6.580918680974423e-7", "7.896534898621902e-7"),
      ("8.282960346522432e-8", "2.892054846438966e-7"),
      ("8.863053201997362e-9", "3.960469677880701e-8"),
      "8.036099814785975e-7", "1.296551643239366e-6")),
    ("E32, Bernoulli", "5.30e-7", None, ".567", ".685",
     (("9.984440912668088e-7", "1.164152150884744e-6"),
      ("7.005606751139353e-7", "7.506369183432166e-7"),
      ("3.744368704963447e-8", "1.344667630634721e-7"),
      ("4.524604328556123e-9", "1.845680835087707e-8"),
      "4.9896843790715e-7", "1.718873084737384e-6")),
    ("E32, VV primero", "9.95e-8", 4, ".558", ".750",
     (("1.422715199713484e-6", "1.630195002210134e-6"),
      ("4.590284636666641e-7", "8.60070328221362e-7"),
      ("5.085891010099629e-8", "2.633837077372463e-7"),
      ("6.584543179309267e-9", "3.619158569867787e-8"),
      "6.903173952614415e-7", "1.408770096814824e-6")),
    ("E32, cualquier orden admitido", "9.08e-8", 7, ".555", ".780",
     (("1.299179587785603e-6", "1.48893298988157e-6"),
      ("6.599490347674885e-7", "7.911293801964269e-7"),
      ("8.225833174221509e-8", "2.859324484708223e-7"),
      ("8.662928046260095e-9", "3.847426921505588e-8"),
      "7.894062778799531e-7", "1.273534243010105e-6")),
)


def check_local():
    assert closure_bound(1, 3) == ((21, 10), (42, 30))
    assert closure_bound(2, 2) == ((18, 12), (36, 36))
    assert multiple_limits() == ((27, 35), (70, 32))
    assert expanded_block(1, 0) == blocks["block_states"](1, 0)
    assert expanded_block(0, 1) == blocks["block_states"](0, 1)
    volume = triple_volume()
    for i, expected in enumerate(((15617916, 66985812), (83295552, 119085888))):
        assert tuple(volume[i][k] for k in FOUR) == expected
    print("Posiciones de la terna y el cuarto dato, raíz CC:", *(volume[1][k] for k in FOUR), flush=True)
    return volume


def check_global(volume):
    engine["check_local"]()
    data, circuit = engine["v21"]["previous_tables"]()
    assert data[2] == circuit[3] and engine["F5"]["RM"] == Q(".001")
    discount, pair = cross["coverage_discounts"](engine, circuit[4], data[2])
    data = (*data[:2], discount, data[3])
    circuit = (*circuit[:3], discount, *circuit[4:])
    old = blocks["triple_reader"](circuit[1], 4)
    for i in (0, 1):
        assert all(0 <= count <= old[i][k] for k, count in volume[i].items())
    assert [old[1][k] for k in FOUR] == [109867680, 160982640]
    assert [data[3][1][k][5] for k in FOUR] == [154769832, 175231216]
    for conditional, expected in enumerate(((128197704, 133334464), (127333704, 132590272))):
        improved, improved_circuit = improved_tables(data, circuit, volume, conditional)
        assert tuple(improved[3][1][k][5] for k in FOUR) == expected
        limits = multiple_limits() if conditional else None
        for case, prior in zip(CASES[conditional*3:conditional*3+3], engine["CASES"]):
            engine["verify_model"](case, improved, improved_circuit,
                                   floor_fields=cross["overlap_floor"](engine, prior),
                                   disjoint_pair_factors=pair, multiple_limits=limits)
        floor = cross["overlap_floor"](engine, engine["CASES"][0])
        for a, s in ((".567", ".685"), (".569", ".673")):
            def partial_obstruction(compositions):
                return engine["obstruction"](improved, Q(a), Q(s), floor_fields=floor,
                                             multiple_limits=limits, compositions=compositions)
            assert Q(".81") < partial_obstruction(FOUR) < Q(".84")
            assert partial_obstruction(FOUR_AND_THREE) > Q("1.074")
    print("Cotas y obstrucciones comprobadas con aritmética racional.", flush=True)


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para comprobar todas las desigualdades.")
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--global", dest="global_row", action="store_true",
                        help="reconstruir las tablas y verificar las filas completas")
    args = parser.parse_args()
    volume = check_local()
    if args.global_row:
        check_global(volume)

