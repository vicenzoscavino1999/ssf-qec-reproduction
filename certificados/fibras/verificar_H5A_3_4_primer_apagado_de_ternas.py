"""Primeros apagados de 3C y 1V+2C que permiten conectar con otros grupos."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations, product
from math import prod
from pathlib import Path
from runpy import run_path


HERE = Path(__file__).resolve().parent
last = run_path(str(HERE / "verificar_H5A_3_4_cierre_por_expansion.py"))
engine = last["engine"]
v18 = engine["v18"]
iso = v18["previous"]("verificar_H5A_3_4_medidas_aisladas.py")
contact = v18["previous"]("verificar_H5A_3_4_contacto_de_checks.py")
MIXED, THREE_C = (0, 1, 1), (1, 1, 1)
PAIRS = ((0, 1), (0, 2), (1, 2))


@lru_cache(None)
def reader_profiles(width, other_neighbors):
    """Máscaras de los centros que tocan una vecindad, hasta permutar sus sitios."""
    full = (1 << width)-1
    relabel = [tuple(sum(1 << p[i] for i in range(width) if m >> i & 1)
                     for m in range(full+1)) for p in permutations(range(width))]
    out = set()

    def enumerate_masks(mask, left, seen):
        if not mask:
            if not any(left):
                out.add(min(tuple(sorted(rename[x] for x in (full, *seen))) for rename in relabel))
            return
        rows = [i for i in range(width) if mask >> i & 1]
        for count in range(1+min(left[i] for i in rows)):
            new_left = tuple(n-count if i in rows else n for i, n in enumerate(left))
            enumerate_masks(mask-1, new_left, seen+(mask,)*count)

    enumerate_masks(full, (other_neighbors,)*width, ())
    return tuple(sorted(out))


def mixed_owners(edges):
    """Celdas impares de cada franja para 1V+2C. La primera etiqueta posee la triple."""
    vc1, vc2, cc = edges
    hit1, hit2 = int(vc1 == 4), int(vc2 == 4)
    for overlap in (range(1, 5) if cc == 4 else (0,)):
        shared = int(cc in (3, 4) and hit1 and hit2)
        if (shared and not overlap) or (overlap == 4 and hit1 != hit2):
            continue
        yield 3-hit1-hit2+2*shared, 4-overlap-hit1+shared, 4-overlap-hit2+shared


def first_moves():
    profiles, stable, reps = contact["geometry"]()
    perfect = sorted({phi for lit, phi in engine["FAR"]["winners"](reps) if lit == phi})
    assert len(perfect) == 63
    for phi in perfect:
        supports = [sum(1 << i for i, row in enumerate(iso["ROWS"]) if row & col & phi)
                    for col in iso["COLS"]]
        if all(mask.bit_count() == 2 for mask in supports):
            assert all(mask in (supports[0], 15^supports[0]) for mask in supports)
    mixed, overlapping = {}, set()
    for types, threshold in ((MIXED, 3), (THREE_C, 4)):
        options = [tuple((0, *iso["positions"](types[i], types[j]))) for i, j in PAIRS]
        canonical = sorted({v18["canonical"](types, e)[1] for e in product(*options)
                            if sum(bool(x) for x in e) >= 2})
        for edges in canonical:
            if types == THREE_C and sum(x >= 3 for x in edges) == 2 and 0 not in edges:
                continue
            best, debit, perfect_moves = iso["triple_budget"](types, edges, profiles, stable,
                                                              with_debits=True)
            if max(edges) == 1:
                assert best <= 1
            if best < threshold:
                continue
            widths = tuple(3+t for t in types)
            if 4 in edges:
                if types == THREE_C:
                    # Cuando no hay lector común todos los pasos tienen ≤2 donantes.
                    # La tabla v18 elimina estos patrones incluso al permitir tres.
                    assert 0 in edges and v18["limited_budget"](types, edges, profiles, stable,
                                                               (iso["ROWS"], iso["COLS"])) < threshold
                    continue
                for owned in mixed_owners(edges):
                    ranks = product(*(range(max(0, u-w//2), u+1) for u, w in zip(owned, widths)))
                    if any(2*sum(r)-sum(owned)-debit[r] >= threshold for r in ranks):
                        overlapping.add(edges)
                        break
                continue
            assert max(2*sum(r)-sum(widths)-cost for r, cost in debit.items()
                       if r != widths) < threshold
            choices = {}
            for step, first_cost in perfect_moves.items():
                remaining = tuple(a-b for a, b in zip(widths, step))
                largest = sum(widths)-threshold-debit[remaining]
                if first_cost <= largest:
                    choices[step] = largest
            if not choices:
                continue
            assert min(min(a) for a in choices) > 0 and 0 not in edges
            if types == THREE_C:
                assert edges in ((1, 1, 2), (1, 2, 2), (2, 2, 2))
                assert all((a in ((1, 1, 1), (2, 2, 2)) and upper == sum(a))
                           for a, upper in choices.items())
            else:
                mixed[edges] = choices
                # Etiquetas de los dos CC intercambiadas.
                rev = edges[1], edges[0], edges[2]
                swapped = {(a, c, b): bound for (a, b, c), bound in choices.items()}
                assert edges != rev or choices == swapped
                mixed[rev] = swapped
    return mixed, overlapping, perfect, stable


def overlapping_volume(overlapping, root_kind):
    """Árboles de posiciones de v18, sólo en los patrones de solapamiento posibles."""
    types = MIXED if root_kind == 0 else (1, 0, 1)
    options = [(0, *iso["positions"](types[i], types[j])) for i, j in PAIRS]
    trees = tuple(contact["rooted_trees"](3))
    total = 0
    for edges in product(*options):
        if v18["canonical"](types, edges)[1] not in overlapping:
            continue
        rel = dict(zip(PAIRS, edges))
        costs = []
        for tree in trees:
            if all(rel[tuple(sorted(e))] for e in tree):
                costs.append(prod(iso["positions"](types[i], types[j])[rel[tuple(sorted((i, j)))]]
                                  for i, j in tree))
        total += min(costs)
    return Q(total, 2 if root_kind == 0 else 1)


def disjoint_sites(choices, perfect, stable):
    row_profiles, col_profiles = reader_profiles(4, 2), reader_profiles(3, 3)
    assert (len(row_profiles), len(col_profiles)) == (23, 16)
    # Los perfiles incluyen el centro y permiten centros distintos de igual máscara.
    rows, cols = iso["ROWS"], iso["COLS"]

    @lru_cache(None)
    def mixed_event(vv, c1, c2):
        if (vv & c1) or (vv & c2) or (c1 & c2):
            return False
        row = next(mask for mask in rows if vv & mask)
        col1, col2 = (next(mask for mask in cols if x & mask) for x in (c1, c2))
        e1 = 2 if (vv & col1 or c1 & row) else 1
        e2 = 2 if (vv & col2 or c2 & row) else 1
        # Un vecino común de filas de los CC puede quedar fuera de este lector.
        same_row = any((r & c1) and (r & c2) for r in rows)
        cc = (3,) if col1 == col2 else ((2,) if same_row else (1, 2))
        allowed = [choices[e1, e2, e] for e in cc if (e1, e2, e) in choices]
        if not allowed:
            return False
        for phi in perfect:
            used = tuple((phi & x).bit_count() for x in (vv, c1, c2))
            if min(used) and stable[phi & ~(vv|c1|c2)]:
                if any(phi.bit_count() <= table.get(used, 0) for table in allowed):
                    return True
        return False

    @lru_cache(None)
    def double_event(a, b, c):
        # Dos filas donadas en cada columna: las tres parejas son iguales o complementarias.
        # Un vecino común en las tres filas ya está cubierto por el recuento directo.
        if a & b & c:
            return False
        for split in (3, 5, 6):
            complement = 15 ^ split
            if all((x & split) == split or (x & complement) == complement for x in (a, b, c)):
                return True
        return False

    def by_anchor(peak, profiles, transverse_degree):
        return transverse_degree * max(sum(peak[mask.bit_count()] for mask in p) for p in profiles)

    multi = Counter()
    mixed_peak = (Counter(), Counter())
    for row_profile in row_profiles:
        for root in row_profile:
            # Dos columnas restantes distintas y ordenadas por su posición en el lector.
            count = sum(double_event(root, u, w) for u, w in product(row_profile, repeat=2))
            multi[root.bit_count()] = max(multi[root.bit_count()], count)
        c_masks = tuple(sum(1 << (3*i+j) for i in range(4) if mask >> i & 1)
                        for mask in row_profile for j in range(3))
        for col_profile in col_profiles:
            v_masks = tuple(mask << (3*i) for i in range(4) for mask in col_profile)
            v_counts, c_counts = [0]*len(v_masks), [0]*len(c_masks)
            for j, k in combinations(range(len(c_masks)), 2):
                x, y = sorted((c_masks[j], c_masks[k]))
                if x & y:
                    continue
                for i, vv in enumerate(v_masks):
                    if mixed_event(vv, x, y):
                        v_counts[i] += 1
                        c_counts[j] += 1
                        c_counts[k] += 1
            for peak, counts, masks in zip(mixed_peak, (v_counts, c_counts), (v_masks, c_masks)):
                for n, mask in zip(counts, masks):
                    peak[mask.bit_count()] = max(peak[mask.bit_count()], n)
    triple_c = 4*4*3*3 + by_anchor(multi, row_profiles, 4)
    mixed_volume = by_anchor(mixed_peak[0], col_profiles, 3), by_anchor(mixed_peak[1], row_profiles, 4)
    return triple_c, mixed_volume


def geometry_caps():
    choices, overlapping, perfect, stable = first_moves()
    ccc, disjoint = disjoint_sites(choices, perfect, stable)
    overlap = tuple(overlapping_volume(overlapping, i) for i in (0, 1))
    total = tuple(a+b for a, b in zip(overlap, disjoint))
    print("3C a saldo cuatro:", ccc, "  1V+2C a saldo tres: disjuntos", disjoint,
          "solapados", tuple(map(int, overlap)), "total", tuple(map(int, total)), flush=True)
    return {((1, 0, 3), 4): ccc, ((0, 1, 2), 3): total[0], ((1, 1, 2), 3): total[1]}


def tighten(data, circuit, caps):
    # Tres marcas con alguna paridad cero dejan a lo sumo dos datos reales y saldo ≤2.
    apply = engine["v18"]["apply_caps"]
    hist, universal, real = (apply(old, caps) for old in (data[3], circuit[6], circuit[7]))
    return (*data[:3], hist), (*circuit[:6], universal, real)



def monomial_coefficients(values):
    """Interpolación exacta en 0,1,... en la base de Newton."""
    coefficients = [Q(0)]*len(values)
    basis = [Q(1)]
    for j in range(len(values)):
        for t, coefficient in enumerate(basis):
            coefficients[t] += coefficient*values[0]
        values = [b-a for a, b in zip(values, values[1:])]
        basis = [(a-j*b)/(j+1) for a, b in zip([Q(0), *basis], [*basis, Q(0)])]
    assert min(coefficients) >= 0
    return coefficients


@lru_cache(None)
def counts_by_strong_edges(prices):
    blocks = engine["v19"]
    def scaled(table, m):
        return tuple(tuple(Q(m*x) for x in row) for row in table)
    weak = scaled(blocks["WEAK"], 1)
    values = [engine["polynomials"](prices, (weak, scaled(blocks["STRONG_DISJOINT"], x),
                                            scaled(blocks["OVERLAP"], x)), 10)
              for x in range(10)]
    out = [{}, {}]
    for i in (0, 1):
        for key in values[0][i]:
            degree = sum(key)
            if degree < 3:
                continue
            coefficients = monomial_coefficients([table[i][key] for table in values[:degree]])
            assert all(sum(c*x**j for j, c in enumerate(coefficients)) == table[i][key]
                       for x, table in enumerate(values))
            running = Q(0)
            cumulative = []
            for c in coefficients:
                running += c
                cumulative.append(running)
            out[i][key] = cumulative
    return out


@lru_cache(None)
def edge_gains(a, s, fields):
    v21 = engine["v21"]
    regular, marked, cc_volume = {}, {}, {}
    for ij in ((0, 0), (0, 1), (1, 1)):
        def features(x, y):
            geom = v21["initial_overlap"](*ij, x, y)
            return (*geom, *(v21["link_overlap"](*ij, x, y, t, a**6, s) for t in range(3)))
        weak = features(*v21["DISTANCES"][ij[0]][ij[1]])
        gains = [tuple(value-reference for value, reference in zip(features(x, y), weak))
                 for x, y in v21["displacements"](*ij) if engine["edge_kind"](*ij, x, y)]
        assert all(min(row) >= 0 for row in gains)
        def cost(row, k):
            return sum(x*y for x, y in zip(row[:k], fields[:k]))
        regular[ij] = min(cost(row, 8) for row in gains)
        marked[ij] = min(cost(row, 5) for row in gains)
        cc_volume[ij] = min(row[6] for row in gains)
    return regular, marked, cc_volume


def adjust_heads(regular, marked, discounts, data, circuit, rate, first, a, s, fields):
    """Cruzar el número de aristas fuertes con el coste de los sucesores exteriores."""
    geometrical, physical, _ = edge_gains(a, s, tuple(fields))
    prices = engine["rounded_prices"](engine["source_prices"](a, rate if first is not None else None, first))
    cum = counts_by_strong_edges(prices)
    hist = data[3] if first is None else circuit[6]
    y = a**6
    for i in (0, 1):
        for key, by_edge in cum[i].items():
            v, c = key
            n = v+c
            types = [p for p in geometrical if (p != (0, 0) or v >= 2)
                     and (p != (1, 1) or c >= 2) and (p != (0, 1) or v*c)]
            factor = engine["v21"]["upper_factor"](min(geometrical[p] for p in types))
            raw = regular[i][key]
            discount = discounts.get(key, Q(1))
            bound = factor**(n-1)*raw
            for j, sub in enumerate(by_edge[:-1]):
                limited = min(raw, discount*sub)
                target_triple = key in ((0, 3), (1, 2))
                if first is None or (j == 0 and target_triple):
                    if j == 0:
                        cap = 1 if target_triple else 4*n//6
                    else:
                        cap = (7*v+8*c+2*min(c, 2*j))//6
                    limited = min(limited, discount*sum(count/y**min(k, cap)
                                                        for k, count in hist[i][key].items()))
                bound += (1-factor)*factor**j*limited
            regular[i][key] = min(raw, bound)
            if key in ((0, 3), (1, 2)):
                assert engine["ANCHOR"](engine["BALL"], i, 1) == 0
                marked[i][key] *= engine["v21"]["upper_factor"](min(physical[p] for p in types))


CASES = (
    ("Bernoulli", "5.40e-7", None, ".568", ".678",
     (("1.018863461563783e-6", "1.18847576096107e-6"),
      ("7.146407184110793e-7", "7.663140985710223e-7"),
      ("3.82659701597073e-8", "1.314458104221317e-7"),
      ("4.894474882845396e-9", "1.960304117424008e-8"),
      "5.41912734408362e-7", "1.949873246352219e-6")),
    ("Circuito, VV primero", "1.005e-7", 4, ".559", ".745",
     (("1.436956340250215e-6", "1.646544704287874e-6"),
      ("4.634411415673427e-7", "8.664768413714867e-7"),
      ("4.966776923092021e-8", "2.455532712949687e-7"),
      ("6.454357188433469e-9", "3.392340581324243e-8"),
      "6.89169849509594e-7", "1.437963471916807e-6")),
    ("Circuito, cualquier orden admitido", "9.18e-8", 7, ".556", ".771",
     (("1.315344624823385e-6", "1.507931954091591e-6"),
      ("6.682455477868437e-7", "8.021741585389381e-7"),
      ("8.370643723694925e-8", "2.831153329273515e-7"),
      ("9.351161378705897e-9", "4.118914525211162e-8"),
      "8.567439531484622e-7", "1.438806511965248e-6")),
)


def reference_fields(case):
    _, u, h, physical, visible, root = case[-1]
    g = 40*Q(physical[0])+27*Q(physical[1])+Q(root)
    fields = (2*Q(u[0]), 2*Q(u[1]), 2*g, Q(h[0]), Q(h[1]), Q(visible))
    return (*engine["field_coefficients"](), *(Q(9, 10)*x for x in fields))


def check_diagnostics(data, floor, a, s):
    from contextlib import redirect_stdout
    from io import StringIO

    def lower(compositions):
        with redirect_stdout(StringIO()):
            return engine["obstruction"](data, Q(a), Q(s), floor_fields=floor,
                                         compositions=compositions, adjust_heads=adjust_heads)
    triples = lower(((0, 3), (1, 2)))
    if (a, s) == (".567", ".685"):
        _, _, saved_volume = edge_gains(Q(a), Q(s), tuple(floor))
        assert min(saved_volume.values()) >= 59190
    four_terms = lower(last["FOUR_AND_THREE"])
    independent = lower(((1, 4),))
    degrees = [lower(tuple((v, n-v) for v in range(n))) for n in range(5, 11)]
    print(f"G3 con (a,s)=({a},{s}). Terna de tipos 3C y 1V+2C {float(triples):.10f},"
          f" cuatro términos {float(four_terms):.10f}, sólo 1V+4C {float(independent):.10f}.")
    print("  Minorantes por orden de cinco a diez:", *(f"{float(x):.8f}" for x in degrees), flush=True)
    assert 0 < triples < four_terms < 1 < independent and min(degrees) > 1
    return triples, four_terms, independent


def check_global(caps):
    engine["check_local"]()
    data, circuit = engine["v21"]["previous_tables"]()
    assert data[2] == circuit[3] and engine["F5"]["RM"] == Q(".001")
    discounts, pair = last["cross"]["coverage_discounts"](engine, circuit[4], data[2])
    data = (*data[:2], discounts, data[3])
    circuit = (*circuit[:3], discounts, *circuit[4:])
    data, circuit = last["improved_tables"](data, circuit, last["check_local"](), False)
    assert (data[3][1][0, 3][4], data[3][0][1, 2][3], data[3][1][1, 2][3]) == (2916, 17577, 62496)
    data, circuit = tighten(data, circuit, caps)
    for case, old in zip(CASES, last["CASES"][:3]):
        assert Q(case[1]) > Q(old[1]) and case[2] == old[2]
        engine["verify_model"](case, data, circuit, floor_fields=reference_fields(old),
                               disjoint_pair_factors=pair, adjust_heads=adjust_heads)
    old_floor = last["cross"]["overlap_floor"](engine, engine["CASES"][0])
    first = check_diagnostics(data, old_floor, ".567", ".685")
    second = check_diagnostics(data, old_floor, ".569", ".673")
    assert first[1] < Q(".88618") and first[2] > Q("1.65150")
    assert second[1] < Q(".90481") and second[2] > Q("1.60553")
    check_diagnostics(data, reference_fields(last["CASES"][0]), CASES[0][3], CASES[0][4])
    print("Filas completas y diagnósticos comprobados con fracciones.", flush=True)


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para mantener las comprobaciones.")
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--global", dest="global_row", action="store_true",
                        help="reconstruir las tablas de entrada y verificar los tres modelos")
    args = parser.parse_args()
    caps = geometry_caps()
    assert caps == {((1, 0, 3), 4): 432, ((0, 1, 2), 3): 5022, ((1, 1, 2), 3): 17904}
    if args.global_row:
        check_global(caps)

