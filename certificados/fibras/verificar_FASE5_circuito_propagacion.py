"""HGP (3,4), una ronda. Transferencia circuital uniforme a r_c=4.2e-8."""

from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations, product
from pathlib import Path
from runpy import run_path


HERE = Path(__file__).resolve().parent
RC, RM = Q("4.2e-8"), Q(".001")
REFERENCE_P, REFERENCE_Q = Q("3.2e-7"), Q(".001")
PAIRS = ((0, 1), (0, 2), (1, 2))


def previous(name):
    return run_path(str(HERE / f"verificar_H5A_3_4_{name}.py"))


def circuit_weights():
    # El bit independiente t se agrega al check físico de esta CNOT.
    for rate in (Q(0), RC, Q("3.75e-7"), Q(1)):
        control = 8*rate/15
        target = 4*rate/(15-8*rate)
        law = Counter()
        for a, t, u in product(range(2), repeat=3):
            mass = (control if a else 1-control)*(target if t else 1-target)/2
            law[a, t ^ (a*u)] += mass
        assert law[0, 0] == 1-4*rate/5
        assert all(law[bits] == 4*rate/15 for bits in ((0, 1), (1, 0), (1, 1)))

    target = 4*RC/(15-8*RC)
    data_odds = (1-2*RC/3)**(-9)*(1-8*RC/15)**(-3)-1
    check_odds = (1-2*RC/3)**(-5)/(1-RM)/(1-target)**7-1
    assert data_odds < Q("3.192001e-7") < REFERENCE_P/(1-2*REFERENCE_P)
    assert check_odds < Q(".001001220") < REFERENCE_Q/(1-2*REFERENCE_Q)

    proposed = Q("3.75e-7")
    odd = (1-(1-4*proposed/3)**9*(1-16*proposed/15)**3)/2
    assert odd > Q("2.84999e-6")
    assert odd/(1-odd) > Q("7.59")*proposed/(1-2*proposed)
    return data_odds, check_odds, odd


def partitions(size):
    """Palabras de crecimiento restringido, una por partición de vértices."""
    if size == 1:
        yield (0,)
        return
    for prefix in partitions(size-1):
        for block in range(max(prefix)+2):
            yield prefix+(block,)


def block_profile(types, labels):
    return tuple(sorted(
        (sum(t == 0 and b == i for t, b in zip(types, labels)),
         sum(t == 1 and b == i for t, b in zip(types, labels)))
        for i in set(labels)))


def check_block_subsets(contacts, profiles, stable, prior_caps):
    @lru_cache(None)
    def cap(groups):
        nv = sum(v for v, c in groups)
        nc = sum(c for v, c in groups)
        if nv+nc == 2:
            # Arista fuerte o relación relajada dentro de un bloque. Entre
            # bloques la pareja es débil, con techos 0, 1, 0 según su tipo.
            return 2 if len(groups) == 1 else (0, 1, 0)[nc]
        types, labels = [], []
        for block, (v, c) in enumerate(groups):
            types.extend((0,)*v + (1,)*c)
            labels.extend((block,)*(v+c))
        gv, gc = groups.count((1, 0)), groups.count((0, 1))
        bulk = 7*nv+10*nc-2*gc
        prior = 84*nv+118*nc-11*gv-28*gc
        direct = contacts["group_bound"](types, labels, profiles, stable)
        sixths = min(prior, 6*prior_caps[nv, nc], 6*direct)
        return min(bulk//6, max(0, (sixths-1)//66))

    seen, count = set(), 0
    for size in (3, 4):
        for types, labels in product(product((0, 1), repeat=size), partitions(size)):
            group = block_profile(types, labels)
            if group in seen:
                continue
            seen.add(group)
            for length in range(2, size):
                for keep in combinations(range(size), length):
                    for sublabels in partitions(length):
                        if any(sublabels[i] == sublabels[j] and labels[keep[i]] != labels[keep[j]]
                               for i, j in combinations(range(length), 2)):
                            continue
                        sub = block_profile([types[i] for i in keep], sublabels)
                        assert cap(sub) <= cap(group), (sub, group)
                        count += 1
    assert (len(seen), count) == (47, 643)
    return count


def pair_cap(i, j, relation):
    if relation == 0:
        return 0
    if i != j:
        return {1: 1, 2: 1, 4: 2}[relation]
    rows = ((0, 0, 1, 2), (0, 1, 2, 2))
    return rows[i][relation-1]


def check_triple_subsets(isolated, profiles, stable):
    seen = set()
    for types in product((0, 1), repeat=3):
        options = [(0,)+tuple(isolated["positions"](types[i], types[j])) for i, j in PAIRS]
        for edges in product(*options):
            if sum(bool(edge) for edge in edges) < 2:
                continue
            rel = dict(zip(PAIRS, edges))
            key = min((tuple(types[i] for i in order),
                       tuple(rel[tuple(sorted((order[i], order[j])))] for i, j in PAIRS))
                      for order in permutations(range(3)))
            if key in seen:
                continue
            seen.add(key)
            bound = isolated["triple_budget"](types, edges, profiles, stable)
            for (i, j), edge in rel.items():
                assert pair_cap(types[i], types[j], edge) <= bound, (types, edges, i, j)
    assert len(seen) == 144
    return len(seen)


def check_donors_and_coverage(far, geometry):
    @lru_cache(None)
    def raw(nv, nc, donors, coverage):
        if nv+nc == 0:
            return 0
        return far["group_budget"](nv, nc, geometry, donors, guard=coverage)

    def cap(nv, nc, donors):
        if nv+nc < 2:
            return 0
        return raw(nv, nc, min(donors, nv+nc), False)

    def cover(nv, nc):
        return min(246*nv+340*nc, raw(nv, nc, nv+nc, True))

    n_caps, n_covers = 0, 0
    for size in range(1, far["LIMIT"]+1):
        for nv in range(size+1):
            nc = size-nv
            for v, c in product(range(nv+1), range(nc+1)):
                if v+c == size:
                    continue
                assert cover(v, c) <= cover(nv, nc), (v, c, nv, nc)
                n_covers += 1
                if size >= 3:
                    for donors in range(2, size+1):
                        assert cap(v, c, donors) <= cap(nv, nc, donors), (v, c, nv, nc, donors)
                        n_caps += 1
    assert (n_caps, n_covers) == (707, 182)
    return n_caps, n_covers


def main():
    data, checks, odd = circuit_weights()
    contacts, stripes, base, history, isolated, far = map(previous, (
        "contacto_de_checks", "franjas", "medidas_y_un_dato", "medidas_por_historia",
        "medidas_aisladas", "propagacion_lejana"))
    profiles, stable, reps = contacts["geometry"]()
    assert max(k.bit_count() for k, ok in enumerate(stable) if ok) == 3
    contacts["universal_discounts"](profiles, stable)
    history["far_credits"](profiles, stable)
    halos = contacts["exposure"](reps, base["budgets"])
    assert isolated["guard_cost"](reps, base["budgets"]) == (504, 628)
    _, caps = stripes["cost_tables"](stripes["light_profiles"](), previous("fondo_pasivo")["step_costs"],
                                      previous("pares_adyacentes")["stripe_profiles"])
    n_blocks = check_block_subsets(contacts, profiles, stable, caps)
    n_triples = check_triple_subsets(isolated, profiles, stable)
    geometry = far["cover_profiles"](reps)
    n_caps, n_covers = check_donors_and_coverage(far, geometry)
    print(f"Subconjuntos: {n_blocks} comparaciones de bloques, {n_triples} patrones de tríos, "
          f"{n_caps} comparaciones de saldo y {n_covers} de cobertura.", flush=True)

    # Reconstruir la fila completa con los testigos publicados en las notas.
    hist = contacts["positions"](profiles, stable, caps, stripes["coefficients"](4))
    tables = history["data_coefficients"](hist)
    triples = isolated["count_triples"](profiles, stable, contacts["rooted_trees"])
    series = isolated["group_series"](10)
    pairs = far["pair_counts"](far["winners"](reps), stable)
    trees = stripes["coefficients"](far["LIMIT"])
    a, s = Q(".5488178589126"), Q(".86")
    kernels = far["kernels"](history["ball"], a**6, s)
    poly, marked = far["group_polynomials"](a**6, history["ball"], tables, triples, series, pairs,
                                             geometry, trees)
    far["certify"](a, s, halos, series, poly, marked, kernels, base["exp_upper"])
    print(f"Circuito uniforme rc<=4.2e-8, rm<=1e-3. Razones superiores por fibra "
          f"{float(data):.12g} y {float(checks):.12g}.")
    print(f"En rc=3.75e-7 la paridad X impar en una fibra VV completa vale {float(odd):.12g}.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Ejecutar sin -O para verificar todas las desigualdades.")
    main()

