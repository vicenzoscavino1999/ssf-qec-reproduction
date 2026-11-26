"""Expand the published tie-relaxed (3,1) census and restore the physical SSF order."""

from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from math import comb, factorial
from runpy import run_path


def independent_decoder(compare, checks):
    audit = run_path(str(compare.COMP/"auditar_esquinas_HGP_225.py"))
    _, independent_checks, generators = audit["incidence"]()
    assert independent_checks == checks
    tags, decode, _ = audit["decoder"](independent_checks, generators)
    return tags, decode


def single_readout_triples(compare, h, checks, tags, ssf):
    affine = run_path(str(compare.COMP / "comparar_HGP_225.py"))
    maps = affine["affine_maps"](checks[:12])
    classes = compare.CIRC["terminal_classes"](compare.CIRC["all_winners"](h, checks, tags))

    @lru_cache(None)
    def transforms(nvv, labels, axis, readout_coordinate):
        return {(tuple(mapping[int(i >= nvv)][v] for i, v in enumerate(labels)),
                 mapping[axis][readout_coordinate]) for mapping in maps}

    candidates = set()
    all_triples = all_possible = 0
    for nvv in range(4):
        axis = affine["axis_orbits"](nvv, maps)
        ordered = possible = 0
        for (left, nl), (right, nr) in product(axis, repeat=2):
            def sites(rows, cols):
                return tuple(12*r+c if i < nvv else 144+9*r+c
                             for i, (r,c) in enumerate(zip(rows,cols)))
            u,v,w = sites(left, right)
            if len({u,v,w}) < 3:
                continue
            ordered += nl*nr
            clean, tag = checks[u]^checks[v]^checks[w], tags[u]^tags[v]^tags[w]
            for j in range(108):
                d = 1 << j
                if not any(t != tag for t in classes(clean ^ d, d)):
                    continue
                possible += nl*nr
                row, col = divmod(j, 9)
                for (a, r), (b, c) in product(transforms(nvv,left,0,row),
                                             transforms(nvv,right,1,col)):
                    candidates.add((*sorted(sites(a,b)), 9*r+c))
        permutation = factorial(nvv)*factorial(3-nvv)
        assert ordered % permutation == possible % permutation == 0
        all_triples += ordered//permutation
        all_possible += possible//permutation
    assert all_triples == comb(225,3) and all_possible == len(candidates) == 101088
    silent = Counter()
    independent_tags, independent_ssf = independent_decoder(compare, checks)
    for u,v,w,j in sorted(candidates):
        d = 1 << j
        observed = checks[u]^checks[v]^checks[w]^d
        terminal, correction = ssf(observed)
        is_silent = (terminal == d and correction != tags[u]^tags[v]^tags[w])
        end, residual = independent_ssf(observed)
        independent = (end == d and residual != independent_tags[u]^independent_tags[v]^independent_tags[w])
        assert independent == is_silent
        if is_silent:
            silent[observed] += 1
    return silent


def two_data_three_readouts(patterns, checks, tags, ssf):
    # The row lemma leaves one prescribed inversion in the data row;
    # the other two readout faults must be at distinct exterior checks.
    outside = [tuple(a ^ b for a,b in combinations(
        [1 << j for j in range(108) if j//9 != r], 2)) for r in range(12)]
    candidates = hits = 0
    for r,u,v,d in patterns:
        original = checks[u]^checks[v]^d
        tag = tags[u]^tags[v]
        for ext in outside[r]:
            terminal, correction = ssf(original^ext)
            hits += (terminal == (d^ext) and correction != tag)
        candidates += len(outside[r])
    assert candidates == 324*comb(99,2)
    return candidates, hits
