"""Bernoulli input catalogues and the published BP+LSD cluster assignments."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import ceil, comb
from pathlib import Path
from time import perf_counter
import gzip
import hashlib
import json
import sys

import numpy as np
from scipy.sparse import csc_matrix
from ldpc_post_selection.bplsd_decoder import SoftOutputsBpLsdDecoder
from ldpc_post_selection.cluster_tools import compute_cluster_norm_fraction
from extra_silent import single_readout_triples

N, M = 225, 108
P, Q = Fraction(1, 10**6), Fraction(1, 10**4)
BITS = tuple(1 << j for j in range(M))
LEE_CUT = 3368
RECORD_BYTES = 14 + 29 + 1


def load_comparison(repo):
    package = repo.resolve() / "comparaciones" / "postselection_v2"
    if not (package / "compare_rules.py").is_file():
        raise ValueError("Pass the root of the reproduction package with --repo.")
    sys.path.insert(0, str(package))
    import compare_rules
    return compare_rules


def input_catalogues(compare):
    original = compare.CIRC["pilot"]
    seed_rows, checks, tags, _, _, ssf = original["build"]()
    patterns = []
    assert original["row_sectors"](checks, tags, ssf, patterns=patterns) == 324
    triples, _ = compare.CIRC["physical_triples"](seed_rows, checks, tags, ssf)
    assert len(triples) == 1656
    silent = {(2, 1): Counter(), (2, 2): Counter(), (3, 0): Counter()}
    for row, u, v, wrong in patterns:
        word = checks[u] ^ checks[v] ^ wrong
        silent[2, 1][word] += 1
        for j, b in enumerate(BITS):
            if j // 9 != row:
                silent[2, 2][word ^ b] += 1
    for (u, v, w), wrong in triples:
        if not wrong:
            silent[3, 0][checks[u] ^ checks[v] ^ checks[w]] += 1
    assert [c.total() for c in silent.values()] == [324, 32076, 36]
    silent[3, 1] = single_readout_triples(compare, seed_rows, checks, tags, ssf)
    marginal = {
        (0, 0): Counter([0]),
        (0, 1): Counter(BITS),
        (0, 2): Counter(a ^ b for a, b in combinations(BITS, 2)),
        (0, 3): Counter(a ^ b ^ c for a, b, c in combinations(BITS, 3)),
        (1, 0): Counter(checks),
        (1, 1): Counter(v ^ b for v in checks for b in BITS),
        (2, 0): Counter(u ^ v for u, v in combinations(checks, 2)),
    }
    for (a, b), hist in marginal.items():
        assert hist.total() == comb(N, a) * comb(M, b)
    return seed_rows, checks, tags, marginal, silent, triples, ssf, patterns


def words_for(marginal, silent):
    return sorted(set().union(*marginal.values(), *silent.values()))


def sector_weight(a, b):
    return P**a * (1-P)**(N-a) * Q**b * (1-Q)**(M-b)


def make_features(compare, checks, words, cache):
    priors = compare.detector_priors(checks, Q, rate=P)
    columns = list(priors)
    assert columns[:M] == list(BITS)
    assert columns[M:] == [(1 << (M+j)) | h for j, h in enumerate(checks)]
    matrix = np.array([[(word >> bit) & 1 for word in columns] for bit in range(M)],
                      dtype=np.uint8)
    decoder = SoftOutputsBpLsdDecoder(
        H=csc_matrix(matrix), p=list(priors.values()), detector_time_coords=[0]*M,
        max_iter=30, bp_method="product_sum", lsd_method="LSD_0", lsd_order=0,
        omp_thread_count=1, schedule="parallel")
    started = perf_counter()
    temp = cache.with_name(cache.name + ".part")
    with gzip.open(temp, "wb", compresslevel=2) as f:
        for idx, word in enumerate(words, 1):
            encoded = word.to_bytes(14, "little")
            if word:
                s = np.unpackbits(np.frombuffer(encoded, dtype=np.uint8),
                                  bitorder="little")[:M]
                correction, _, _, soft = decoder.decode(s)
                assert np.array_equal(decoder.H @ correction % 2, s)
                assert not np.any(correction & (soft["clusters"] == 0))
                memberships = soft["clusters"][M:] != 0
                erased = np.packbits(memberships, bitorder="little").tobytes()
                score = compute_cluster_norm_fraction(soft["cluster_llrs"], 2)
                lee_discards = ceil(10**6 * score) > LEE_CUT
            else:
                erased, lee_discards = bytes(29), False
            f.write(encoded + erased + bytes((lee_discards,)))
            if idx % 50000 == 0:
                print(f"{idx}/{len(words)} input syndromes; {perf_counter()-started:.1f} s",
                      flush=True)
    temp.replace(cache)
    return perf_counter() - started


def collect(repo, cache):
    started = perf_counter()
    compare = load_comparison(repo)
    _, checks, _, marginal, silent, *_ = input_catalogues(compare)
    words = words_for(marginal, silent)
    geometry_seconds = perf_counter() - started
    feature_seconds = make_features(compare, checks, words, cache)
    meta = dict(unique_syndromes=len(words), geometry_seconds=geometry_seconds,
                feature_seconds=feature_seconds, elapsed_seconds=perf_counter()-started,
                sha256=hashlib.sha256(cache.read_bytes()).hexdigest())
    cache.with_name(cache.name + ".json").write_text(json.dumps(meta, indent=2)+"\n")
    print(json.dumps(meta), flush=True)


def read_features(cache, words):
    meta = json.loads(cache.with_name(cache.name + ".json").read_text())
    assert hashlib.sha256(cache.read_bytes()).hexdigest() == meta["sha256"]
    assert len(words) == meta["unique_syndromes"]
    with gzip.open(cache, "rb") as f:
        data = f.read()
    assert len(data) == len(words)*RECORD_BYTES
    features = {}
    for offset, word in enumerate(words):
        i = offset*RECORD_BYTES
        assert int.from_bytes(data[i:i+14], "little") == word
        erased = int.from_bytes(data[i+14:i+43], "little")
        lee = data[i+43]
        assert erased.bit_length() <= N and lee in (0, 1)
        features[word] = (erased, lee)
    return features, meta


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    args = parser.parse_args()
    collect(args.repo, args.cache)
