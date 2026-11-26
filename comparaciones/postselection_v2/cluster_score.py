"""Lee--English--Bartlett's published score for the observed input syndrome."""
from collections import defaultdict
from functools import lru_cache
from math import expm1, log1p, sqrt

import numpy as np
from scipy.sparse import csc_matrix
from ldpc_post_selection.bplsd_decoder import SoftOutputsBpLsdDecoder
from ldpc_post_selection.cluster_tools import compute_cluster_norm_fraction


def detector_priors(checks, h, schedule=None, rate=None):
    logs = defaultdict(float)

    def add(mask, probability):
        if mask:
            logs[mask] += log1p(-2*float(probability))

    for bit in range(108):
        add(1 << bit, h)
    if schedule is None:
        for v, mask in enumerate(checks):
            add(mask | (1 << (108+v)), rate)
    else:
        r = float(rate)
        gate = (8*r/15)/(1+sqrt(1-16*r/15))
        for v, row in enumerate(schedule):
            data = 1 << (108+v)
            future = sum(1 << bit for bit in row.values())
            for time in range(12):
                if time in row:
                    bit = 1 << row[time]
                    future ^= bit
                    for word in (data | future, bit, data | (future ^ bit)):
                        add(word, gate)
                else:
                    add(data | future, 2*r/3)
            assert not future
    return {mask: -expm1(total)/2 for mask, total in sorted(logs.items())}


def scorer(priors):
    words = list(priors)
    matrix = np.array([[(word >> bit) & 1 for word in words] for bit in range(108)],
                      dtype=np.uint8)
    decoder = SoftOutputsBpLsdDecoder(H=csc_matrix(matrix), p=list(priors.values()),
                detector_time_coords=[0]*108, max_iter=30, bp_method="product_sum",
                lsd_method="LSD_0", lsd_order=0, omp_thread_count=1, schedule="parallel")

    @lru_cache(maxsize=60000)
    def score(word):
        if not word:
            return 0.0
        syndrome = np.unpackbits(np.frombuffer(word.to_bytes(14, "little"), dtype=np.uint8),
                                 bitorder="little")[:108]
        prediction, _, _, result = decoder.decode(syndrome)
        assert np.array_equal(decoder.H @ prediction % 2, syndrome)
        return float(compute_cluster_norm_fraction(result["cluster_llrs"], 2))
    return score


def check_circuit_priors(priors, schedule, rc, rm):
    """Compare every (record, final-X) column with an independently propagated Stim circuit."""
    import stim
    from ldpc_post_selection.stim_tools import dem_to_parity_check
    n,m = 225,108
    circuit = stim.Circuit()
    circuit.append("R", range(n+m))
    for layer in range(12):
        gates = [(v,n+row[layer]) for v,row in enumerate(schedule) if layer in row]
        active = [q for gate in gates for q in gate]
        if active:
            circuit.append("CX",active)
            circuit.append("DEPOLARIZE2",active,float(rc))
        circuit.append("DEPOLARIZE1",sorted(set(range(n+m))-set(active)),float(rc))
    circuit.append("M",range(n,n+m),float(rm))
    circuit.append("M",range(n))
    for j in range(m):
        circuit.append("DETECTOR",[stim.target_rec(j-n-m)])
    for v in range(n):
        circuit.append("OBSERVABLE_INCLUDE",[stim.target_rec(v-n)],v)
    record, final, probabilities = dem_to_parity_check(circuit.detector_error_model())
    rebuilt = {}
    for i,p in enumerate(probabilities):
        mask = sum(1<<int(j) for j in record.indices[record.indptr[i]:record.indptr[i+1]])
        mask |= sum(1<<(m+int(j)) for j in final.indices[final.indptr[i]:final.indptr[i+1]])
        assert mask not in rebuilt
        rebuilt[mask] = p
    assert rebuilt.keys() == priors.keys()
    assert max(abs(rebuilt[s]-v) for s,v in priors.items()) < 5e-15
