"""Exact minimum number of paid data sites in a nontrivial X logical."""

from collections import Counter
from functools import reduce
from itertools import combinations
from operator import xor


def set_bits(word):
    while word:
        bit = word & -word
        yield bit.bit_length()-1
        word ^= bit


class LogicalGap:
    def __init__(self, seed_rows, checks, tags):
        self.columns = [(s << 9) | t for s, t in zip(checks, tags, strict=True)]
        assert all(t.bit_length() <= 9 for t in tags)
        classical = [w for w in range(1, 1 << 12)
                     if all((w & row).bit_count() % 2 == 0 for row in seed_rows)]
        assert Counter(w.bit_count() for w in classical) == {6: 6, 12: 1}
        for r in range(12):
            for w in classical:
                labelled = reduce(xor, (self.columns[12*r+j] for j in set_bits(w)), 0)
                assert 0 < labelled < 512
        self.row_cost = [min((w & ~free).bit_count() for w in classical)
                         for free in range(1 << 12)]
        self.witness = next(w for w in classical if w.bit_count() == 6)
        self.counts = Counter()

    def __call__(self, free):
        # For a logical (A,B), some seed kernel word b has b^T A a nonzero
        # seed codeword. Each column of that word needs at least one VV site.
        union = 0
        upper = 6
        for r in range(12):
            row = (free >> (12*r)) & 4095
            union |= row
            upper = min(upper, self.row_cost[row])
        # If B=0 each nonzero row of A is a seed codeword and the best
        # single row attains the minimum. If B!=0 the total weight is >=7
        # and the CC part costs >=1 unless some CC site is free.
        nonzero_cc_lower = max(self.row_cost[union] + int(free >> 144 == 0),
                               7-free.bit_count())
        lower = min(upper, nonzero_cc_lower)
        if lower == upper:
            self.counts["geometric_bounds_agree"] += 1
            return upper
        assert upper <= 5
        self.counts["quotient_search"] += 1
        return self.quotient_gap(free, upper, lower)

    def quotient_gap(self, free, upper, lower=0):
        # Divide the syndrome and its associated logical label by the
        # free columns. A zero syndrome with a nonzero tag is a logical.
        pivots = {}
        for j in set_bits(free):
            word = self.columns[j]
            while word >= 512:
                bit = 1 << (word.bit_length()-1)
                prev = pivots.get(bit)
                if prev is None:
                    pivots[bit] = word
                    break
                word ^= prev
            else:
                if word:
                    return 0
        if upper <= 1:
            return upper
        basis = sorted(pivots.items(), reverse=True)
        reduced = []
        for j, word in enumerate(self.columns):
            if (free >> j) & 1:
                continue
            for bit, pivot in basis:
                if word & bit:
                    word ^= pivot
            if 0 < word < 512:
                return 1
            if word:
                reduced.append(divmod(word, 512))
        if upper <= 2:
            return upper
        singles = {}
        for syndrome, tag in reduced:
            if singles.setdefault(syndrome, tag) != tag:
                return 2
        if upper <= 3:
            return upper
        assert upper in (4, 5)
        # A pair colliding with a single (or another pair) with a different
        # tag gives a logical of cost at most 3 (or 4). Repeated sites cancel.
        pair_tags = {}
        best = upper
        for (s, t), (other_s, other_t) in combinations(singles.items(), 2):
            synd, label = s ^ other_s, t ^ other_t
            single_tag = singles.get(synd)
            if single_tag is not None and single_tag != label:
                return 3
            if best > 4 and pair_tags.setdefault(synd, label) != label:
                best = 4
                if lower >= 4:
                    return best
        return best
