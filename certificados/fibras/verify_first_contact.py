"""Exact checks for Lemma 14. The degree-uniform proofs are in the manuscript."""

from fractions import Fraction
from functools import cache


def check_degrees():
    pairs = 0
    for a in range(1, 33):
        for b in range(max(a, 2), 2 * a + 9):
            minimum = min(2 * a, a + b - 2)
            multiple = []
            for x in range(b + 1):
                for y in range(a + 1):
                    weight = x + y
                    area = x * (a - y) + y * (b - x)
                    if area == 0 or 2 * area < a * weight:
                        continue
                    assert area >= a
                    if weight > 1:
                        multiple.append(area)
                    if a >= 3 and y > 0:
                        assert area >= b
            assert min(multiple) == minimum
            if a >= 3:
                check_short_rows(a, b)
            pairs += 1
    return pairs


def check_short_rows(a, b):
    # Any admitted area below b uses rows only. A positive drop has its
    # area's parity. Compare integer drops with the untouched CC singleton.
    for t in range(1, b // 2 + 1):
        balance = b - 2 * t
        for area in range(a, b, a):
            for drop in range(1 if area % 2 else 2, area + 1, 2):
                if b * drop < area * balance:
                    continue
                if b <= a + 2:
                    assert drop >= balance
                in_wider_band = 2 * b < 3 * a or (a % 2 and b < 2 * a)
                if in_wider_band:
                    assert area == a
                    max_lights_before_contact = 1 + a // 2
                    if (area + drop) // 2 <= max_lights_before_contact:
                        assert drop >= balance


def candidates(a, b):
    # In a b by a parent rectangle the first b sites select rows, the next
    # a sites select columns. An all-ones seed has this same rectangle at every parent.
    rows = [((1 << a) - 1) << (a * i) for i in range(b)]
    cols = [sum(1 << (a * i + j) for i in range(b)) for j in range(a)]
    footprints = [0]
    for line in rows + cols:
        footprints += [old ^ line for old in footprints]
    return {action: f for action, f in enumerate(footprints)
            if f and 2 * f.bit_count() >= a * action.bit_count()}


def score(syndrome, footprint):
    return Fraction(syndrome.bit_count() - (syndrome ^ footprint).bit_count(),
                    footprint.bit_count())


def check_witness(a, b, background_cells, scores, expected_balance):
    moves = candidates(a, b)
    row, column = 1, 1 << b
    mixed = row | column
    # One fixed order. Every physical generator has the same choices of
    # footprints, so a preferred representative of each first move suffices.
    order = [row, mixed] + [f for f in moves if f not in (row, mixed)]
    d = sum(1 << (a * i + j) for i, j in background_cells)
    s = d ^ moves[column]
    balance = b - 2 * (moves[column] & d).bit_count()
    assert balance == expected_balance
    first_drop = s.bit_count() - (s ^ moves[row]).bit_count()
    assert first_drop < balance
    beta = Fraction(1, 100)
    for expected, value in zip((row, mixed), scores[:2]):
        winner = max(order, key=lambda f: score(s, moves[f]))
        assert winner == expected and score(s, moves[winner]) == value >= beta
        s ^= moves[winner]
    terminal = max(score(s, f) for f in moves.values())
    assert s == d and terminal == scores[2] < beta
    print(f"({a},{b}) scores {scores[0]}, {scores[1]}, {terminal}. "
          f"First drop {first_drop} < {balance}. Final syndrome equals D.")


def band_seed(k):
    """Column supports of the actual binary seed, including its regular completion."""
    a, b = 6 + k, 9 + k
    common = ("c", "d1", "d2", "d3", "e1", "e2")
    far = [tuple(f"p{j}_{i}" for i in range(a - 2)) for j in range(k)]
    raw = {"L": common + tuple(f"u{j}" for j in range(k)),
           "M": common + tuple(f"w{j}" for j in range(k))}
    remote = []
    for j, block in enumerate(far):
        raw[f"H{j}"] = block + (f"u{j}", f"w{j}")
        raw[f"R{j}"] = block + ("z1", "z2")
        remote.append(f"R{j}")
    if k == 1:
        shared = tuple(f"t{i}" for i in range(a - 1))
        for j in (1, 2):
            raw[f"T{j}"] = shared + (f"z{j}",)
            remote.append(f"T{j}")
    cells = sorted(set().union(*raw.values()))
    assert len(raw) == 6 and len(cells) == 3 * a
    assert all(sum(cell in edge for edge in raw.values()) == 2 for cell in cells)
    copies = range(2 * a)
    seed = [tuple((t, cell) for cell in edge) for t in copies for edge in raw.values()]
    # A distinct pair of complementary intervals for each shift, plus even/odd copies.
    halves = [{(start + j) % (2 * a) for j in range(a)} for start in range(a)]
    halves.append(set(copies[::2]))
    for cell in cells:
        for half in halves:
            for part in (half, set(copies) - half):
                seed.append(tuple((t, cell) for t in sorted(part)))
    rows = {}
    for index, column in enumerate(seed):
        for cell in column:
            rows.setdefault(cell, set()).add(index)
    assert all(len(set(col)) == a for col in seed)
    assert all(len(row) == b for row in rows.values())
    assert len({frozenset(col) for col in seed}) == len(seed)
    assert len({frozenset(row) for row in rows.values()}) == len(rows)
    assert (len(rows), len(seed)) == (6 * a * a, 6 * a * b)
    return seed, rows, {name: i for i, name in enumerate(raw)}, remote


@cache
def footprints_by_area(a, b):
    by_area = {}
    for f in candidates(a, b).values():
        by_area.setdefault(f.bit_count(), set()).add(f)
    return by_area


@cache
def best_score(a, b, row_masks):
    """Maximize over all admitted subsets in one B-by-A generator rectangle."""
    top, bottom = -1, 1
    for columns in range(1 << a):
        y = columns.bit_count()
        inside = [(row & columns).bit_count() for row in row_masks]
        lights = sum(inside)
        gains = sorted((row.bit_count() - 2 * n for row, n in zip(row_masks, inside)),
                       reverse=True)
        for x in range(b + 1):
            area = a * x + b * y - 2 * x * y
            if area and 2 * area >= a * (x + y):
                drop = 2 * lights - area
                if drop * bottom > top * area:
                    top, bottom = drop, area
            if x < b:
                lights += gains[x]
    best = Fraction(top, bottom)
    # Independently enumerate every eligible subset's actual XOR footprint.
    state = sum(mask << (i * a) for i, mask in enumerate(row_masks))
    direct = max(Fraction(2 * max((state & f).bit_count() for f in shapes) - area, area)
                 for area, shapes in footprints_by_area(a, b).items())
    assert best == direct
    return best


def global_score(state, seed, rows, a, b):
    # Visit every physical parent up to identical projections of its row and
    # column indices on the nonzero syndrome. Pad omitted rows/columns with zeros.
    active_i = {i for i, _ in state}
    active_j = {j for _, j in state}
    first = {tuple(sorted(active_i & neighbors)) for neighbors in rows.values()}
    second = {tuple(j for j in col if j in active_j) for col in seed}
    result = Fraction(-1)
    for indices in first:
        for cells in second:
            masks = [sum(1 << t for t, j in enumerate(cells) if (i, j) in state)
                     for i in indices]
            masks.extend([0] * (b - len(indices)))
            result = max(result, best_score(a, b, tuple(sorted(masks))))
    return result


def check_band(k):
    a, b = 6 + k, 9 + k
    seed, rows, index, remote = band_seed(k)
    c = (0, "c")
    i = index["L"]
    assert i in rows[c]
    column = {(j, c) for j in rows[c]}

    def row(name):
        return {(i, cell) for cell in seed[index[name]]}

    # Any parent that can contain or read v and a remote fault must meet both
    # footprints. Its column set would contain c and a remote check.
    near = set().union(*(col for col in seed if c in col))
    assert all(near.isdisjoint(seed[index[name]]) for name in remote)
    d = {(i, (0, f"d{j}")) for j in (1, 2, 3)}
    d.update((j, c) for j in sorted(rows[c] - {i})[:3])
    s0 = d ^ column
    for name in remote:
        s0 ^= row(name)
    balance = len(column) - 2 * len(column & d)
    prefix = [row(f"H{j}") for j in range(k)]
    target = row("L")
    mixed = row("M") ^ column
    # Each preferred action belongs to the parent (c, its column). The first
    # k+1 are singletons. The last selects one row and the c column together.
    assert c in seed[index["M"]] and len(mixed) == a + b - 2
    assert all(f.isdisjoint(column) for f in prefix)
    assert target & column == {(i, c)}
    order = prefix + [target, mixed]
    s = s0.copy()
    weights, maxima = [len(s)], []
    for pos, move in enumerate(order):
        best = global_score(s, seed, rows, a, b)
        drop = len(s) - len(s ^ move)
        assert Fraction(drop, len(move)) == best >= Fraction(1, b)
        assert all(Fraction(len(s) - len(s ^ old), len(old)) < best for old in order[:pos])
        if pos == k:
            new_lights = s - s0
            assert new_lights & target == {(i, (0, f"u{j}")) for j in range(k)}
            assert max(sum((j, cell) in new_lights for cell in col)
                       for j in rows[c] for col in seed if c in col) == k
            assert drop == a - 4 < balance == b - 6
            assert drop == (a * balance + b - 1) // b
        maxima.append(best)
        s ^= move
        weights.append(len(s))
    assert s == d
    last = global_score(d, seed, rows, a, b)
    assert last < 0
    expected = {1: ([15, 12, 9, 6], Fraction(1, 5), Fraction(-1, 7)),
                2: ([23, 19, 15, 11, 6], Fraction(5, 17), Fraction(-1, 4))}
    sizes, mixed_score, final_score = expected[k]
    assert weights == sizes
    assert maxima == [Fraction(a - 4, a)] * (k + 1) + [mixed_score] and last == final_score
    print(f"({a},{b}) seed {len(rows)} x {len(seed)}, weights {weights}. "
          f"Exact maxima {', '.join(str(x) for x in maxima + [last])}. "
          f"First contact {a-4} < {balance} after {k} new check(s) on the winning row.")


if __name__ == "__main__":
    if not __debug__:
        raise RuntimeError("Run without -O to enable the exact checks.")
    count = check_degrees()
    print(f"Area and integer bounds checked on {count} degree pairs.")
    check_witness(2, 3, [(0, 1)], (1, Fraction(1, 3), 0), 3)
    check_witness(2, 4, [(0, 1)], (1, Fraction(1, 2), 0), 4)
    check_witness(3, 6, [(4, 0), (5, 0), (0, 1)],
                  (Fraction(1, 3), Fraction(1, 7), Fraction(-1, 7)), 2)
    check_band(1)
    check_band(2)
