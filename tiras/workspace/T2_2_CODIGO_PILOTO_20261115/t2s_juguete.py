#!/usr/bin/env python3
"""T2-2: modelos de juguete para las autopruebas y la prueba de humo. Sin decoder ni datos del proyecto.

Un juguete tiene NS sitios, cada uno de un tipo de una lista corta; el tipo fija sus clases y su ley
  l_r(s, c) = alfa_c r^e_c (c >= 1),   l_r(s, 0) = 1 - sum_c alfa_c r^e_c,
con alfa_c racional y e_c en {1, 2}. Cada clase tiene una respuesta (e, sindrome) pseudoaleatoria fijada por la
semilla del juguete. F (las "silenciosas"): al menos 3 sitios con fallo y h(e, sindrome) = 0 mod `modulo`, con h
un resumen entero fijo. Asi el atajo (menos de 3 sitios nunca es silenciosa) vale por construccion, como V5.
P_sil(r) se calcula exactamente por enumeracion (fracciones) cuando el espacio es pequeno.
"""
import hashlib
import itertools
import random
from fractions import Fraction as F

ALPHAS = (F(1, 2), F(3, 10), F(1, 5), F(7, 10), F(2, 5))


def _rng(seed, tag):
    return random.Random(int.from_bytes(hashlib.sha256(("juguete|%s|%s" % (seed, tag)).encode()).digest(), "big"))


class Juguete:
    def __init__(self, spec):
        self.spec = dict(spec)
        NS, ntypes, seed = spec["NS"], spec["tipos"], spec["semilla"]
        self.mod = spec.get("modulo", 5)
        self.maxcls = spec.get("max_clases", 3)
        g = _rng(seed, "tipos")
        types = []
        if spec.get("tipos_fijos"):
            types = [[(F(a), int(e)) for a, e in t] for t in spec["tipos_fijos"]]
            ntypes = len(types)
        else:
            for _ in range(ntypes):
                n = g.randint(1, self.maxcls)
                types.append([(ALPHAS[g.randrange(len(ALPHAS))], g.choice((1, 1, 2))) for _ in range(n)])
        self.types = types
        self.site_type = [g.randrange(ntypes) for _ in range(NS)]
        self.NS = NS
        self.ncls = [len(types[t]) for t in self.site_type]
        bits_e, bits_s = spec.get("bits_e", 4), spec.get("bits_s", 6)
        self.resp = [[(0, 0)] + [(g.getrandbits(bits_e) | 1, g.getrandbits(bits_s)) for _ in range(n)]
                     for n in self.ncls]
        self.calls = 0
        self.explicit = None
        if spec.get("silenciosas") is not None:
            self.explicit = set()
            for pairs in spec["silenciosas"]:
                cls = [0] * NS
                for s, c in pairs:
                    cls[s] = c
                self.explicit.add(tuple(cls))

    def laws(self, r):
        r = F(r)
        out = []
        for t in self.site_type:
            row = [a * r ** e for a, e in self.types[t]]
            out.append([1 - sum(row)] + row)
        return out

    def silent_cls(self, cls):
        nf = sum(1 for c in cls if c)
        if nf < 3:
            return False
        if self.explicit is not None:
            return tuple(cls) in self.explicit
        e = y = 0
        for s, c in enumerate(cls):
            if c:
                ee, yy = self.resp[s][c]
                e ^= ee
                y ^= yy
        return self._h(e, y)

    def _h(self, e, y):
        return ((e * 1000003 + y * 7919 + 12345) % 2147483647) % self.mod == 0

    def oracle(self, e, syn, cls):
        self.calls += 1
        if sum(1 for c in cls if c) < 3:
            return False
        if self.explicit is not None:
            return tuple(cls) in self.explicit
        return self._h(e, syn)

    def anchor_sample(self, rng, Ba, tick):
        """Muestreo directo por sitio en r_0 (umbrales float de la ley exacta)."""
        r0 = F(self.spec["r0"])
        cum = []
        for row in self.laws(r0):
            acc, c = F(0), []
            for p in row:
                acc += p
                c.append(float(acc))
            c[-1] = 1.0
            cum.append(c)
        rnd = rng.random
        for h in range(1, Ba + 1):
            if not h % 4096:
                tick()
            pairs = []
            for s in range(self.NS):
                u = rnd()
                cs = cum[s]
                if u < cs[0]:
                    continue
                for c in range(1, len(cs)):
                    if u < cs[c]:
                        pairs.append((s, c))
                        break
            if pairs:
                yield h, pairs

    # -------------------------------------------------------------- referencias exactas por enumeracion
    def states(self):
        return itertools.product(*[range(n + 1) for n in self.ncls])

    def exact(self, rates, kmax=4):
        """P_sil(r) y masas de F con |E| <= k (k = 3..kmax) exactas."""
        sil = [cls for cls in self.states() if self.silent_cls(cls)]
        out = {}
        for r in rates:
            L = self.laws(r)
            tot = F(0)
            byk = {}
            for cls in sil:
                p = F(1)
                for s, c in enumerate(cls):
                    p *= L[s][c]
                tot += p
                w = sum(1 for c in cls if c)
                byk[w] = byk.get(w, F(0)) + p
            out[r] = {"P": tot, "L_ing": {k: sum((v for w, v in byk.items() if w <= k), F(0))
                                          for k in range(3, kmax + 1)}}
        return out, sil


def make(spec):
    return Juguete(spec)
