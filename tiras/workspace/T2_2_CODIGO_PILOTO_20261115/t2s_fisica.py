#!/usr/bin/env python3
"""T2-2, piloto m = 6 (BORRADOR_PLAN_T2_2.txt v3): modelo fisico TIRA_6x6 para t2s_motor.py.

Codigo preparado para revision (version 1). NO SE HA EJECUTADO: ni sus controles previos, ni el decoder, ni el
muestreo del ancla. Las autopruebas y la prueba de humo no lo importan. Solo correra con el plan sellado y una
orden del usuario, dentro de t2_2.py --corrida (bajo t2s_supervisor.py).

Contrato (1, 2.2, 2.4, 2.6, V3, V4, V5 del plan):
- Sin cambios en b4_tiras.py, t2_calibracion.py, t2_mc.py ni t2_cp.py: se importan tras comprobar su SHA-256
  (t2_mc.py comprueba a su vez t2_calibracion.py, RESULTADO_T2_0.json, b4_tiras.py y RESULTADO_B4.json).
- Instancia TIRA_6x6 de B4 (b4.make_inst), decoder SSF (ce.make_decoder), calendario, ubicaciones y clases de
  sitio de T2-0 (t2.schedule, t2.locations, t2.site_classes con restricted = False). Sitio s: clases no nulas
  CL[s] (respuesta global (e, sindrome)) con claves canonicas CK[s].
- Ley de sitio en la tasa r (r_m = r_c = r): t2.law_of(LTK[s], L, r, r), leida en la clave nula y en CK[s];
  su soporte debe ser exactamente {0} U CK[s] con probabilidades > 0 (V3; si no, se detiene).
- Evento (oracle): como t2_calibracion.classify y t2_mc.decode_event: cor = dec(sindrome), R = e ^ cor;
  silenciosa si ce.lin(cols, R) == 0 y R no esta en S_X (ce.inside(R, bz) falso). Cache del decoder vaciada al
  pasar de 300000 entradas, como en T2-1.
- Ancla: muestreo por ubicaciones de t2_mc.w_block (saltos geometricos por grupos de igual probabilidad, y
  categoria por ubicacion con fallo en orden creciente), con las marginales X de t2_mc.location_laws en
  r_c = r_m = r_0. Cada historia se traduce a sitios por XOR de las claves canonicas de las respuestas de su
  fibra (V4). Mismos supuestos de simulacion que D-mc (MT19937, 53 bits, inversion geometrica en coma
  flotante, umbrales en coma flotante): no se acota el error de ley.
- Controles previos (detienen): dependencias; matrices iguales a B4; 527 ubicaciones; NS = 91 y 409 clases;
  marginales X iguales a las formas cerradas; leyes de sitio desde las marginales = law_of en r_0 y en las
  cuatro tasas objetivo; propagate de cada fallo unico = resp_formula; V4 en historias aleatorias del ancla
  (propagate = XOR de las respuestas de sitio) con un flujo propio; fuentes de V5 (exhaustivo de T2-0 y filas
  k = 2 de T2-1). Ninguno llama al decoder.
"""
import json
import math
import os
import sys
from fractions import Fraction as F

import t2s_motor as tmo

T21 = "T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/"
T20 = "T2_0_CALIBRACION_CIRCUITO_TIRAS_20261113/"
B4D = "B4_TIRAS_ESCALADO_N_20261113/"
RES_T21 = "T2_1_CORRIDA_2_EJECUCION_20261115/corrida_2/RESULTADO_T2_1.json"
SHA = {
    T21 + "t2_mc.py": "fcf0323c280381c19c807f95e2f475ae670902d3254e8dfbf9e48aff3634a957",
    T21 + "t2_cp.py": "7e1a913b84ddb11941bcd8c09d26a7a37b23ced7173b51ccd0f9cd0c098f9565",
    T20 + "t2_calibracion.py": "195bf3317f7af114f4cdc1daad714e828534ddd57abab29c81daa7126d17cc61",
    T20 + "RESULTADO_T2_0.json": "fb56c4324034549476813bd1a0d2bacd6f46c27d552fbf71ba2d0813d62df38e",
    B4D + "b4_tiras.py": "b1a2cb8f93affcc3713eadd4ff860a2d74756b91752011cd6b16da44def2dd15",
    B4D + "RESULTADO_B4.json": "640db376441656cfb72e16aa00e0e99e2fc2eb6047905dfe031f5d994634b5f7",
    RES_T21: "9769c76eae7d54a5377b0397d3ad0e0ee6f68207cef8f3c70ec4971fe9f6df69",
}
M_STRIP = 6
N_LOCS = 527
NS_PLAN = 91
NCLS_PLAN = 409
RES_MAX = 300000
V4_HISTORIAS = 20000


def check_deps():
    import t2s_reg as tr
    for p, s in SHA.items():
        if tr.sha_file(p) != s:
            raise ValueError("SHA-256 distinto del registrado: " + p)   # queda registrado en la preparacion


class Fisico:
    def __init__(self, spec):
        check_deps()
        sys.path.insert(0, os.path.abspath(T21))
        try:
            import t2_mc as mcm  # noqa: E402  (comprueba t2_calibracion, T2-0, b4_tiras y B4)
        except SystemExit as ex:                      # t2_mc sale con SystemExit ante un SHA-256 distinto
            raise ValueError("t2_mc.py no se pudo importar: %s" % ex)
        self.mcm, self.t2, self.b4 = mcm, mcm.t2, mcm.b4
        t2, b4 = self.t2, self.b4
        ce, rect, _ = b4.load_all(with_b1=False)
        inst = b4.make_inst(ce, rect, "tira", M_STRIP)
        self.ce, self.inst = ce, inst
        self.dec, self.RES = ce.make_decoder(inst)
        self.layer, self.L = t2.schedule(inst)
        self.locs = t2.locations(inst, self.layer, self.L)
        CL, CK, LTK = t2.site_classes(inst, self.layer, self.L, self.locs, False)
        self.CL, self.CK, self.LTK = CL, CK, LTK
        self.NS = inst["NS"]
        self.ncls = [len(x) for x in CL]
        self.resp = [[(0, 0)] + [tuple(es) for es in CL[s]] for s in range(self.NS)]
        self.zero = [tuple(0 for _ in CK[s][0]) for s in range(self.NS)]
        self.ck_index = [{k: c + 1 for c, k in enumerate(CK[s])} for s in range(self.NS)]
        self.r0 = F(spec["r0"])
        self._anchor_tables()

    # -------------------------------------------------------------- leyes (V3)
    def laws(self, r):
        out = []
        for s in range(self.NS):
            dist = {k: v for k, v in self.t2.law_of(self.LTK[s], self.L, r, r).items() if v != 0}
            want = {self.zero[s]} | set(self.CK[s])
            if set(dist) != want:
                raise ValueError("V3: soporte del sitio %d distinto de {0} U C_s en r = %s" % (s, r))
            out.append([dist[self.zero[s]]] + [dist[k] for k in self.CK[s]])
        return out

    # -------------------------------------------------------------- evento (una llamada al decoder)
    def oracle(self, e, syn, cls):
        if len(self.RES) > RES_MAX:
            self.RES.clear()
        cor, _ = self.dec(syn)
        r = e ^ cor
        return self.ce.lin(self.inst["cols"], r) == 0 and not self.ce.inside(r, self.inst["bz"])

    # -------------------------------------------------------------- ancla (t2_mc.w_block, a sitios)
    def _site_of_loc(self):
        inst, layer, t2 = self.inst, self.layer, self.t2
        N = inst["N"]
        site, keyfun = [], []
        for loc in self.locs:
            if loc[0] in ("dwait", "cnot"):
                v = loc[1]
                chk = sorted(self.b4.bits(inst["cols"][v]), key=lambda c, v=v: layer[(c, v)])
                site.append(v)
                keyfun.append(lambda e, s, chk=chk: (1 if e else 0,) + tuple(s >> c & 1 for c in chk))
            else:
                c0 = loc[1]
                site.append(N + c0)
                keyfun.append(lambda e, s, c0=c0: (1 if s >> c0 & 1 else 0,))
        return site, keyfun

    def _anchor_tables(self):
        mcm, t2 = self.mcm, self.t2
        r0 = self.r0
        laws = mcm.location_laws(self.locs, r0, r0)
        self.loc_laws = laws
        groups = {}
        for i, law in enumerate(laws):
            groups.setdefault(sum(w for _, w in law), []).append(i)
        self.groups = [(float(q), math.log1p(-float(q)), groups[q], q) for q in sorted(groups)]
        self.cats = []
        for law in laws:
            tot = sum(w for _, w in law)
            acc, cum = 0, []
            for c, w in law:
                acc += w / tot
                cum.append((float(acc), c))
            cum[-1] = (1.0, cum[-1][1])
            self.cats.append(cum)
        site, keyfun = self._site_of_loc()
        self.loc_site = site
        self.loc_key = []
        for i, law in enumerate(laws):
            d = {}
            for c, _ in law:
                e, s = t2.resp_formula(self.inst, self.layer, self.locs[i], c)
                d[c] = keyfun[i](e, s)
            self.loc_key.append(d)

    def history_faults(self, rng, Ba, tick):
        """Como t2_mc.w_block: (h, [(indice de ubicacion, categoria)]) en orden creciente de h (0..Ba-1)."""
        rnd, log = rng.random, math.log
        hits, it = {}, 0
        for _, lq, idxs, _ in self.groups:
            ng = len(idxs)
            T = Ba * ng
            pos = -1
            while True:
                pos += 1 + int(log(1.0 - rnd()) / lq)
                if pos >= T:
                    break
                h, r = divmod(pos, ng)
                lst = hits.get(h)
                if lst is None:
                    hits[h] = [idxs[r]]
                else:
                    lst.append(idxs[r])
                it += 1
                if not it % 65536:
                    tick()
        for h in sorted(hits):
            faults = []
            for i in sorted(hits[h]):
                cum = self.cats[i]
                if len(cum) == 1:
                    faults.append((i, cum[0][1]))
                else:
                    u = rnd()
                    for th, c in cum:
                        if u < th:
                            faults.append((i, c))
                            break
            yield h, faults

    def to_sites(self, faults):
        keys = {}
        for i, c in faults:
            b = self.loc_site[i]
            k = self.loc_key[i][c]
            old = keys.get(b)
            keys[b] = k if old is None else tuple(x ^ y for x, y in zip(old, k))
        pairs = []
        for b in sorted(keys):
            k = keys[b]
            if k != self.zero[b]:
                if k not in self.ck_index[b]:
                    raise ValueError("V4: clave de sitio %r fuera de C_%d" % (k, b))
                pairs.append((b, self.ck_index[b][k]))
        return pairs

    def anchor_sample(self, rng, Ba, tick):
        for h, faults in self.history_faults(rng, Ba, tick):
            pairs = self.to_sites(faults)
            if pairs:
                yield h + 1, pairs

    # -------------------------------------------------------------- controles previos (sin decoder)
    def prechecks(self, rates_check, master, check):
        mcm, t2, b4, inst = self.mcm, self.t2, self.b4, self.inst
        RB4 = json.load(open(B4D + "RESULTADO_B4.json", encoding="utf-8"))["instancias"]
        check("T2-2 TIRA_6x6: matrices iguales a las de B4", b4.mats_sha(inst) == RB4["TIRA_6x6"]["matrices_sha256"])
        check("T2-2: %d ubicaciones (plan: %d)" % (len(self.locs), N_LOCS), len(self.locs) == N_LOCS)
        check("T2-2: NS = %d y %d clases no nulas (plan: %d y %d)" % (self.NS, sum(self.ncls), NS_PLAN, NCLS_PLAN),
              self.NS == NS_PLAN and sum(self.ncls) == NCLS_PLAN)
        fd, fc = t2.fibres(inst, self.locs)
        for r in [F(x) for x in rates_check]:
            laws = mcm.location_laws(self.locs, r, r)
            lbl = {loc: law for loc, law in zip(self.locs, laws)}
            w = 2 * r / 3
            m1, mc, mr = mcm.marginal_1q(r), mcm.marginal_cnot(r), mcm.marginal_read(r, r)
            check("T2-2 r = %s: marginales X = formas cerradas" % r,
                  m1[1] == w and all(mc[c] == 4 * r / 15 for c in ((1, 0), (0, 1), (1, 1)))
                  and mr[1] == w * (1 - r) + (1 - w) * r)
            bad = 0
            for b in range(self.NS):
                is_data = b < inst["N"]
                idx = b if is_data else b - inst["N"]
                mine = mcm.site_law_from_locations(inst, self.layer, self.L, idx, is_data, fd, fc, lbl)
                ref = {k: v for k, v in t2.law_of(self.LTK[b], self.L, r, r).items() if v != 0}
                bad += int(mine != ref)
            check("T2-2 r = %s: leyes de sitio desde las marginales = law_of en los %d sitios" % (r, self.NS),
                  bad == 0, bad)
        bad = 0
        for i, law in enumerate(self.loc_laws):
            for c, _ in law:
                bad += int(t2.propagate(inst, self.layer, self.L, [(self.locs[i], c)])
                           != t2.resp_formula(inst, self.layer, self.locs[i], c))
        check("T2-2: propagate de cada fallo unico = resp_formula", bad == 0, bad)
        rng = tmo.stream(tmo.label_controles(master))
        bad = n = 0
        for h, faults in self.history_faults(rng, V4_HISTORIAS, lambda: None):
            e = y = 0
            for s, c in self.to_sites(faults):
                ee, yy = self.resp[s][c]
                e ^= ee
                y ^= yy
            ref = t2.propagate(inst, self.layer, self.L, [(self.locs[i], c) for i, c in faults])
            bad += int((e, y) != ref)
            n += 1
        check("T2-2 V4: propagate = XOR de las respuestas de sitio en %d historias del ancla" % n, bad == 0 and n > 0,
              {"historias": n, "distintas": bad})
        R0 = json.load(open(T20 + "RESULTADO_T2_0.json", encoding="utf-8"))["instancias"]["TIRA_6x6"]
        ex = R0["exhaustivo"]
        check("T2-2 V5 (T2-0): exhaustivo con 0 silenciosas con 1 y 2 sitios, 78 con 3; f_s_exhaustivo = 3",
              ex["1"]["silenciosas"] == 0 and ex["2"]["silenciosas"] == 0 and ex["3"]["silenciosas"] == 78
              and ex["1"]["historias"] == 409 and ex["2"]["historias"] == 82371 and R0["f_s_exhaustivo"] == 3)
        R1 = json.load(open(RES_T21, encoding="utf-8"))["instancias"]["TIRA_6x6"]
        rows = [c["por_k"]["2"] for c in R1["celdas"].values()]
        check("T2-2 V5 (T2-1): 8 filas k = 2 con H_comp COMPLETO, L_ing = 0 y L_comp = [0, 0]",
              len(rows) == 8 and all(r["H_comp"] == "COMPLETO" and r["L_ing"] == "0"
                                     and r["L_comp"] == ["0", "0"] for r in rows))


def make(spec):
    return Fisico(spec)
