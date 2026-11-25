#!/usr/bin/env python3
"""T2-1 (BORRADOR_PLAN_T2_1.txt, N1): rival ingenuo con clases (eps_ing y L_ing por el lema P3C), en fracciones
exactas.

Codigo preparado para revision (version 2, tras la revision de N1 v1: R1 y R2). No se ha ejecutado sobre los
catalogos de circuito. Solo se ejecutan su autoprueba (--autoprueba: matematica y lector con grupos sinteticos en
el formato del motor) y, como prueba de desarrollo sin decoder ni DP, el control fenomenologico (--control-fen).
- R1: en la lectura, las silenciosas de 3 sitios se agrupan por (seccion, indice local), no por el indice local
  solo, que se repite entre secciones de un mismo grupo (grupos de 2 y 3 anclas en m = 24 y 48).
- R2: la espera de los grupos (plazo de POLL s) y los bucles largos consultan los techos del conductor (GUARD) y
  el plazo; el avance de la lectura queda en PROGRESS como evidencia si hay una parada.

Contrato (N1 del borrador):
- eps_ing(k) = Pr[|F| >= k + 1], con |F| = numero de sitios no nulos. Es la cola exacta de una suma de Bernoulli
  independientes, una por sitio, con Pr[no nulo] = 1 - z_b, donde z_b es la probabilidad de la clase nula de la
  ley de su fibra (t2_calibracion.law_of). En el control fenomenologico es p para los datos y q para los checks.
  Se calcula como producto de binomiales por tipo de ley, truncado en el grado k, en fracciones.
- L_ing(k) por el lema P3C, con los catalogos de T2-1M (trabajo/, sin cambios, SHA-256 comprobado al leer cada
  grupo) y sin llamar al decoder:
  (a) decoraciones validas silenciosas con <= k sitios: n_minus por (Y, firma);
  (b) con k = f_s + 1 = 4: una decoracion silenciosa g de 3 sitios (n_minus) y un singlete x que pasa con clase
      trivial (n_plus), con Y(x) disjunto de N_K[Y(g)] (campo NB de la seccion; YR y NB en el orden del motor).
  L_ing(k) = B [ suma_(a) n_minus W(firma) + suma_(b) n_minus(g) W(g) n_plus(x) W(x) ]. W es el producto de los
  odds de los simbolos (t2_motor.odds_circ u odds_fen) y B la probabilidad de que todas las bases sean nulas.
  Con k <= 2, L_ing = 0. Es aplicable si no hay ninguna silenciosa con < 3 sitios (se comprueba).
- H-comp: L_ing es exacto solo si la completitud del catalogo esta certificada en esa instancia y ese corte.
  Procedencia previa: T2-0 y C-T20 de T2-1M (m = 6 y 12 hasta k = 4; m = 24 hasta k = 3). Para m = 24 con k = 4
  y para m = 48, decide la regla de suma de T2-1 (N3). COMPLETO en el corte k certifica todos los cortes <= k,
  porque H-comp en k pide todas las decoraciones validas con <= k sitios. Sin certificado, la suma del catalogo,
  L_cat, es solo cota inferior de L_ing y se etiqueta asi.
- Coste del rival: construccion de las leyes de fibra (calendario, ubicaciones, respuestas y XOR por fibra con
  t2_calibracion.site_classes; leyes por celda), en segundos y operaciones contadas; la cola exacta y P3C, en
  segundos; y la truncacion ingenua de la literatura, suma_(f<=k) e_f(n) decodificaciones.

Controles (detienen):
- SHA-256 de las dependencias.
- Lectura completa del catalogo:
  - validas por f, polimeros y entradas por k, numero de firmas, grupos y anclas por grupo iguales a
    RESULTADO_T2_1M.json;
  - todos los simbolos con odds definidos;
  - todas las decoraciones de un sitio validas.
- P3C aplicable: ninguna decoracion silenciosa con < 3 sitios.
- Fenomenologico (m = 6, 12, 24 y 48):
  - silenciosas conexas por f y (na, nb), singletes, pares compatibles y silenciosas totales = P3 de B4;
  - eps_ing y L_ing con k = 1..4 identicos (fracciones) a RESULTADO_B4.json en las 12 celdas;
  - truncacion ingenua = llamadas_ingenua_convencional de B4.
- S2: eps_ing con k = 1..3 identico a RESULTADO_CONTROL.json (D1-CE) y a RESULTADO_B4.json en las 12 celdas.
- Circuito (los usa N3):
  - silenciosas por f = T2-0 (m = 6 y 12 con f <= 4; m = 24 con f <= 3);
  - truncacion ingenua = T2-0;
  - L_ing(3) = termino dominante de T2-0 (misma fraccion: el float coincide) en las 8 celdas;
  - B = Pr[|F| = 0].

Uso (desde D:\\paper ia cuantica\\qec_revision_matematica\\05_COMPROBACIONES):
  python -I -B T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/t2_rival.py --autoprueba
  python -I -B T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/t2_rival.py --control-fen --output <nuevo.json>
"""
import argparse
import bisect
import hashlib
import json
import math
import multiprocessing as mp
import os
import pickle
import random
import sys
import tempfile
import time
import zlib
from fractions import Fraction as F
from itertools import combinations, compress

sys.set_int_max_str_digits(0)
T0 = time.perf_counter()
T0_WALL = time.time()
T21M = "T2_1M_MOTOR_ANCLAS_PUNTO_FIJO_20261113/"
T20 = "T2_0_CALIBRACION_CIRCUITO_TIRAS_20261113/"
B4D = "B4_TIRAS_ESCALADO_N_20261113/"
RCE = "D1_CONTROL_EXACTO_Y_METODOLOGIA_20261113/RESULTADO_CONTROL.json"
WDIR = T21M + "trabajo/"
SHA = {
    T21M + "t2_motor.py": "c0d297805a55ecfe514fd623c21ac5a1f272dd2b26a2f52a6f10a87205e6cf3a",
    T21M + "RESULTADO_T2_1M.json": "573e4e0f39b853a7f40acc4b4d2d7dce2155aa3ff521e5ae165e2c62719b6358",
    WDIR + "MANIFIESTO.json": "4671bab0a5eec923b1fd7b9fdd7af499bd64dddf9efce4aedd372495d8f35730",
    T20 + "RESULTADO_T2_0.json": "fb56c4324034549476813bd1a0d2bacd6f46c27d552fbf71ba2d0813d62df38e",
    B4D + "RESULTADO_B4.json": "640db376441656cfb72e16aa00e0e99e2fc2eb6047905dfe031f5d994634b5f7",
    RCE: "334d6485968c914c299bbc7e5362cf49f4b975b740d47d78bccfaec47255024d",
}


def sha_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 24), b""):
            h.update(blk)
    return h.hexdigest()


for _p, _s in SHA.items():
    if sha_file(_p) != _s:
        raise SystemExit("SHA-256 distinto del registrado: " + _p)
sys.path.insert(0, os.path.abspath(T21M))
import t2_motor as mot  # noqa: E402  (tras comprobar su SHA-256; el motor comprueba los suyos)

t2 = mot.t2
b4 = mot.b4
bits = b4.bits
check = b4.check

KMAX = 4
F_S = 3                      # orden silencioso con circuito (T2-0); P3C usa k <= F_S + 1
PROCS = 14
STRIPS = (6, 12, 24, 48)
C6, C12, C18, C24 = 1 << 6, 1 << 12, 1 << 18, 1 << 24
# H-comp previo a T2-1: (corte maximo certificado, fuente). Se comprueba contra los registros en prov_checks.
HCOMP_PREVIO = {
    "circ_TIRA_6x6": (4, "T2-0 (regla de suma con k = 4, COMPLETO), reproducida por C-T20 de T2-1M"),
    "circ_TIRA_12x6": (4, "T2-0 (regla de suma con k = 4, COMPLETO), reproducida por C-T20 de T2-1M"),
    "circ_TIRA_24x6": (3, "T2-0 (regla de suma con k = 3, COMPLETO), reproducida por C-T20 de T2-1M"),
    "circ_TIRA_48x6": (0, "sin regla de suma anterior a T2-1"),
}
T20_KMAX = {6: 4, 12: 4, 24: 3}
# R2 (revision de N1 v1): techos durante la lectura y los bucles largos de N1. GUARD[0] lo fija el conductor (N3:
# plazo, memoria del principal, memoria total marcada por t2_motor.Monitor y disco); sin conductor solo se mira el
# plazo. La parada es cooperativa: como mucho POLL segundos en espera, o un grupo, o GUARD_EVERY iteraciones.
POLL = 5.0
GUARD = [None]
GUARD_EVERY = 1 << 12
PROGRESS = {}


def _guard(where, deadline=None):
    if deadline is not None and time.time() > deadline:
        raise b4.Stop("techo de tiempo de la unidad (%s)" % where)
    if GUARD[0] is not None:
        GUARD[0](where)


def nsym(s):
    """Numero de simbolos de 6 bits de una firma (= sitios de la decoracion)."""
    return (s.bit_length() + 5) // 6


def show(x):
    """Vista de 25 cifras (cota inferior y superior coinciden en una fraccion exacta)."""
    return mot.dec(x, False)


def load_refs():
    j = lambda p: json.load(open(p, encoding="utf-8"))  # noqa: E731
    return {"B4": j(B4D + "RESULTADO_B4.json")["instancias"], "T20": j(T20 + "RESULTADO_T2_0.json")["instancias"],
            "T21M": j(T21M + "RESULTADO_T2_1M.json"), "CE": j(RCE)["instancias"],
            "MAN": j(WDIR + "MANIFIESTO.json")}


def prov_checks(refs):
    """La procedencia de HCOMP_PREVIO es la de los registros: regla de suma COMPLETO en T2-0 y en C-T20 de T2-1M."""
    RT, ctl = refs["T20"], refs["T21M"]["controles"]["C-T20"]
    ok = refs["T21M"]["status"] == "PASS" and refs["T21M"]["decision_D_m1"] == "VALIDADO"
    for m, kk in T20_KMAX.items():
        rs = RT["TIRA_%dx6" % m]["regla_de_suma"]
        c = ctl.get("circ_TIRA_%dx6_suma_k%d" % (m, kk), {})
        ok = ok and rs["k"] == kk and rs["veredicto"] == "COMPLETO" and c.get("veredicto") == "COMPLETO"
        ok = ok and HCOMP_PREVIO["circ_TIRA_%dx6" % m][0] == kk
    check("N1 procedencia de H-comp: regla de suma COMPLETO en T2-0 y en C-T20 de T2-1M (m = 6 y 12 con k = 4; "
          "m = 24 con k = 3); T2-1M PASS y D-m1 VALIDADO", ok)


def hcomp(tag, k, sum_rules):
    """Veredicto de H-comp en (instancia, corte k). sum_rules: {corte: veredicto} de la regla de suma de T2-1."""
    prev, src = HCOMP_PREVIO.get(tag, (0, "sin regla de suma"))
    if k <= prev:
        return "COMPLETO", src
    for kk in sorted(sum_rules):
        if kk >= k and sum_rules[kk] == "COMPLETO":
            return "COMPLETO", "regla de suma de T2-1 con k = %d" % kk
    if k in sum_rules:
        return sum_rules[k], "regla de suma de T2-1 con k = %d" % k
    return "SIN_REGLA", "sin regla de suma que cubra este corte"


# ------------------------------------------------------------------ instancia, clases y leyes (sin decoder)
def strip_context(ce, rect, m, restricted, RB4):
    name = "TIRA_%dx6" % m
    tag = ("fen_" if restricted else "circ_") + name
    inst = b4.make_inst(ce, rect, "tira", m)
    check("N1 %s: matrices = B4" % tag, b4.mats_sha(inst) == RB4[name]["matrices_sha256"], verbose=False)
    t = time.perf_counter()
    layer, L = t2.schedule(inst)
    locs = t2.locations(inst, layer, L)
    CL, CK, LTK = t2.site_classes(inst, layer, L, locs, restricted)
    t_cls = time.perf_counter() - t
    lts = sorted(set(LTK))
    LT = [lts.index(x) for x in LTK]
    NS = inst["NS"]
    check("N1 %s: simbolos (tipo de ley, clase) caben en 6 bits" % tag,
          len(lts) <= 7 and max(len(c) for c in CL) <= 8, verbose=False)
    CKrep = [CK[LT.index(i)] for i in range(len(lts))]
    check("N1 %s: las clases de cada sitio son las de su tipo de ley" % tag,
          all(CK[b] == CKrep[LT[b]] for b in range(NS)), verbose=False)
    check("N1 %s: k_log = 1" % tag, len(ce.logical_reps(inst)) == 1, verbose=False)
    _, _, gsize, ng = mot.layout(inst)
    nclass = [len(c) for c in CL]
    cx = {"lts": lts, "cnt": [LT.count(i) for i in range(len(lts))], "CKrep": CKrep, "L": L, "N": inst["N"],
          "M": inst["M"], "nclass": nclass}
    cost = {"segundos_clases": round(t_cls, 6), "ubicaciones": len(locs),
            "respuestas_por_ubicacion": sum(len(t2.VISIBLE[loc[0]]) for loc in locs), "fibras": NS,
            "clases_no_nulas": sum(nclass), "tipos_de_ley": len(lts),
            "nota": "respuestas y XOR por fibra de t2_calibracion.site_classes; segundos de reloj"}
    return {"tag": tag, "name": name, "m": m, "restricted": restricted, "NS": NS, "gsize": gsize, "ngroups": ng,
            "cx": cx, "coste_clases": cost}


def law_cost(cx, cells):
    """Coste de las leyes de fibra por celda (t2_calibracion.law_of), con las convoluciones contadas."""
    t = time.perf_counter()
    conv = 0
    for rc, rm in cells:
        for lt in cx["lts"]:
            t2.law_of(lt, cx["L"], rc, rm)
            conv += (cx["L"] + 2) if lt[0] == "d" else (lt[1] + 2)
    return {"segundos": round(time.perf_counter() - t, 6), "convoluciones": conv, "celdas": len(cells),
            "tipos_de_ley": len(cx["lts"])}


def groups_circ(cx, rc, rm):
    g = {}
    for t, lt in enumerate(cx["lts"]):
        law = t2.law_of(lt, cx["L"], rc, rm)
        z = law[tuple([0] * len(next(iter(law))))]
        g[1 - z] = g.get(1 - z, 0) + cx["cnt"][t]
    return sorted(((n, pi) for pi, n in g.items()), key=lambda x: x[1])


def groups_fen(cx, p, q):
    nd = sum(c for t, c in enumerate(cx["cnt"]) if cx["lts"][t][0] == "d")
    nc = sum(cx["cnt"]) - nd
    check("N1: tipos de dato y de check = N y M", nd == cx["N"] and nc == cx["M"], verbose=False)
    return [(cx["N"], p), (cx["M"], q)]


def pb_terms(groups, K):
    """P(|F| = j), j = 0..K, para una suma de Bernoulli independientes: groups = [(cantidad, Pr[no nulo])]."""
    poly = [F(1)] + [F(0)] * K
    for n, pi in groups:
        pi = F(pi)
        if not 0 <= pi <= 1:
            raise ValueError("probabilidad fuera de [0, 1]")
        if pi == 0 or n == 0:
            continue
        q = 1 - pi
        bt = [math.comb(n, a) * pi ** a * q ** (n - a) for a in range(min(n, K) + 1)]
        new = [F(0)] * (K + 1)
        for i, x in enumerate(poly):
            if x:
                for a, y in enumerate(bt):
                    if i + a > K:
                        break
                    new[i + a] += x * y
        poly = new
    return poly


def eps_from_terms(terms, k):
    return 1 - sum(terms[:k + 1], F(0))


def esym_counts(nclass, K):
    """Truncacion ingenua de la literatura: suma_(f<=k) e_f(n) decoraciones (una decodificacion cada una)."""
    e = t2.esym(nclass, K)
    return {str(k): sum(e[1:k + 1]) for k in range(1, K + 1)}


# ------------------------------------------------------------------ lectura del catalogo (procesos)
def w_scan(task):
    try:
        return _scan(task)
    except Exception as ex:
        return {"error": "%s: %s" % (type(ex).__name__, ex), "g": task[0]}


def _scan(task):
    """Un grupo: comprueba su SHA-256 (read_group del motor), recorre sus secciones y devuelve recuentos de control,
    las silenciosas por firma, los singletes y las silenciosas de 3 sitios con su soporte y su NB."""
    g, path, sha, deadline = task
    blob, secs = mot.read_group(path, sha)
    mv = memoryview(blob)
    va = [0] * (KMAX + 1)        # suma de n_all con firma < 2^(6k) (acumulado)
    ent = [0] * (KMAX + 1)
    pol = [0] * (KMAX + 1)
    tot = [0, 0, 0]
    silent = {}
    single = []
    sil3 = {}
    sigs = set()
    over = 0
    npoly = 0
    for a in sorted(secs):
        if time.time() > deadline:
            raise RuntimeError("techo de tiempo de la unidad (lectura del catalogo)")
        p, n = secs[a]
        d = pickle.loads(zlib.decompress(mv[p:p + n]))
        if d["a"] != a:
            raise RuntimeError("seccion con ancla distinta")
        YR, NB, off, sig, nA, nP, nM = (d[x] for x in ("YR", "NB", "off", "sig", "nA", "nP", "nM"))
        ns = len(sig)
        if len(off) != len(YR) + 1 or off[0] != 0 or off[-1] != ns or not len(nA) == len(nP) == len(nM) == ns \
                or len(NB) != len(YR):
            raise RuntimeError("seccion mal formada en el ancla %d" % a)
        npoly += len(YR)
        sigs.update(sig)
        if ns and max(sig) >= C24:
            over += 1
        tot[0] += sum(nA)
        tot[1] += sum(nP)
        tot[2] += sum(nM)
        firsts = list(map(sig.__getitem__, off[:-1]))
        m1 = None
        for k, cap in ((1, C6), (2, C12), (3, C18)):
            msk = list(map(cap.__gt__, sig))
            if k == 1:
                m1 = msk
            ent[k] += sum(msk)
            va[k] += sum(compress(nA, msk))
            pol[k] += sum(map(cap.__gt__, firsts))
        ent[4] += ns
        va[4] += sum(nA)
        pol[4] += len(YR)
        for j in compress(range(ns), nM):
            s, v = sig[j], nM[j]
            silent[s] = silent.get(s, 0) + v
            if C12 <= s < C18:
                i = bisect.bisect_right(off, j) - 1
                # R1 (revision de N1 v1): i es local a la seccion a; la clave es (a, i), no i, para no mezclar
                # polimeros de secciones distintas del mismo grupo (YR y NB del primero con firmas del segundo)
                e = sil3.get((a, i))
                if e is None:
                    sil3[(a, i)] = e = [YR[i], NB[i], []]
                e[2].append((s, v))
        for j in compress(range(ns), m1):
            i = bisect.bisect_right(off, j) - 1
            single.append((YR[i], sig[j], nP[j], nA[j], nM[j]))
        del d, YR, NB, off, sig, nA, nP, nM, firsts, m1
    mv.release()
    return {"g": g, "va": va, "ent": ent, "pol": pol, "tot": tot, "silent": silent, "single": single,
            "sil3": list(sil3.values()), "sigs": sorted(sigs), "sobre_k4": over, "polimeros": npoly,
            "secciones": len(secs), "bytes": len(blob)}


def scan_catalogue(pool, sc, MAN, R1M, deadline, wdir=WDIR, worker=w_scan):
    """Lee el catalogo de T2-1M de la instancia (todos los grupos, con su SHA-256) y comprueba la lectura completa.
    R2: espera los resultados con plazo de POLL segundos y consulta los techos (_guard) en cada espera y tras cada
    grupo recibido; el avance queda en PROGRESS[tag] como evidencia si hay una parada."""
    tag = sc["tag"]
    man = MAN[tag]
    t = time.perf_counter()
    tasks = [(int(g), os.path.join(wdir, tag, e["archivo"]), e["sha256"], deadline) for g, e in man.items()]
    va = [0] * (KMAX + 1)
    ent = [0] * (KMAX + 1)
    pol = [0] * (KMAX + 1)
    tot = [0, 0, 0]
    silent, single, sil3, sigs = {}, [], [], set()
    over = nbytes = nsec = 0
    seen = set()
    prog = PROGRESS[tag] = {"estado": "leyendo", "grupos": len(tasks), "grupos_leidos": 0, "grupos_ids": [],
                            "secciones": 0, "bytes": 0, "segundos": 0.0, "esperas": 0}
    it = pool.imap_unordered(worker, tasks, chunksize=1)
    while True:
        try:
            r = it.next(timeout=POLL)
        except mp.TimeoutError:
            prog["esperas"] += 1
            prog["segundos"] = round(time.perf_counter() - t, 3)
            _guard("N1, catalogo de %s, en espera" % tag, deadline)
            continue
        except StopIteration:
            break
        if "error" in r:
            prog["estado"] = "detenido"
            if "techo de tiempo" in r["error"]:
                raise b4.Stop("N1 %s: %s" % (tag, r["error"]))
            raise b4.Stop("N1 %s grupo %d: %s" % (tag, r["g"], r["error"]))
        seen.add(r["g"])
        prog.update(grupos_leidos=len(seen), grupos_ids=sorted(seen), secciones=prog["secciones"] + r["secciones"],
                    bytes=prog["bytes"] + r["bytes"], segundos=round(time.perf_counter() - t, 3))
        _guard("N1, catalogo de %s" % tag, deadline)
        for i in range(KMAX + 1):
            va[i] += r["va"][i]
            ent[i] += r["ent"][i]
            pol[i] += r["pol"][i]
        for i in range(3):
            tot[i] += r["tot"][i]
        for s, v in r["silent"].items():
            silent[s] = silent.get(s, 0) + v
        single += r["single"]
        sil3 += r["sil3"]
        sigs.update(r["sigs"])
        over += r["sobre_k4"]
        nbytes += r["bytes"]
        nsec += r["secciones"]
    prog["estado"] = "leido"
    secs = round(time.perf_counter() - t, 3)
    rc = R1M["instancias"][tag]
    cat = rc["catalogo"]
    validas = [va[f] - va[f - 1] for f in range(1, KMAX + 1)]
    exp = [cat["por_f"][str(f)]["validas"] for f in range(1, KMAX + 1)]
    check("N1 %s: grupos leidos = manifiesto (%d), sitios, grupos y anclas por grupo = T2-1M" % (tag, len(man)),
          seen == {int(g) for g in man} and sc["NS"] == rc["sitios"] and sc["ngroups"] == rc["grupos"]
          and sc["gsize"] == rc["anclas_por_grupo"], verbose=False)
    check("N1 %s: validas por f = T2-1M %s" % (tag, exp), validas == exp, validas)
    check("N1 %s: polimeros y entradas por k y numero de firmas = T2-1M" % tag,
          [pol[k] for k in range(1, KMAX + 1)] == [cat["polimeros_por_k"][str(k)] for k in range(1, KMAX + 1)]
          and [ent[k] for k in range(1, KMAX + 1)] == [cat["entradas_por_k"][str(k)] for k in range(1, KMAX + 1)]
          and len(sigs) == cat["firmas"] and nbytes == cat["bytes"] and over == 0,
          {"polimeros": pol[1:], "entradas": ent[1:], "firmas": len(sigs), "bytes": nbytes})
    cx = sc["cx"]
    okodds = all(all((sy >> 3) < len(cx["lts"]) and (sy & 7) < len(cx["CKrep"][sy >> 3]) for sy in mot.sig_syms(s))
                 for s in sigs)
    check("N1 %s: todos los simbolos de las firmas tienen odds definidos" % tag, okodds, verbose=False)
    check("N1 %s: todas las decoraciones de un sitio son validas (%d = suma de clases)" % (tag, validas[0]),
          validas[0] == sum(cx["nclass"]), validas[0])
    check("N1 %s: n_plus + n_minus <= n_all en el total del catalogo" % tag, tot[1] + tot[2] <= tot[0], tot)
    sil_f = {f: 0 for f in range(1, KMAX + 1)}
    for s, v in silent.items():
        sil_f[nsym(s)] += v
    fs = min((f for f, v in sil_f.items() if v), default=None)
    check("N1 %s: P3C aplicable, ninguna decoracion silenciosa con < %d sitios (f_s = %s)" % (tag, F_S, fs),
          fs is None or fs >= F_S, sil_f)
    t1 = time.perf_counter()
    pairs, npairs, nben, wben = p3_pairs(single, sil3)
    t_pairs = time.perf_counter() - t1
    return {"tag": tag, "silenciosas": silent, "silenciosas_por_f": {str(f): v for f, v in sil_f.items()},
            "f_s": fs, "pares": pairs, "pares_compatibles": npairs,
            "singletes_que_pasan_con_clase_trivial": sum(x[2] for x in single),
            "singletes_que_no_pasan": sum(x[3] - x[2] - x[4] for x in single),
            "singletes_silenciosos": sum(x[4] for x in single), "entradas_de_singletes_benignos": nben,
            "polimeros_silenciosos_de_3_sitios": len(sil3), "firmas": sorted(sigs), "validas_por_f": validas,
            "secciones": nsec, "bytes": nbytes,
            "coste": {"segundos_lectura": secs, "segundos_pares": round(t_pairs, 3)}}


def p3_pairs(single, sil3):
    """P3C (b): pares (g silenciosa de 3 sitios, singlete x que pasa con clase trivial) con YR(x) & NB(g) = 0.
    Devuelve {(firma de g, firma de x): numero de pares} y el total. Indice por bit para los singletes que cortan NB."""
    ben = [(Y, s, nP) for (Y, s, nP, nA, nM) in single if nP]
    TOT, idx = {}, {}
    for t, (Y, s, n) in enumerate(ben):
        if t % GUARD_EVERY == 0:
            _guard("N1, indice de singletes de P3C")
        TOT[s] = TOT.get(s, 0) + n
        for b in bits(Y):
            idx.setdefault(b, []).append(t)
    pairs = {}
    total = 0
    for u, (Y, NB, ents) in enumerate(sil3):
        if u % GUARD_EVERY == 0:
            _guard("N1, pares de P3C")
        C = dict(TOT)
        hit = set()
        for b in bits(NB):
            hit.update(idx.get(b, ()))
        for t in hit:
            C[ben[t][1]] -= ben[t][2]
        for sg, nm in ents:
            for sx, c in C.items():
                if c:
                    pairs[(sg, sx)] = pairs.get((sg, sx), 0) + nm * c
                    total += nm * c
    return pairs, total, len(ben), sum(TOT.values())


# ------------------------------------------------------------------ valores del rival por celda
def rival_values(scan, od, B, terms, kmax=KMAX):
    """eps_ing(k) y L_cat(k) exactos para k = 1..kmax. L_cat = L_ing si H-comp esta certificado (lo decide hcomp)."""
    W = {}

    def w(s):
        x = W.get(s)
        if x is None:
            x = F(1)
            for sym in mot.sig_syms(s):
                x *= od[sym]
            W[s] = x
        return x

    La = [F(0)] * (KMAX + 1)
    for u, (s, c) in enumerate(scan["silenciosas"].items()):
        if u % GUARD_EVERY == 0:
            _guard("N1, silenciosas del rival")
        La[nsym(s)] += c * w(s)
    Lb = F(0)
    for u, ((sg, sx), c) in enumerate(scan["pares"].items()):
        if u % GUARD_EVERY == 0:
            _guard("N1, pares del rival")
        Lb += c * w(sg) * w(sx)
    out = {}
    for k in range(1, kmax + 1):
        a = sum(La[1:k + 1], F(0))
        b = Lb if k == F_S + 1 else F(0)
        out[k] = {"eps_ing": eps_from_terms(terms, k), "L_cat": B * (a + b), "L_a": B * a, "L_b": B * b}
    return out


def nanb(s, lts):
    na = sum(1 for sy in mot.sig_syms(s) if lts[sy >> 3][0] == "d")
    return na, nsym(s) - na


# ------------------------------------------------------------------ controles de N1
def fen_control(pool, ce, rect, refs, deadline, grid12):
    """Clases fenomenologicas en m = 6, 12, 24 y 48: identidad con B4 (fracciones)."""
    RB4 = refs["B4"]
    out = {}
    for m in STRIPS:
        _guard("N1, control fenomenologico m = %d" % m, deadline)
        sc = strip_context(ce, rect, m, True, RB4)
        tag, name = sc["tag"], sc["name"]
        scan = scan_catalogue(pool, sc, refs["MAN"], refs["T21M"], deadline)
        cat = RB4[name]["catalogo"]
        P3 = cat["P3"]
        lts = sc["cx"]["lts"]
        conn = {}
        for s, v in scan["silenciosas"].items():
            f = str(nsym(s))
            kk = "%d,%d" % nanb(s, lts)
            conn.setdefault(f, {})
            conn[f][kk] = conn[f].get(kk, 0) + v
        exp = {f: d for f, d in cat["silenciosas_conexas_por_f"].items() if d}
        check("N1 %s: silenciosas conexas por f y (na, nb) = B4 %s" % (tag, exp), conn == exp, conn)
        tot = {}
        for f, d in conn.items():
            for kk, v in d.items():
                tot[kk] = tot.get(kk, 0) + v
        for (sg, sx), c in scan["pares"].items():
            a1, b1 = nanb(sg, lts)
            a2, b2 = nanb(sx, lts)
            kk = "%d,%d" % (a1 + a2, b1 + b2)
            tot[kk] = tot.get(kk, 0) + c
        got = [scan["f_s"], scan["pares_compatibles"], scan["singletes_que_pasan_con_clase_trivial"],
               scan["singletes_que_no_pasan"]]
        expP = [P3["f_s"], P3["pares_compatibles"], P3["singletes_que_pasan_con_clase_trivial"],
                P3["singletes_que_no_pasan"]]
        check("N1 %s: f_s, pares compatibles y singletes = P3 de B4 %s" % (tag, expP), got == expP, got)
        check("N1 %s: silenciosas totales con <= 4 sitios por (na, nb) = P3 de B4" % tag,
              tot == P3["silenciosas_todas"], tot)
        tri = esym_counts(sc["cx"]["nclass"], KMAX)
        c0 = RB4[name]["celdas"]
        first = c0[next(iter(c0))]["por_k"]
        check("N1 %s: truncacion ingenua = llamadas_ingenua_convencional de B4" % tag,
              all(tri[str(k)] == first[str(k)]["llamadas_ingenua_convencional"] for k in range(1, KMAX + 1)), tri)
        same, nb_ok = 0, True
        for p, q in grid12:
            _guard("N1, control fenomenologico de %s" % tag, deadline)
            key = "%s|%s" % (p, q)
            od, B = mot.odds_fen(sc["cx"], p, q)
            terms = pb_terms(groups_fen(sc["cx"], p, q), KMAX)
            nb_ok = nb_ok and terms[0] == B
            vals = rival_values(scan, od, B, terms)
            for k in range(1, KMAX + 1):
                rk = c0[key]["por_k"][str(k)]
                same += int(vals[k]["eps_ing"] == b4.frac(rk["eps_ing"]["exacto"])
                            and vals[k]["L_cat"] == b4.frac(rk["L_ing"]["exacto"]))
        check("N1 %s: B = Pr[|F| = 0] en las 12 celdas" % tag, nb_ok)
        check("N1 %s: eps_ing y L_ing con k = 1..4 identicos (fracciones) a B4 en las 12 celdas (%d de 48)"
              % (tag, same), same == 12 * KMAX, same)
        out[tag] = {"identicas": same, "f_s": scan["f_s"], "pares_compatibles": scan["pares_compatibles"],
                    "silenciosas_por_f": scan["silenciosas_por_f"], "coste": scan["coste"]}
    return out


def s2_eps(refs, grid12, kmax=3):
    """S2_34_6x8: eps_ing con k <= 3 identico a D1-CE y a B4 (12 celdas). Devuelve las fracciones."""
    RC, RB4 = refs["CE"]["S2_34_6x8"], refs["B4"]["S2_34_6x8"]
    N, M = RC["N"], RC["M"]
    check("N1 S2_34_6x8: N y M de D1-CE = B4 (%d, %d)" % (N, M), (N, M) == (RB4["N"], RB4["M"]), verbose=False)
    out, same = {}, 0
    for p, q in grid12:
        key = "%s|%s" % (p, q)
        terms = pb_terms([(N, p), (M, q)], kmax)
        out[key] = {}
        for k in range(1, kmax + 1):
            e = eps_from_terms(terms, k)
            out[key][k] = e
            same += int(e == b4.frac(RC["celdas"][key]["por_k"][str(k)]["eps_ing"]["exacto"])
                        == b4.frac(RB4["celdas"][key]["por_k"][str(k)]["eps_ing"]["exacto"]))
    check("N1 S2_34_6x8: eps_ing con k = 1..%d identico a D1-CE y a B4 en las 12 celdas (%d de %d)"
          % (kmax, same, 12 * kmax), same == 12 * kmax, same)
    return out


def circuit_rival(pool, ce, rect, m, refs, cells, deadline):
    """Circuito: contexto, lectura del catalogo y eps_ing, L_cat por celda y k; controles con T2-0."""
    _guard("N1, rival de m = %d" % m, deadline)
    sc = strip_context(ce, rect, m, False, refs["B4"])
    tag, name = sc["tag"], sc["name"]
    lc = law_cost(sc["cx"], cells)
    scan = scan_catalogue(pool, sc, refs["MAN"], refs["T21M"], deadline)
    tri = esym_counts(sc["cx"]["nclass"], KMAX)
    RT = refs["T20"].get(name)
    if RT is not None:
        kT = T20_KMAX[m]
        got = [scan["silenciosas_por_f"][str(f)] for f in range(1, kT + 1)]
        exp = [RT["silenciosas_por_f"][str(f)] for f in range(1, kT + 1)]
        check("N1 %s: silenciosas por f (f <= %d) = T2-0 %s" % (tag, kT, exp), got == exp, got)
        check("N1 %s: truncacion ingenua (k <= %d) = T2-0" % (tag, kT),
              all(tri[str(k)] == RT["truncacion_ingenua_decodificaciones"][str(k)] for k in range(1, kT + 1)), tri)
        check("N1 %s: f_s = %d = T2-0" % (tag, F_S), scan["f_s"] == F_S == RT["f_s_circuito"], scan["f_s"])
    t = time.perf_counter()
    vals, lead_ok, b_ok = {}, 0, True
    for rc, rm in cells:
        _guard("N1, celdas de %s" % tag, deadline)
        key = "%s|%s" % (rc, rm)
        od, B = mot.odds_circ(sc["cx"], rc, rm)
        terms = pb_terms(groups_circ(sc["cx"], rc, rm), KMAX)
        b_ok = b_ok and terms[0] == B
        v = rival_values(scan, od, B, terms)
        vals[key] = {"B": B, "terminos_PB": terms, "por_k": v}
        if RT is not None:
            lead_ok += int(float(v[3]["L_cat"]) == RT["P_sil_dominante"][key]["P_sil_dominante"])
    check("N1 %s: B = Pr[|F| = 0] en las %d celdas" % (tag, len(cells)), b_ok)
    if RT is not None:
        check("N1 %s: L_ing(3) = termino dominante de T2-0 (float de la misma fraccion) en las %d celdas"
              % (tag, len(cells)), lead_ok == len(cells), lead_ok)
    scan["coste"]["segundos_cola_y_P3C_celdas"] = round(time.perf_counter() - t, 3)
    return {"contexto": sc, "leyes": lc, "lectura": scan, "valores": vals, "truncacion_ingenua": tri}


def scan_summary(scan):
    """Datos del catalogo que bastan para recalcular el rival sin releerlo (evidencia para --verificar)."""
    return {"silenciosas_por_firma": {str(s): v for s, v in sorted(scan["silenciosas"].items())},
            "pares_por_firmas": [[sg, sx, c] for (sg, sx), c in sorted(scan["pares"].items())],
            **{x: scan[x] for x in ("silenciosas_por_f", "f_s", "pares_compatibles",
                                    "singletes_que_pasan_con_clase_trivial", "singletes_que_no_pasan",
                                    "singletes_silenciosos", "entradas_de_singletes_benignos",
                                    "polimeros_silenciosos_de_3_sitios", "validas_por_f", "secciones", "bytes",
                                    "coste")},
            "firmas": len(scan["firmas"])}


def scan_from_summary(sm):
    return {"silenciosas": {int(s): v for s, v in sm["silenciosas_por_firma"].items()},
            "pares": {(sg, sx): c for sg, sx, c in sm["pares_por_firmas"]}}


# ------------------------------------------------------------------ R1: regresion del lector con datos sinteticos
def section_dict(a, polys):
    """Seccion con el formato de t2_motor: polys = [(YR, NB, [(firma, nA, nP, nM), ...])] por polimero."""
    YR, NB, off, sig, nA, nP, nM = [], [], [0], [], [], [], []
    for Y, N, ents in polys:
        YR.append(Y)
        NB.append(N)
        for s, x, p, q in ents:
            sig.append(s)
            nA.append(x)
            nP.append(p)
            nM.append(q)
        off.append(len(sig))
    return {"a": a, "YR": YR, "NB": NB, "off": off, "sig": sig, "nA": nA, "nP": nP, "nM": nM}


def write_group(path, sections):
    """Escribe un grupo como t2_motor (cabecera: ancla en 4 bytes y longitud en 8; seccion con pickle y zlib) y
    devuelve el SHA-256 del fichero, que read_group comprueba."""
    h = hashlib.sha256()
    with open(path, "wb") as fo:
        for d in sections:
            blob = zlib.compress(pickle.dumps(d, protocol=5), 1)
            head = d["a"].to_bytes(4, "little") + len(blob).to_bytes(8, "little")
            fo.write(head)
            fo.write(blob)
            h.update(head)
            h.update(blob)
    return h.hexdigest()


def regresion_r1(scan_fn=None, tmp=None):
    """R1 (revision de N1 v1): dos secciones de un grupo, anclas 0 y 1, con una silenciosa de 3 sitios en el mismo
    indice local i = 0 y NB distintos (7 y 26). Caso directo: un singlete benigno con Y = 8 solo es compatible con la
    primera (un par, no dos). Caso inverso: Y = 1 solo con la segunda (un par, no cero). Caso conjunto: los dos.
    Mascaras de prueba del lector, no historias fisicas. Se compara el lector real (scan_fn, por defecto _scan, con
    read_group del motor) sobre el grupo completo, sobre cada seccion en un grupo aparte y el conteo directo."""
    scan_fn = scan_fn or _scan
    g1, g2 = (1 << 12) | (2 << 6) | 3, (4 << 12) | (5 << 6) | 6     # firmas silenciosas de 3 simbolos
    x8, x1 = 9, 10                                                     # firmas de los singletes (1 simbolo)
    sil_a = (7, 7, [(g1, 1, 0, 1)])
    sil_b = (26, 26, [(g2, 1, 0, 1)])
    cases = {
        "directo": ([section_dict(0, [sil_a]), section_dict(1, [sil_b]), section_dict(3, [(8, 8, [(x8, 1, 1, 0)])])],
                    {(g1, x8): 1}),
        "inverso": ([section_dict(0, [sil_a, (1, 1, [(x1, 1, 1, 0)])]), section_dict(1, [sil_b])], {(g2, x1): 1}),
        "conjunto": ([section_dict(0, [sil_a, (1, 1, [(x1, 1, 1, 0)])]), section_dict(1, [sil_b]),
                      section_dict(3, [(8, 8, [(x8, 1, 1, 0)])])], {(g1, x8): 1, (g2, x1): 1}),
    }
    out = {}
    with (tempfile.TemporaryDirectory(prefix="t2_rival_R1_") if tmp is None else _Keep(tmp)) as d:
        for name, (secs, exp) in cases.items():
            p = os.path.join(d, "%s_grupo.bin" % name)
            r = scan_fn((0, p, write_group(p, secs), time.time() + 600))
            grouped = p3_pairs(r["single"], r["sil3"])[0]
            single, sil3 = [], []
            for u, s in enumerate(secs):
                q = os.path.join(d, "%s_seccion_%d.bin" % (name, s["a"]))
                rq = scan_fn((u, q, write_group(q, [s]), time.time() + 600))
                single += rq["single"]
                sil3 += rq["sil3"]
            apart = p3_pairs(single, sil3)[0]
            direct = {}
            for s in secs:
                for i in range(len(s["YR"])):
                    for j in range(s["off"][i], s["off"][i + 1]):
                        if not (C12 <= s["sig"][j] < C18 and s["nM"][j]):
                            continue
                        for t in secs:
                            for i2 in range(len(t["YR"])):
                                for j2 in range(t["off"][i2], t["off"][i2 + 1]):
                                    if t["sig"][j2] < C6 and t["nP"][j2] and t["YR"][i2] & s["NB"][i] == 0:
                                        kk = (s["sig"][j], t["sig"][j2])
                                        direct[kk] = direct.get(kk, 0) + s["nM"][j] * t["nP"][j2]
            out[name] = {"esperado": exp, "agrupado": grouped, "secciones_por_separado": apart, "directo": direct,
                         "polimeros_silenciosos_de_3_sitios": len(r["sil3"]), "secciones": r["secciones"]}
    return out


class _Keep:
    """Directorio dado (se conserva): misma interfaz que TemporaryDirectory."""

    def __init__(self, d):
        self.d = d

    def __enter__(self):
        os.makedirs(self.d, exist_ok=True)
        return self.d

    def __exit__(self, *a):
        return False


def regresion_r1_ok(rr):
    return all(v["agrupado"] == v["secciones_por_separado"] == v["directo"] == v["esperado"] for v in rr.values())


# ------------------------------------------------------------------ autoprueba (solo matematica)
def autoprueba():
    rng = random.Random(20261114)
    ok = []
    # 1. cola de Poisson-binomial agrupada frente a la enumeracion de todos los subconjuntos
    gs = [(3, F(1, 7)), (2, F(2, 9)), (1, F(0)), (2, F(1)), (1, F(5, 11))]
    flat = [pi for n, pi in gs for _ in range(n)]
    brute = [F(0)] * (len(flat) + 1)
    for mask in range(1 << len(flat)):
        pr = F(1)
        for i, pi in enumerate(flat):
            pr *= pi if mask >> i & 1 else 1 - pi
        brute[bin(mask).count("1")] += pr
    good = all(pb_terms(gs, K) == brute[:K + 1] for K in range(0, len(flat) + 1))
    good = good and sum(brute) == 1 and eps_from_terms(pb_terms(gs, 3), 3) == sum(brute[4:], F(0))
    ok.append(("cola exacta agrupada = enumeracion de los 2^9 subconjuntos (K = 0..9) y eps = cola", good))
    # 2. mismo resultado que naive de control_exacto (formula fenomenologica de B4) con p != q y q = 0
    good = True
    for n, m_, p, q in ((5, 3, F(1, 10), F(1, 20)), (4, 6, F(1, 1000), F(0)), (7, 2, F(1, 3), F(1, 3))):
        for k in range(0, 5):
            le = F(0)
            for a in range(k + 1):
                for b in range(k + 1 - a):
                    if q == 0 and b > 0:
                        continue
                    le += math.comb(n, a) * p ** a * (1 - p) ** (n - a) * math.comb(m_, b) * q ** b * (1 - q) ** (m_ - b)
            good = good and eps_from_terms(pb_terms([(n, p), (m_, q)], k), k) == 1 - le
    ok.append(("eps_ing = formula fenomenologica de B4 (control_exacto.naive), con q = 0 incluido", good))
    # 3. pares de P3C frente a la fuerza bruta en un catalogo sintetico
    NSY = 24
    adj = [0] * NSY
    for v in range(NSY):
        for w in range(NSY):
            if v != w and rng.random() < 0.15:
                adj[v] |= 1 << w
                adj[w] |= 1 << v

    def nbr(Y):
        nb = Y
        for v in bits(Y):
            nb |= adj[v]
        return nb

    single, sil3 = [], []
    for v in range(NSY):
        for c in range(rng.randint(1, 3)):
            Y = (1 << v) | (adj[v] & rng.getrandbits(NSY))
            s = (rng.randrange(6) * 8 + c) + 1
            nA = rng.randint(1, 3)
            nP = rng.randint(0, nA)
            single.append((Y, s, nP, nA, 0))
    for _ in range(40):
        vs = rng.sample(range(NSY), 3)
        Y = sum(1 << v for v in vs)
        ents = []
        for _ in range(rng.randint(1, 3)):
            syms = sorted(rng.randrange(40) for _ in range(3))
            sg = 0
            for sy in syms:
                sg = (sg << 6) | (sy + 1)
            ents.append((sg, rng.randint(1, 4)))
        sil3.append([Y, nbr(Y), ents])
    pairs, total, _, _ = p3_pairs(single, sil3)
    bf, bt = {}, 0
    for Y, NB, ents in sil3:
        for Yx, sx, nP, nA, nM in single:
            if nP and Yx & NB == 0:
                for sg, nm in ents:
                    bf[(sg, sx)] = bf.get((sg, sx), 0) + nm * nP
                    bt += nm * nP
    ok.append(("pares de P3C con indice por bit = fuerza bruta (%d pares)" % bt, pairs == bf and total == bt))
    # 4. valores del rival: (a) por numero de simbolos y (b) solo con k = 4
    od = {i: F(rng.randint(1, 9), 1000) for i in range(64)}
    sil = {}
    for f in (3, 4):
        for _ in range(5):
            syms = sorted(rng.randrange(40) for _ in range(f))
            s = 0
            for sy in syms:
                s = (s << 6) | (sy + 1)
            sil[s] = sil.get(s, 0) + rng.randint(1, 5)
    B = F(9, 10)
    terms = pb_terms([(10, F(1, 100))], KMAX)
    v = rival_values({"silenciosas": sil, "pares": pairs}, od, B, terms)

    def wt(s):
        r = F(1)
        while s:
            r *= od[(s & 63) - 1]
            s >>= 6
        return r

    a3 = sum((c * wt(s) for s, c in sil.items() if s < C18), F(0))
    a4 = sum((c * wt(s) for s, c in sil.items()), F(0))
    bb = sum((c * wt(sg) * wt(sx) for (sg, sx), c in pairs.items()), F(0))
    good = v[1]["L_cat"] == 0 and v[2]["L_cat"] == 0 and v[3]["L_cat"] == B * a3 and v[4]["L_cat"] == B * (a4 + bb)
    good = good and all(nsym(s) in (3, 4) for s in sil) and nsym(1) == 1 and nsym(C6) == 2 and nsym(C18 - 1) == 3
    ok.append(("L_cat por P3C: 0 con k <= 2, (a) con f <= k, (b) solo con k = 4; nsym de las firmas", good))
    # 5. H-comp: procedencia previa y regla de suma de T2-1
    good = (hcomp("circ_TIRA_6x6", 4, {})[0] == "COMPLETO" and hcomp("circ_TIRA_24x6", 3, {})[0] == "COMPLETO"
            and hcomp("circ_TIRA_24x6", 4, {})[0] == "SIN_REGLA"
            and hcomp("circ_TIRA_24x6", 4, {4: "COMPLETO"})[0] == "COMPLETO"
            and hcomp("circ_TIRA_48x6", 2, {4: "COMPLETO"})[0] == "COMPLETO"
            and hcomp("circ_TIRA_48x6", 4, {4: "NO_DECIDIBLE"})[0] == "NO_DECIDIBLE"
            and hcomp("circ_TIRA_48x6", 3, {4: "INCOMPLETO"})[0] == "SIN_REGLA"
            and hcomp("circ_TIRA_48x6", 4, {4: "INCOMPLETO"})[0] == "INCOMPLETO")
    ok.append(("H-comp: procedencia, COMPLETO cubre los cortes menores, INCOMPLETO y NO_DECIDIBLE no", good))
    # 6. truncacion ingenua de la literatura = combinaciones con clases
    nclass = [7, 7, 3, 1, 1]
    e = esym_counts(nclass, 3)
    bf2 = {str(k): sum(math.prod(nclass[i] for i in c) for f in range(1, k + 1) for c in combinations(range(5), f))
           for k in range(1, 4)}
    ok.append(("truncacion ingenua = suma de productos de clases por subconjunto", e == bf2))
    # 7. R1: lector real con dos secciones del mismo grupo y el mismo indice local
    rr = regresion_r1()
    ok.append(("R1: lector real (read_group y _scan) con dos secciones de un grupo, mismo indice local y NB "
               "distintos; casos directo, inverso y conjunto: agrupado = secciones por separado = conteo directo = "
               "esperado (%s)" % {k: sum(v["agrupado"].values()) for k, v in rr.items()}, regresion_r1_ok(rr)))
    # 8. R2: la espera de la lectura y los bucles de P3C y del rival consultan el techo (GUARD), sin procesos
    calls = []

    def g_stop(where):
        calls.append(where)
        if len(calls) >= 2:
            raise b4.Stop("techo sintetico (%s)" % where)

    class _Wait:
        def next(self, timeout=None):
            raise mp.TimeoutError

    class _NoPool:
        def imap_unordered(self, fn, tasks, chunksize=1):
            return _Wait()

    sc0 = {"tag": "sintetica", "NS": 1, "ngroups": 1, "gsize": 1}
    many3 = [(a << 12) | (b << 6) | c for a in (1, 2) for b in range(1, 64) for c in range(1, 64)][:GUARD_EVERY + 1]
    stops = []
    saved = GUARD[0]
    try:
        GUARD[0] = g_stop
        for fn in (lambda: scan_catalogue(_NoPool(), sc0, {"sintetica": {"0": {"archivo": "x", "sha256": "0"}}},
                                          None, time.time() + 600),
                   lambda: p3_pairs(single, sil3 * (GUARD_EVERY + 1)),
                   lambda: rival_values({"silenciosas": {s: 1 for s in many3}, "pares": {}},
                                        {i: F(1) for i in range(64)}, B, terms)):
            calls.clear()
            try:
                fn()
                stops.append(None)
            except b4.Stop as ex:
                stops.append(str(ex))
    finally:
        GUARD[0] = saved
    good = (stops[0] is not None and "en espera" in stops[0] and PROGRESS["sintetica"]["esperas"] == 2
            and PROGRESS["sintetica"]["grupos_leidos"] == 0 and all(s is not None for s in stops))
    try:
        scan_catalogue(_NoPool(), sc0, {"sintetica": {"0": {"archivo": "x", "sha256": "0"}}}, None, time.time() - 1)
        good = False
    except b4.Stop as ex:
        good = good and "techo de tiempo" in str(ex)
    ok.append(("R2: la espera de la lectura (sin resultados), los pares de P3C y los valores del rival consultan el "
               "techo y se detienen; el plazo vencido detiene la espera", good))
    for desc, g in ok:
        check("autoprueba N1: " + desc, g)
    return len(ok)


# ------------------------------------------------------------------ programa
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--autoprueba", action="store_true")
    ap.add_argument("--control-fen", action="store_true", help="control de B4 y de S2 (sin decoder ni DP)")
    ap.add_argument("--output")
    a = ap.parse_args()
    if a.autoprueba:
        n = autoprueba()
        print("AUTOPRUEBA N1: %d grupos de comprobaciones, todos PASS" % n, flush=True)
        raise SystemExit(0)
    if not a.control_fen or not a.output:
        raise SystemExit("uso: --autoprueba | --control-fen --output <nuevo.json>")
    if os.path.exists(a.output):
        raise SystemExit("--output ya existe")
    out = {"unidad": "T2-1 N1 (control fenomenologico, prueba de desarrollo)",
           "programa_sha256": sha_file(os.path.abspath(__file__)), "dependencias": SHA}
    code = 1
    deadline = T0_WALL + 3600.0
    try:
        ce, rect, _ = b4.load_all(with_b1=False)
        refs = load_refs()
        prov_checks(refs)
        grid12 = [(p, q) for p in ce.GRID_P for q in (p, p / 2, F(0))]
        with mp.get_context("spawn").Pool(PROCS) as pool:
            out["fen"] = fen_control(pool, ce, rect, refs, deadline, grid12)
        e = s2_eps(refs, grid12)
        out["S2_eps_ing"] = {key: {str(k): show(x) for k, x in d.items()} for key, d in e.items()}
        code = 0
    except (b4.CheckFail, b4.Stop) as ex:
        out["detenido_en"] = str(ex)
    except Exception as ex:  # nunca cuenta como PASS
        out["detenido_en"] = "error: %s: %s" % (type(ex).__name__, ex)
    out["llamadas_decoder"] = b4.CALLS[0]
    out["controles"] = b4.CHECKS
    out["segundos"] = round(time.perf_counter() - T0, 3)
    with open(a.output, "w", encoding="ascii", newline="\n") as fo:
        fo.write(json.dumps(out, ensure_ascii=True, indent=1) + "\n")
    print("N1 control fenomenologico:", "PASS" if code == 0 else "DETENIDO: " + out["detenido_en"], flush=True)
    raise SystemExit(code)


if __name__ == "__main__":
    main()
