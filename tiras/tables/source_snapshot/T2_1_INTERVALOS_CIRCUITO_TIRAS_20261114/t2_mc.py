#!/usr/bin/env python3
"""T2-1 (BORRADOR_PLAN_T2_1.txt, N4/U4/D-mc): Monte Carlo directo de control en TIRA_6x6 con ruido de circuito.

Codigo preparado para revision (version 2, con las reparaciones R1, R3 y R4 del dictamen DICTAMEN_CODIGO_T2_1.txt).
NO se ha ejecutado, ni siquiera en prueba: el Monte Carlo requiere el plan sellado y la orden del usuario.

Contrato (U4 del borrador):
- Celda m = 6, r_c = r_m = 1e-3; N = 7e9 historias fijadas antes, en bloques de 1e6 con un flujo
  pseudoaleatorio por bloque (MT19937 de Python, semilla derivada por SHA-256 de "T2-1 D-mc", la semilla maestra
  20261114 y el numero de bloque). Los resultados no dependen del reparto entre procesos.
- Por ubicacion se sortea la marginal X de su ley del apendice C (suma sobre las partes Z), calculada aqui
  por enumeracion de las leyes completas de Pauli: espera o preparacion, 2 r_c/3; CNOT, 4 r_c/15 por cada una
  de las tres partes X no nulas; lectura, w (1 - r_m) + (1 - w) r_m con w = 2 r_c/3. Las ubicaciones son
  independientes: la ley conjunta de las partes X es el producto. Las ubicaciones con fallo visible de cada
  historia se sortean por saltos geometricos sobre la secuencia de ubicaciones de cada grupo de igual
  probabilidad (con uniformes ideales, igual en ley a sortear cada ubicacion).
- Propagacion literal con t2_calibracion.propagate (T2-0) y decodificacion completa con el SSF de B4. Evento: el
  de B4 y T2-0 (H_X R = 0 y R fuera de S_X).
- Atajo: las historias con <= 2 ubicaciones con fallo visible no se decodifican y cuentan como no silenciosas
  (cada ubicacion afecta a una sola fibra y f_s = 3 en m = 6). Comprobacion del atajo, exhaustiva y separada del
  recuento: la historia vacia, todas las historias con una ubicacion con fallo visible (747) y todas las de dos
  ubicaciones distintas (278301), con todas sus categorias visibles: 279049 llamadas al decoder. Una silenciosa
  detiene la unidad. Comparte el decoder con los demas controles: comprueba la implementacion, no es una
  auditoria independiente de todas las dependencias.
- Recuento x de silenciosas sobre las N historias; intervalo de Clopper-Pearson bilateral al 99.9 % con t2_cp.py
  (evaluacion conservadora de las colas, alcance x <= t2_cp.XMAX). D-mc: COMPATIBLE si corta el intervalo
  certificado de la celda con k = 4, [L_lo, H_hi] con H_hi = min(1, L_hi + eps_hi); si no, INCOMPATIBLE. Solo
  significa que no se detecto incompatibilidad. Sin parada temprana ni repeticion.

Alcance estadistico (R3). El 99.9 % es la cobertura nominal del intervalo bajo el modelo binomial ideal: las N
historias independientes con la ley racional de las ubicaciones. El muestreo implementado se aparta de ese
modelo en estos supuestos de simulacion, que se declaran y no se acotan: MT19937 como fuente de uniformes
ideales, con flujos independientes entre bloques; resolucion de 53 bits de random(); inversion geometrica con
log y log1p en coma flotante; umbrales de categoria convertidos a coma flotante. No se reclama una cota del error
de ley. A4 (el certificado) es determinista y no depende de estos supuestos.

Presupuesto (R1). Un contador global compartido por todos los procesos cuenta las llamadas al decoder desde la
comprobacion del atajo: cada llamada se reserva en el contador antes de hacerse, de modo que el total nunca pasa
de CALLS_MAX (techo estricto). El plazo es comun (reloj de pared, desde el inicio del programa) y se consulta antes
de cada llamada al decoder, al empezar cada bloque y cada 2^16 pasos del muestreo. La parada es cooperativa: un
proceso termina su paso en curso y devuelve el bloque como interrumpido; el padre espera a todos, de modo que
todo el trabajo hecho queda contado. La latencia de parada es como mucho una llamada (propagacion y decodificacion)
o 2^16 pasos de muestreo, mas el vaciado de la cola; la salida registra la duracion maxima medida de una llamada
(ya en la comprobacion del atajo, antes del Monte Carlo) y el exceso sobre el plazo. Si se alcanza un techo: la
etapa en curso queda NO_EVALUADA, se guardan el motivo y el coste consumido y no se publica recuento parcial.

Trazabilidad (R4). El certificado entra como archivo JSON con sus metadatos (lo escribe el conductor N3 desde su
resultado) y se comprueba antes de muestrear: unidad T2-1, instancia TIRA_6x6, r_c, r_m, k = 4, SHA-256 del
decoder (b4_tiras.py) y de la propagacion (t2_calibracion.py), L_lo y H_hi como fracciones con
0 <= L_lo <= H_hi <= 1, y el SHA-256 del archivo de origen. La salida guarda el SHA-256 del certificado.
Con la muestra completa, el recuento (x, historias por numero de fallos e identificadores y recuentos de todos
los bloques) se escribe en <output>_recuento.json ANTES de evaluar Clopper-Pearson; un fallo de esa evaluacion
deja D-mc NO_EVALUADO y conserva el recuento.

Controles previos (detienen): SHA-256 de las dependencias; matrices de TIRA_6x6 iguales a B4; numero de
ubicaciones (527) y tamanos de los estratos (747 y 278301); marginales X iguales a las formas cerradas; leyes de
sitio reconstruidas desde las marginales iguales a t2_calibracion.law_of; propagate de cada fallo unico igual a
resp_formula. Control de contabilidad: contador global = llamadas del atajo + llamadas del Monte Carlo.

Uso (desde D:\\paper ia cuantica\\qec_revision_matematica\\05_COMPROBACIONES):
  python -I -B T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/t2_mc.py --output <nuevo.json> --solo-atajo
  python -I -B T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/t2_mc.py --output <nuevo.json> --certificado <cert.json>
"""
import argparse
import hashlib
import json
import math
import multiprocessing as mp
import os
import random
import sys
import time
from fractions import Fraction as F

sys.set_int_max_str_digits(0)
T0 = time.perf_counter()
T0_WALL = time.time()
UDIR = "T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/"
T20 = "T2_0_CALIBRACION_CIRCUITO_TIRAS_20261113/"
B4D = "B4_TIRAS_ESCALADO_N_20261113/"
SHA = {
    T20 + "t2_calibracion.py": "195bf3317f7af114f4cdc1daad714e828534ddd57abab29c81daa7126d17cc61",
    T20 + "RESULTADO_T2_0.json": "fb56c4324034549476813bd1a0d2bacd6f46c27d552fbf71ba2d0813d62df38e",
    B4D + "b4_tiras.py": "b1a2cb8f93affcc3713eadd4ff860a2d74756b91752011cd6b16da44def2dd15",
    B4D + "RESULTADO_B4.json": "640db376441656cfb72e16aa00e0e99e2fc2eb6047905dfe031f5d994634b5f7",
}


def sha_file(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


for _p, _s in SHA.items():
    if sha_file(_p) != _s:
        raise SystemExit("SHA-256 distinto del registrado: " + _p)
sys.path.insert(0, os.path.abspath(T20))
sys.path.insert(0, os.path.abspath(UDIR))
import t2_calibracion as t2  # noqa: E402  (tras comprobar su SHA-256)
import t2_cp  # noqa: E402  (su SHA-256 se registra en la salida; el plan sellado lo fijara)

b4 = t2.b4
bits = b4.bits
check = b4.check

M_STRIP = 6
K_CERT = 4
RC = F(1, 1000)
RM = F(1, 1000)
N_HIST = 7_000_000_000
BLOCK = 1_000_000
SEED = 20261114
PROCS = 14
CALLS_MAX = 200_000_000
SEG_MAX = 1.5 * 3600.0
SAMPLE_CHECK = 1 << 16
N_LOCS = 527
N_STRATUM1 = 747
N_STRATUM2 = 278301
N_SHORTCUT = 1 + N_STRATUM1 + N_STRATUM2
ALPHA = "0.001"
RES_MAX = 300000


# ------------------------------------------------------------------ presupuesto global (R1)
_BUD = {}


class Agotado(Exception):
    """Techo de llamadas o de tiempo: la etapa en curso queda NO_EVALUADA."""


def budget_init(calls, deadline):
    _BUD["calls"] = calls
    _BUD["deadline"] = deadline


def check_deadline():
    if time.time() > _BUD["deadline"]:
        raise Agotado("techo de tiempo")


def check_budget():
    """Sin reservar: plazo y llamadas restantes (al empezar un bloque)."""
    check_deadline()
    if _BUD["calls"].value >= CALLS_MAX:
        raise Agotado("techo de llamadas")


def reserve_call():
    """Reserva una llamada al decoder en el contador global antes de hacerla: el total nunca pasa de CALLS_MAX."""
    check_deadline()
    c = _BUD["calls"]
    with c.get_lock():
        if c.value >= CALLS_MAX:
            raise Agotado("techo de llamadas")
        c.value += 1


# ------------------------------------------------------------------ leyes por ubicacion (marginal X exacta)
PAULI_X = {"I": 0, "X": 1, "Y": 1, "Z": 0}


def marginal_1q(rc):
    law = {"I": 1 - rc, "X": rc / 3, "Y": rc / 3, "Z": rc / 3}
    out = {}
    for p, w in law.items():
        out[PAULI_X[p]] = out.get(PAULI_X[p], 0) + w
    return out


def marginal_cnot(rc):
    out = {}
    for p1 in "IXYZ":
        for p2 in "IXYZ":
            w = 1 - rc if p1 == p2 == "I" else rc / 15
            key = (PAULI_X[p1], PAULI_X[p2])
            out[key] = out.get(key, 0) + w
    return out


def marginal_read(rc, rm):
    one = marginal_1q(rc)
    out = {0: F(0), 1: F(0)}
    for x, w in one.items():
        for flip, wf in ((0, 1 - rm), (1, rm)):
            out[x ^ flip] += w * wf
    return out


def location_laws(locs, rc, rm):
    """Por ubicacion: lista de (categoria visible, probabilidad exacta). Categorias como en VISIBLE de T2-0."""
    m1, mc, mr = marginal_1q(rc), marginal_cnot(rc), marginal_read(rc, rm)
    laws = []
    for loc in locs:
        kind = loc[0]
        if kind == "dwait":
            laws.append([((1, 0), m1[1])])
        elif kind in ("prep", "await"):
            laws.append([((0, 1), m1[1])])
        elif kind == "read":
            laws.append([((0, 1), mr[1])])
        elif kind == "cnot":
            laws.append([(c, mc[c]) for c in ((1, 0), (0, 1), (1, 1))])
        else:
            raise ValueError("ubicacion desconocida " + kind)
    return laws


def site_law_from_locations(inst, layer, L, v_or_c, is_data, fd, fc, laws_by_loc):
    """Ley de un sitio reconstruida como XOR de las respuestas (resp_formula) de las ubicaciones de su fibra."""
    if is_data:
        v = v_or_c
        chk = sorted(bits(inst["cols"][v]), key=lambda c: layer[(c, v)])
        zero = (0,) * (len(chk) + 1)

        def key(e, s):
            return (1 if e else 0,) + tuple(s >> c & 1 for c in chk)
        fib = fd[v]
    else:
        c0 = v_or_c
        zero = (0,)

        def key(e, s):
            return (1 if s >> c0 & 1 else 0,)
        fib = fc[c0]
    dist = {zero: F(1)}
    for loc in fib:
        opts = [(zero, 1 - sum(w for _, w in laws_by_loc[loc]))]
        for cat, w in laws_by_loc[loc]:
            e, s = t2.resp_formula(inst, layer, loc, cat)
            opts.append((key(e, s), w))
        dist = t2.conv(dist, opts)
    return {k: w for k, w in dist.items() if w != 0}


# ------------------------------------------------------------------ procesos
_W = {}


def w_init():
    sys.set_int_max_str_digits(0)
    ce, rect, _ = b4.load_all(with_b1=False)
    inst = b4.make_inst(ce, rect, "tira", M_STRIP)
    dec, RES = ce.make_decoder(inst)
    layer, L = t2.schedule(inst)
    locs = t2.locations(inst, layer, L)
    laws = location_laws(locs, RC, RM)
    groups = {}
    for i, law in enumerate(laws):
        groups.setdefault(sum(w for _, w in law), []).append(i)
    G = []
    for q in sorted(groups):
        G.append((float(q), math.log1p(-float(q)), groups[q], q))
    cats = []
    for law in laws:
        tot = sum(w for _, w in law)
        acc, cum = F(0), []
        for c, w in law:
            acc += w / tot
            cum.append((float(acc), c))
        cum[-1] = (1.0, cum[-1][1])
        cats.append(cum)
    _W.update(ce=ce, inst=inst, dec=dec, RES=RES, layer=layer, L=L, locs=locs, groups=G, cats=cats, tmax=0.0)


def w_init_pool(calls, deadline):
    budget_init(calls, deadline)
    w_init()


def decode_event(faults):
    """faults: lista de (indice de ubicacion, categoria). Devuelve (silenciosa, pasa). Reserva antes la llamada."""
    reserve_call()
    t = time.perf_counter()
    W = _W
    ce, inst = W["ce"], W["inst"]
    E, s = t2.propagate(inst, W["layer"], W["L"], [(W["locs"][i], c) for i, c in faults])
    if len(W["RES"]) > RES_MAX:
        W["RES"].clear()
    cor, _ = W["dec"](s)
    r = E ^ cor
    passes = ce.lin(inst["cols"], r) == 0
    dt = time.perf_counter() - t
    if dt > W["tmax"]:
        W["tmax"] = dt
    return passes and not ce.inside(r, inst["bz"]), passes


def block_seed(b):
    h = hashlib.sha256(("T2-1 D-mc|%d|%d" % (SEED, b)).encode("ascii")).digest()
    return int.from_bytes(h[:8], "big")


def w_block(task):
    b, nh = task
    t = time.perf_counter()
    W = _W
    W["tmax"] = 0.0
    dec_n = sil = 0
    try:
        check_budget()
        rng = random.Random(block_seed(b))
        rnd = rng.random
        log = math.log
        hits = {}
        it = 0
        for _, lq, idxs, _ in W["groups"]:
            ng = len(idxs)
            T = nh * ng
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
                if not it % SAMPLE_CHECK:
                    check_deadline()
        hist = {}
        sil_by = {}
        cats = W["cats"]
        for h in sorted(hits):
            lst = hits[h]
            nf = len(lst)
            hist[nf] = hist.get(nf, 0) + 1
            if nf < 3:
                continue
            faults = []
            for i in sorted(lst):
                cum = cats[i]
                if len(cum) == 1:
                    faults.append((i, cum[0][1]))
                else:
                    u = rnd()
                    for th, c in cum:
                        if u < th:
                            faults.append((i, c))
                            break
            silent, _ = decode_event(faults)
            dec_n += 1
            if silent:
                sil += 1
                sil_by[nf] = sil_by.get(nf, 0) + 1
        hist[0] = nh - len(hits)
    except Agotado as ex:
        return {"b": b, "interrumpido": str(ex), "decodificadas": dec_n,
                "segundos": round(time.perf_counter() - t, 3), "decodificacion_max_s": W["tmax"], "pid": os.getpid()}
    return {"b": b, "historias": nh, "decodificadas": dec_n, "silenciosas": sil,
            "silenciosas_por_fallos": sil_by, "por_fallos": hist, "segundos": round(time.perf_counter() - t, 3),
            "decodificacion_max_s": W["tmax"], "pid": os.getpid()}


# ------------------------------------------------------------------ controles y comprobacion del atajo
def prechecks(out):
    w_init()
    W = _W
    ce, inst, layer, L, locs = W["ce"], W["inst"], W["layer"], W["L"], W["locs"]
    RB4 = json.load(open(B4D + "RESULTADO_B4.json", encoding="utf-8"))["instancias"]
    name = "TIRA_%dx6" % M_STRIP
    check("D-mc %s: matrices iguales a las de B4" % name, b4.mats_sha(inst) == RB4[name]["matrices_sha256"])
    check("D-mc: %d ubicaciones (plan: %d)" % (len(locs), N_LOCS), len(locs) == N_LOCS, len(locs))
    laws = location_laws(locs, RC, RM)
    n1 = sum(len(x) for x in laws)
    n2 = (n1 * n1 - sum(len(x) ** 2 for x in laws)) // 2
    check("D-mc: estratos 1 y 2 con %d y %d historias (plan: %d y %d)" % (n1, n2, N_STRATUM1, N_STRATUM2),
          (n1, n2) == (N_STRATUM1, N_STRATUM2), [n1, n2])
    w = 2 * RC / 3
    m1, mc, mr = marginal_1q(RC), marginal_cnot(RC), marginal_read(RC, RM)
    check("D-mc: marginales X = formas cerradas (2 r_c/3; 4 r_c/15 x 3; w (1 - r_m) + (1 - w) r_m)",
          m1[1] == w and all(mc[c] == 4 * RC / 15 for c in ((1, 0), (0, 1), (1, 1)))
          and mr[1] == w * (1 - RM) + (1 - w) * RM and sum(m1.values()) == 1 and sum(mc.values()) == 1
          and sum(mr.values()) == 1)
    laws_by_loc = {loc: law for loc, law in zip(locs, laws)}
    fd, fc = t2.fibres(inst, locs)
    _, _, LTK = t2.site_classes(inst, layer, L, locs, False)
    bad = 0
    for b in range(inst["NS"]):
        is_data = b < inst["N"]
        idx = b if is_data else b - inst["N"]
        mine = site_law_from_locations(inst, layer, L, idx, is_data, fd, fc, laws_by_loc)
        ref = {k: v for k, v in t2.law_of(LTK[b], L, RC, RM).items() if v != 0}
        bad += int(mine != ref)
    check("D-mc: leyes de sitio desde las marginales = law_of de T2-0 en los %d sitios" % inst["NS"], bad == 0, bad)
    bad = 0
    for i, law in enumerate(laws):
        for c, _ in law:
            bad += int(t2.propagate(inst, layer, L, [(locs[i], c)]) != t2.resp_formula(inst, layer, locs[i], c))
    check("D-mc: propagate de cada fallo unico = resp_formula (%d)" % n1, bad == 0, bad)
    out["parametros"] = {"m": M_STRIP, "r_c": str(RC), "r_m": str(RM), "N": N_HIST, "bloque": BLOCK,
                         "semilla_maestra": SEED, "alfa": ALPHA, "procesos": PROCS, "techo_llamadas": CALLS_MAX,
                         "techo_segundos": SEG_MAX, "pasos_de_muestreo_entre_comprobaciones": SAMPLE_CHECK,
                         "xmax_clopper_pearson": t2_cp.XMAX,
                         "ubicaciones_por_probabilidad_de_fallo_visible": {str(q): len(ix) for _, _, ix, q
                                                                           in W["groups"]}}
    return laws


def shortcut_check(laws, out):
    """Estratos 0, 1 y 2 enteros, con propagacion literal y decodificacion completa (279049 llamadas)."""
    t = time.perf_counter()
    _W["tmax"] = 0.0
    single = [(i, c) for i, law in enumerate(laws) for c, _ in law]
    n = sil = passed = 0
    try:
        s0, p0 = decode_event([])
        n, sil, passed = 1, int(s0), int(p0)
        for f in single:
            s, p = decode_event([f])
            n += 1
            sil += s
            passed += p
        for a in range(len(single)):
            ia = single[a][0]
            for bb in range(a + 1, len(single)):
                if single[bb][0] == ia:
                    continue
                s, p = decode_event([single[a], single[bb]])
                n += 1
                sil += s
                passed += p
    except Agotado as ex:
        out["comprobacion_atajo"] = {"decision": "NO_EVALUADA", "motivo": str(ex), "historias_hechas": n,
                                     "historias": N_SHORTCUT, "segundos": round(time.perf_counter() - t, 3)}
        raise
    out["comprobacion_atajo"] = {"historias": n, "silenciosas": sil, "pasan_filtro": passed,
                                 "segundos": round(time.perf_counter() - t, 3),
                                 "decodificacion_max_s": _W["tmax"]}
    check("D-mc atajo: %d historias con <= 2 ubicaciones con fallo visible, ninguna silenciosa" % n,
          n == N_SHORTCUT and sil == 0, {"historias": n, "silenciosas": sil})
    check("D-mc atajo: contador global = %d llamadas" % N_SHORTCUT, _BUD["calls"].value == N_SHORTCUT,
          _BUD["calls"].value)
    return n


def write_json(path, obj):
    with open(path, "w", encoding="ascii", newline="\n") as fo:
        fo.write(json.dumps(obj, ensure_ascii=True, indent=1) + "\n")
    return sha_file(path)


def run_mc(out, cert, calls, deadline, n_short, rec_path):
    nb, rem = divmod(N_HIST, BLOCK)
    assert rem == 0
    tasks = [(b, BLOCK) for b in range(nb)]
    ctx = mp.get_context("spawn")
    agg = {"historias": 0, "decodificadas": 0, "silenciosas": 0, "por_fallos": {}, "silenciosas_por_fallos": {}}
    blocks = {}
    cut = None
    dec_all = 0
    tmax = 0.0
    with ctx.Pool(PROCS, initializer=w_init_pool, initargs=(calls, deadline)) as pool:
        for r in pool.imap_unordered(w_block, tasks, chunksize=1):
            dec_all += r["decodificadas"]
            tmax = max(tmax, r["decodificacion_max_s"])
            if "interrumpido" in r:
                cut = cut or r["interrumpido"]
                continue
            if r["b"] in blocks:
                raise RuntimeError("bloque repetido %d" % r["b"])
            blocks[r["b"]] = [r["decodificadas"], r["silenciosas"], r["segundos"], r["pid"]]
            for k in ("historias", "decodificadas", "silenciosas"):
                agg[k] += r[k]
            for k, v in r["por_fallos"].items():
                agg["por_fallos"][k] = agg["por_fallos"].get(k, 0) + v
            for k, v in r["silenciosas_por_fallos"].items():
                agg["silenciosas_por_fallos"][k] = agg["silenciosas_por_fallos"].get(k, 0) + v
    used = calls.value
    elapsed = time.time() - T0_WALL
    out["coste"] = {"llamadas_decoder_total": used, "atajo": n_short, "monte_carlo": dec_all,
                    "techo_llamadas": CALLS_MAX, "segundos": round(elapsed, 3), "techo_segundos": SEG_MAX,
                    "exceso_sobre_el_plazo_s": round(max(0.0, elapsed - SEG_MAX), 3),
                    "decodificacion_max_s": tmax, "bloques_completos": len(blocks), "bloques": nb}
    check("D-mc: contador global de llamadas = atajo + Monte Carlo (%d)" % used, used == n_short + dec_all,
          [used, n_short, dec_all])
    finish_mc(out, cert, agg, blocks, nb, cut, rec_path)


def finish_mc(out, cert, agg, blocks, nb, cut, rec_path):
    """Tras recibir todos los bloques: muestra incompleta -> NO_EVALUADO sin recuento; completa -> recuento
    escrito en disco, despues Clopper-Pearson y D-mc."""
    if cut or len(blocks) != nb or agg["historias"] != N_HIST:
        out["D_mc"] = {"decision": "NO_EVALUADO", "motivo": cut or "bloques incompletos",
                       "lectura": "muestra incompleta: no se publica recuento ni intervalo"}
        return
    check("D-mc: suma de historias por numero de fallos = N", sum(agg["por_fallos"].values()) == N_HIST)
    check("D-mc: ninguna silenciosa con < 3 fallos visibles entre las decodificadas",
          all(int(k) >= 3 for k in agg["silenciosas_por_fallos"]))
    x = agg["silenciosas"]
    rec = {"x": x, "N": N_HIST, "decodificadas": agg["decodificadas"],
           "por_fallos": {str(k): v for k, v in sorted(agg["por_fallos"].items())},
           "silenciosas_por_fallos": {str(k): v for k, v in sorted(agg["silenciosas_por_fallos"].items())},
           "bloques_campos": ["bloque", "decodificadas", "silenciosas", "segundos", "pid"],
           "bloques": [[b] + blocks[b] for b in sorted(blocks)]}
    out["recuento"] = {x_: rec[x_] for x_ in ("x", "N", "decodificadas", "por_fallos", "silenciosas_por_fallos")}
    sha_rec = write_json(rec_path, {"unidad": "T2-1 D-mc recuento", "programa_sha256": out["programa_sha256"],
                                    "certificado_sha256": cert["sha256"], "recuento": rec})
    out["recuento_archivo"] = {"archivo": os.path.basename(rec_path), "sha256": sha_rec}
    try:
        cp = t2_cp.cp_interval(x, N_HIST, ALPHA)
        ok = t2_cp.intersects(cp, cert["lo"], cert["hi"])
    except Exception as ex:  # el recuento ya esta guardado; un fallo aqui nunca da COMPATIBLE
        out["D_mc"] = {"decision": "NO_EVALUADO",
                       "motivo": "fallo al evaluar Clopper-Pearson: %s: %s" % (type(ex).__name__, ex),
                       "lectura": "recuento completo conservado en %s" % os.path.basename(rec_path)}
        return
    out["clopper_pearson"] = cp
    out["D_mc"] = {"decision": "COMPATIBLE" if ok else "INCOMPATIBLE",
                   "lectura": ("no se detecto incompatibilidad" if ok else
                               "el intervalo de Monte Carlo no corta el certificado: se detiene la unidad"),
                   "cobertura": "nominal 99.9 % bajo el modelo binomial ideal y los supuestos de simulacion",
                   "certificado": [str(cert["lo"]), str(cert["hi"])]}


# ------------------------------------------------------------------ certificado (R4)
CERT_FIELDS = ("unidad", "instancia", "r_c", "r_m", "k", "b4_tiras_sha256", "t2_calibracion_sha256", "L_lo", "H_hi",
               "origen")


def load_cert(path):
    """Certificado de la celda escrito por N3 desde su resultado; se rechaza si no es el de esta celda."""
    with open(path, "rb") as fh:
        raw = fh.read()
    bad = []
    try:
        c = json.loads(raw.decode("ascii"))
        bad += [x for x in CERT_FIELDS if x not in c]
        if not bad:
            want = {"unidad": "T2-1", "instancia": "TIRA_%dx6" % M_STRIP, "k": K_CERT,
                    "b4_tiras_sha256": SHA[B4D + "b4_tiras.py"],
                    "t2_calibracion_sha256": SHA[T20 + "t2_calibracion.py"]}
            bad += [x for x, v in want.items() if c[x] != v]
            if F(c["r_c"]) != RC or F(c["r_m"]) != RM:
                bad.append("r_c o r_m")
            lo, hi = F(c["L_lo"]), F(c["H_hi"])
            if not 0 <= lo <= hi <= 1:
                bad.append("0 <= L_lo <= H_hi <= 1")
            org = c["origen"]
            if sha_file(org["archivo"]) != org["sha256"]:
                bad.append("SHA-256 del origen")
    except (ValueError, KeyError, TypeError, OSError, ZeroDivisionError) as ex:
        bad.append("%s: %s" % (type(ex).__name__, ex))
    if bad:
        raise SystemExit("certificado rechazado: " + "; ".join(bad))
    return {"lo": lo, "hi": hi, "archivo": path, "sha256": hashlib.sha256(raw).hexdigest(), "contenido": c}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--certificado", help="archivo JSON del certificado de la celda (lo escribe N3)")
    ap.add_argument("--solo-atajo", action="store_true")
    a = ap.parse_args()
    rec_path = os.path.splitext(a.output)[0] + "_recuento.json"
    if os.path.exists(a.output) or os.path.exists(rec_path):
        raise SystemExit("--output (o su _recuento.json) ya existe")
    cert = None
    if not a.solo_atajo:
        if not a.certificado:
            raise SystemExit("falta --certificado")
        cert = load_cert(a.certificado)
    out = {"unidad": "T2-1 D-mc", "programa_sha256": sha_file(os.path.abspath(__file__)),
           "t2_cp_sha256": sha_file(os.path.join(os.path.abspath(UDIR), "t2_cp.py")), "dependencias": SHA,
           "modo": "solo atajo" if a.solo_atajo else "Monte Carlo"}
    if cert:
        out["certificado"] = {"archivo": cert["archivo"], "sha256": cert["sha256"], "contenido": cert["contenido"]}
    calls = mp.get_context("spawn").Value("q", 0)
    deadline = T0_WALL + SEG_MAX
    budget_init(calls, deadline)
    code = 1
    try:
        laws = prechecks(out)
        n_short = shortcut_check(laws, out)
        if not a.solo_atajo:
            run_mc(out, cert, calls, deadline, n_short, rec_path)
        code = 0 if out.get("D_mc", {}).get("decision") in (None, "COMPATIBLE") else 1
    except Agotado as ex:
        out["detenido_en"] = "techo alcanzado: " + str(ex)
        out["coste"] = {"llamadas_decoder_total": calls.value, "techo_llamadas": CALLS_MAX,
                        "segundos": round(time.time() - T0_WALL, 3), "techo_segundos": SEG_MAX}
    except b4.CheckFail as ex:
        out["detenido_en"] = str(ex)
    except Exception as ex:  # se guarda lo que haya; nunca cuenta como COMPATIBLE
        out["detenido_en"] = "error: %s: %s" % (type(ex).__name__, ex)
    out["controles"] = b4.CHECKS
    out["segundos"] = round(time.perf_counter() - T0, 3)
    write_json(a.output, out)
    print("D-mc:", out.get("D_mc", {}).get("decision", "solo atajo" if a.solo_atajo and code == 0 else "detenido"),
          flush=True)
    raise SystemExit(code)


if __name__ == "__main__":
    main()
