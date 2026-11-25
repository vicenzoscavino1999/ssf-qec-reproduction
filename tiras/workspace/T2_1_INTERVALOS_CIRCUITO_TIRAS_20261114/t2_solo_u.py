#!/usr/bin/env python3
"""T2-1 (BORRADOR_PLAN_T2_1.txt, N2): adaptador SOLO_U para S2_34_6x8 (k_log = 4) y su enlace en TIRA_6x6.

Codigo preparado para revision (version 2, con las reparaciones R1 y R2 del dictamen DICTAMEN_CODIGO_T2_1.txt).
NO se ha ejecutado: se prueba solo con orden del usuario.

Contrato (N2 del borrador, con la precision del revisor del 14/11/2026):
- Catalogo por crecimiento (regla de B2, como build_catalogue/w_cat de T2-1M) con f <= 3, sin filtro fisico ni
  logico: no calcula clases, caracteres ni chis. Conserva la firma que da el peso: por soporte extendido Y y por
  firma (simbolos ordenados de 6 bits, tipo de ley * 8 + clase + 1), el numero n_all de decoraciones validas.
  En S2 hay dos tipos de ley, dato y check, asi que la firma es exactamente (a, b): n[Y,a,b].
- Se escribe en el formato en disco del motor (secciones por ancla con YR, NB, off, sig, nA; nP y nM vacios) y se
  evalua con t2_motor.dp_job SIN CAMBIOS, en su modo de un canal ("suma"), que solo lee nA: U con los odds
  fenomenologicos de cada celda (t2_motor.odds_fen), punto fijo de 384 bits, cotas inferior y superior. Cada
  componente cumple el corte por si mismo (firma con <= k simbolos); no hay corte global de la familia.
- Salidas por celda y k = 1..3: U en [U_lo, U_hi] y eps_comp = 1 - B U con extremos hacia fuera,
  B = (1-p)^N (1-q)^M. R_3 de referencia con eps_ing de RESULTADO_B4.json (el conductor lo recalcula con N1).
- S2: orden de la fase A (control_exacto.phase_a_order), como el control C- de B4.
- Enlace (control positivo, TIRA_6x6 fenomenologico, k_log = 1): mismo constructor con los tipos de ley del motor
  (t2_calibracion.site_classes restringido) y el orden del motor (t2_motor.layout). Las secciones con f <= 3
  deben coincidir entrada a entrada (Y, NB, firma, n_all) con trabajo/fen_TIRA_6x6 restringido a firmas con
  <= 3 simbolos, y los intervalos de eps_comp con k = 3 deben contener las fracciones de RESULTADO_B4.json.
  No se exige identidad byte a byte de extremos con el motor.
- No desactiva el guard k_log == 1 del motor (ese guard esta en catalogue_instance, que no se usa aqui).

Controles (detienen): SHA-256 de las dependencias; matrices de S2 = FASE_A.json y RESULTADO_B4.json; k_log de S2 =
4; validas por f = B4 y D1-CE; soportes distintos = B4 y D1-CE; polimeros, estados y terminos por k = gas_por_k
de B4 y D1-CE; fracciones exactas de eps_comp de B4 y de D1-CE dentro de los intervalos (12 celdas, k = 1..3);
anchura relativa de eps_comp <= 1e-30; ninguna ancla mala en la DP; el enlace de TIRA_6x6.

Presupuesto (R1): un contador global de llamadas al decoder para los dos catalogos, que se reserva antes de cada
llamada (techo estricto CALLS_MAX), y un plazo comun de reloj de pared desde el inicio del programa, consultado
antes de cada llamada, antes de lanzar la DP y, dentro de la DP sellada, tras cada ancla (su campo plazo recibe
el mismo plazo, sin minimo). La latencia de parada es una llamada al decoder o el proceso de un ancla en la DP. Si
se alcanza un techo, la unidad se detiene, sin resultados parciales, y guarda motivo y coste.

Evidencia (R2): por celda y corte se guardan U_lo y U_hi como enteros de punto fijo (U = valor / 2^PREC, con
PREC), B como fraccion exacta y eps_lo, eps_hi como fracciones exactas, ademas de la vista legible de 25 cifras;
tambien el manifiesto de los catalogos generados y los datos de cada trabajo de la DP. --verificar <salida.json>
recalcula desde esos datos, sin DP, los extremos, las contenciones de B4 y D1-CE y la anchura relativa.

Uso (desde D:\\paper ia cuantica\\qec_revision_matematica\\05_COMPROBACIONES):
  python -I -B T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/t2_solo_u.py --output <nuevo.json> --trabajo <dir nuevo>
  python -I -B T2_1_INTERVALOS_CIRCUITO_TIRAS_20261114/t2_solo_u.py --verificar <salida.json>
"""
import argparse
import hashlib
import json
import multiprocessing as mp
import os
import pickle
import sys
import time
import zlib
from array import array
from fractions import Fraction as F
from itertools import product

sys.set_int_max_str_digits(0)
T0 = time.perf_counter()
T0_WALL = time.time()
T21M ="T2_1M_MOTOR_ANCLAS_PUNTO_FIJO_20261113/"
B4D = "B4_TIRAS_ESCALADO_N_20261113/"
RCE = "D1_CONTROL_EXACTO_Y_METODOLOGIA_20261113/RESULTADO_CONTROL.json"
FASEA = "D1_VIABILIDAD_ALCANCE_CERTIFICADO_20261112/FASE_A.json"
SHA = {
    T21M + "t2_motor.py": "c0d297805a55ecfe514fd623c21ac5a1f272dd2b26a2f52a6f10a87205e6cf3a",
    T21M + "trabajo/MANIFIESTO.json": "4671bab0a5eec923b1fd7b9fdd7af499bd64dddf9efce4aedd372495d8f35730",
    B4D + "RESULTADO_B4.json": "640db376441656cfb72e16aa00e0e99e2fc2eb6047905dfe031f5d994634b5f7",
    RCE: "334d6485968c914c299bbc7e5362cf49f4b975b740d47d78bccfaec47255024d",
    FASEA: "e68435581f661be720cc03be80aae62261b8810fdc0aa09ea5b3af934b87de2c",
}


def sha_file(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


for _p, _s in SHA.items():
    if sha_file(_p) != _s:
        raise SystemExit("SHA-256 distinto del registrado: " + _p)
sys.path.insert(0, os.path.abspath(T21M))
import t2_motor as mot  # noqa: E402  (tras comprobar su SHA-256; el motor comprueba los suyos)

t2 = mot.t2
b4 = mot.b4
bits = b4.bits
check = b4.check

KMAX_S2 = 3
K_LINK = 3
PROCS = 14
SEG_MAX = 3600.0
CALLS_MAX = 5_000_000
RES_MAX = 300000
S2_VALID = [148, 10682, 526848]
S2_KLOG = 4
DEADLINE = T0_WALL + SEG_MAX
_COST = {"llamadas": 0}


def frac(v):
    return F(v["exacto"] if isinstance(v, dict) else v)


def check_deadline(where):
    if time.time() > DEADLINE:
        raise b4.Stop("techo de tiempo (%s)" % where)


def reserve_call(where):
    """Reserva una llamada al decoder en el contador global antes de hacerla: el total nunca pasa de CALLS_MAX."""
    check_deadline(where)
    if _COST["llamadas"] >= CALLS_MAX:
        raise b4.Stop("techo de llamadas (%s)" % where)
    _COST["llamadas"] += 1


# ------------------------------------------------------------------ catalogo sin filtro (n_all por Y y firma)
def build_catalogue_u(ce, inst, CL, LT, order, kmax, wdir, name):
    NS, adj = inst["NS"], inst["adjK"]
    dec, RES = ce.make_decoder(inst)
    pos = [0] * NS
    for i, v in enumerate(order):
        pos[v] = i
    check("%s: el orden es una permutacion de los %d sitios" % (name, NS), sorted(order) == list(range(NS)))
    adjR = [0] * NS
    for v in range(NS):
        adjR[pos[v]] = mot.relabel(adj[v], pos)
    gsize = -(-NS // mot.NGROUPS)
    ng = -(-NS // gsize)
    poly, CONN, info = {}, {}, {}
    rec_prev = None
    calls = 0
    for f in range(1, kmax + 1):
        t = time.perf_counter()
        if f == 1:
            cand = [1 << i for i in range(NS)]
        else:
            cs = set()
            for A, yu in rec_prev.items():
                nb = yu
                for v in bits(yu):
                    nb |= adj[v]
                nb &= ~A
                for x in bits(nb):
                    cs.add(A | (1 << x))
            cand = sorted(cs)
        rec_new = {}
        nval = nasg = 0
        for sm in cand:
            sites = list(bits(sm))
            yu = 0
            for idxs in product(*[range(len(CL[b])) for b in sites]):
                s = 0
                for b, i in zip(sites, idxs):
                    s ^= CL[b][i][1]
                if len(RES) > RES_MAX:
                    RES.clear()
                reserve_call("catalogo de " + name)
                _, touched = dec(s)
                calls += 1
                nasg += 1
                Y = sm | touched
                if not t2.connected(Y, adj, CONN):
                    continue
                nval += 1
                yu |= Y
                sig = 0
                for sy in sorted(LT[b] * 8 + i for b, i in zip(sites, idxs)):
                    sig = (sig << 6) | (sy + 1)
                slot = poly.get(Y)
                if slot is None:
                    poly[Y] = {sig: 1}
                else:
                    slot[sig] = slot.get(sig, 0) + 1
            if yu and f < kmax:
                rec_new[sm] = yu
        info[str(f)] ={"candidatos": len(cand), "asignaciones": nasg, "validas": nval,
                        "segundos": round(time.perf_counter() - t, 3)}
        print("  %s f=%d: %d candidatos, %d validas (%.1fs)" % (name, f, len(cand), nval, time.perf_counter() - t),
              flush=True)
        rec_prev = rec_new
    os.makedirs(wdir, exist_ok=False)
    byanc = {}
    for Y, slot in poly.items():
        YR = mot.relabel(Y, pos)
        byanc.setdefault(mot.anchor_of(YR), []).append((YR, slot))
    man, firmas = {}, set()
    npk = [0] * (kmax + 1)
    for g in range(ng):
        anchors = [a for a in sorted(byanc) if a // gsize == g]
        if not anchors:
            continue
        path = os.path.join(wdir, "grupo_%03d.bin" % g)
        hsh = hashlib.sha256()
        nbytes = npol = 0
        with open(path, "wb") as fo:
            for a in anchors:
                recs = sorted(byanc[a], key=lambda r: r[0])
                YRs, NBs = [], []
                off, sig, nA = array("Q", [0]), array("I"), array("I")
                for YR, slot in recs:
                    nb = YR
                    for i in bits(YR):
                        nb |= adjR[i]
                    YRs.append(YR)
                    NBs.append(nb)
                    sg = sorted(slot)
                    sig.extend(sg)
                    nA.extend(slot[s] for s in sg)
                    off.append(len(sig))
                    firmas.update(sg)
                    for k in range(1, kmax + 1):
                        if sg[0] < 1 << (6 * k):
                            npk[k] += 1
                blob = zlib.compress(pickle.dumps({"a": a, "YR": YRs, "NB": NBs, "off": off, "sig": sig, "nA": nA,
                                                   "nP": array("I"), "nM": array("I")}, protocol=5), 1)
                head = a.to_bytes(4, "little") + len(blob).to_bytes(8, "little")
                fo.write(head)
                fo.write(blob)
                hsh.update(head)
                hsh.update(blob)
                nbytes += len(head) + len(blob)
                npol += len(YRs)
        man[str(g)] = {"archivo": "grupo_%03d.bin" % g, "sha256": hsh.hexdigest(), "bytes": nbytes,
                       "polimeros": npol}
    return {"nombre": name, "dir": wdir, "NS": NS, "gsize": gsize, "ngroups": ng, "manifiesto": man,
            "firmas": sorted(firmas), "por_f": info, "polimeros_por_k": {str(k): npk[k] for k in range(1, kmax + 1)},
            "soportes": len(poly), "llamadas": calls, "solo_u": True}


# ------------------------------------------------------------------ DP de un canal del motor
def u_spec(cat, lab, k, odds):
    return {"etiqueta": lab, "dir": cat["dir"], "manifiesto": cat["manifiesto"], "NS": cat["NS"],
            "gsize": cat["gsize"], "ngroups": cat["ngroups"], "k": k, "modo": "suma",
            "celdas": [{"clave": lab, "odds": odds}], "exacto": False, "limite_grupos": None, "quitar": None,
            "firmas": cat["firmas"], "plazo": DEADLINE, "prioridad": 0.0}


def run_u_jobs(specs):
    for s in specs:
        if s["modo"] != "suma":
            raise ValueError("el catalogo SOLO_U solo admite el modo de un canal")
    check_deadline("antes de la DP")
    ctx = mp.get_context("spawn")
    R = {}
    with ctx.Pool(min(PROCS, len(specs))) as pool:
        for r in pool.imap_unordered(mot.dp_job, specs, chunksize=1):
            if "error" in r:
                raise b4.Stop("DP %s NO_MEDIDA: %s" % (r["etiqueta"], r["error"]))
            R[r["etiqueta"]] = r
            check_deadline("DP")
    return R


def eps_interval(B, vals):
    S = 1 << mot.PREC
    Ul, Uh = vals
    return 1 - B * F(Uh, S), 1 - B * F(Ul, S)


def evidence(B, vals, lo, hi):
    """R2: datos exactos para comprobar el intervalo sin repetir la DP."""
    return {"U_punto_fijo": [vals[0], vals[1]], "B": str(B), "eps_comp_exacto": [str(lo), str(hi)],
            "eps_comp": [mot.dec(lo, False), mot.dec(hi, True)]}


def job_data(R):
    return {lab: {x: v for x, v in r.items() if x != "valores"} for lab, r in sorted(R.items())}


def cat_data(cat):
    return {x: cat[x] for x in ("dir", "NS", "gsize", "ngroups", "manifiesto", "firmas", "por_f", "polimeros_por_k",
                                "soportes", "llamadas")}


# ------------------------------------------------------------------ S2
def run_s2(ce, rect, wroot, RB4, RC, grid12, out):
    name = "S2_34_6x8"
    print("== SOLO_U %s (k <= %d)" % (name, KMAX_S2), flush=True)
    inst = b4.make_inst(ce, rect, "s2", None)
    fa = json.load(open(FASEA, encoding="utf-8"))["instancias"][name]["matrices_sha256"]
    check("SOLO_U %s: matrices = FASE_A de D1 y = B4" % name,
          b4.mats_sha(inst) == fa == RB4[name]["matrices_sha256"])
    check("SOLO_U %s: k_log = %d (por eso solo U)" % (name, S2_KLOG), len(ce.logical_reps(inst)) == S2_KLOG)
    N, M, NS = inst["N"], inst["M"], inst["NS"]
    CL = [[(1 << v, inst["cols"][v])] for v in range(N)] + [[(0, 1 << c)] for c in range(M)]
    LTK = [("d",)] * N + [("c",)] * M
    lts = sorted(set(LTK))
    LT = [lts.index(x) for x in LTK]
    order = ce.phase_a_order(inst)
    cat = build_catalogue_u(ce, inst, CL, LT, order, KMAX_S2, os.path.join(wroot, "u_" + name), name)
    got = [cat["por_f"][str(f)]["validas"] for f in range(1, KMAX_S2 + 1)]
    rcf = RC[name]["catalogo"]["por_f"]
    check("SOLO_U %s: validas por f = B4 = D1-CE %s" % (name, S2_VALID),
          got == S2_VALID == [RB4[name]["catalogo"]["por_f"][str(f)]["validas"] for f in range(1, 4)]
          == [rcf[str(f)]["decoraciones_validas"] for f in range(1, 4)], got)
    check("SOLO_U %s: soportes distintos = B4 = D1-CE" % name,
          cat["soportes"] == RB4[name]["catalogo"]["polimeros_distintos"] == RC[name]["catalogo"]["polimeros_distintos"],
          cat["soportes"])
    cx = {"lts": lts, "N": N, "M": M}
    specs, Bs = [], {}
    for k in range(1, KMAX_S2 + 1):
        for p, q in grid12:
            key = "%s|%s" % (p, q)
            od, B = mot.odds_fen(cx, p, q)
            Bs[key] = B
            specs.append(u_spec(cat, "s2_k%d_%s" % (k, key), k, od))
    R = run_u_jobs(specs)
    bad = sum(r["anclas_malas"] for r in R.values())
    check("SOLO_U %s: C-anc en los %d trabajos (anclas malas = 0)" % (name, len(R)), bad == 0, bad)
    o = {"sitios": NS, "N": N, "M": M, "orden": "fase A", "precision_bits": mot.PREC, "catalogo": cat_data(cat),
         "trabajos_dp": job_data(R), "estructura": {}, "celdas": {}}
    for k in range(1, KMAX_S2 + 1):
        st = set()
        for p, q in grid12:
            r = R["s2_k%d_%s|%s" % (k, p, q)]
            st.add((r["polimeros"], r["estados"], r["terminos"]))
        ref = RB4[name]["gas_por_k"][str(k)]
        refc = RC[name]["gas_por_k"][str(k)]
        exp = (ref["polimeros"], ref["estados"], ref["terminos"])
        check("SOLO_U %s k=%d: polimeros, estados y terminos = B4 = D1-CE %s" % (name, k, list(exp)),
              st == {exp} and exp == (refc["polimeros"], refc["estados"], refc["terminos"]), sorted(st))
        o["estructura"][str(k)] = list(exp)
    worst = F(0)
    nin = 0
    for p, q in grid12:
        key = "%s|%s" % (p, q)
        cb, cr = RB4[name]["celdas"][key]["por_k"], RC[name]["celdas"][key]["por_k"]
        oc = {}
        for k in range(1, KMAX_S2 + 1):
            r = R["s2_k%d_%s" % (k, key)]
            lo, hi = eps_interval(Bs[key], r["valores"][0])
            xb, xr = frac(cb[str(k)]["eps_comp"]), frac(cr[str(k)]["eps_comp"])
            ok = lo <= xb <= hi and lo <= xr <= hi
            nin += int(ok)
            w = mot.relw(lo, hi)
            worst = max(worst, w)
            oc[str(k)] = dict(evidence(Bs[key], r["valores"][0], lo, hi), contiene_B4_y_D1CE=ok,
                              anchura_relativa=float(w))
            if k == 3:
                ei = frac(cb["3"]["eps_ing"])
                if lo <= 0:
                    raise b4.Stop("eps_comp(3) sin cota inferior positiva en " + key)
                R3 = (ei / hi, ei / lo)
                oc["R_3_referencia_eps_ing_B4"] = [mot.dec(R3[0], False), mot.dec(R3[1], True)]
                oc["R_3_referencia_exacto"] = [str(R3[0]), str(R3[1])]
                oc["R_3_contiene_B4"] = R3[0] <= frac(cb["3"]["R"]) <= R3[1]
        o["celdas"][key] = oc
    check("SOLO_U %s: los intervalos de eps_comp contienen las fracciones de B4 y D1-CE (%d de %d)"
          % (name, nin, 12 * KMAX_S2), nin == 12 * KMAX_S2, nin)
    check("SOLO_U %s: anchura relativa de eps_comp <= 1e-30 (%.3e)" % (name, float(worst)), worst <= mot.REL_W,
          float(worst))
    check("SOLO_U %s: R_3 de referencia contiene el R de B4 en las 12 celdas" % name,
          all(o["celdas"][c]["R_3_contiene_B4"] for c in o["celdas"]))
    o["anchura_relativa_max"] = float(worst)
    r3max = max(F(o["celdas"][c]["R_3_referencia_eps_ing_B4"][1]) for c in o["celdas"])
    o["D_t4_referencia"] = {"decision": "PASA" if r3max < 2 else "FALLA", "R_3_max_cota_superior": str(r3max),
                            "nota": "con eps_ing de B4; el conductor de T2-1 decide con eps_ing de N1"}
    o["segundos"] = round(time.perf_counter() - T0, 3)
    out["instancias"][name] = o


# ------------------------------------------------------------------ enlace en TIRA_6x6 (fenomenologico)
def load_sections(dirpath, man):
    secs = {}
    for g, ent in man.items():
        blob, idx = mot.read_group(os.path.join(dirpath, ent["archivo"]), ent["sha256"])
        mv = memoryview(blob)
        for a, (p, n) in idx.items():
            secs[a] = pickle.loads(zlib.decompress(mv[p:p + n]))
        mv.release()
    return secs


def restrict(sec, cap):
    out = {}
    for i, YR in enumerate(sec["YR"]):
        a, b = sec["off"][i], sec["off"][i + 1]
        ent = [(s, n) for s, n in zip(sec["sig"][a:b], sec["nA"][a:b]) if s < cap]
        if ent:
            out[YR] = (sec["NB"][i], ent)
    return out


def run_link(ce, rect, wroot, RB4, grid12, out):
    m = 6
    name = "TIRA_%dx6" % m
    tag = "fen_" + name
    print("== SOLO_U enlace %s (k <= %d)" % (tag, K_LINK), flush=True)
    inst = b4.make_inst(ce, rect, "tira", m)
    check("SOLO_U enlace %s: matrices = B4" % tag, b4.mats_sha(inst) == RB4[name]["matrices_sha256"])
    layer, L = t2.schedule(inst)
    locs = t2.locations(inst, layer, L)
    CL, _, LTK = t2.site_classes(inst, layer, L, locs, True)
    lts = sorted(set(LTK))
    LT = [lts.index(x) for x in LTK]
    check("SOLO_U enlace %s: simbolos caben en 6 bits" % tag, len(lts) <= 7 and max(len(c) for c in CL) <= 8)
    order = mot.layout(inst)[0]
    cat = build_catalogue_u(ce, inst, CL, LT, order, K_LINK, os.path.join(wroot, "u_" + tag), tag)
    man_m = json.load(open(T21M + "trabajo/MANIFIESTO.json", encoding="utf-8"))[tag]
    mine = load_sections(cat["dir"], cat["manifiesto"])
    ref = load_sections(os.path.join(T21M, "trabajo", tag), man_m)
    cap = 1 << (6 * K_LINK)
    diff = n_ent = 0
    for a in set(mine) | set(ref):
        x = restrict(mine[a], cap) if a in mine else {}
        y = restrict(ref[a], cap) if a in ref else {}
        n_ent += sum(len(v[1]) for v in x.values())
        if x != y:
            diff += 1
    check("SOLO_U enlace %s: secciones f <= %d = trabajo/%s entrada a entrada (Y, NB, firma, n_all; %d entradas)"
          % (tag, K_LINK, tag, n_ent), diff == 0, diff)
    cx = {"lts": lts, "N": inst["N"], "M": inst["M"]}
    specs, Bs = [], {}
    for p, q in grid12:
        key = "%s|%s" % (p, q)
        od, B = mot.odds_fen(cx, p, q)
        Bs[key] = B
        specs.append(u_spec(cat, "enlace_k%d_%s" % (K_LINK, key), K_LINK, od))
    R = run_u_jobs(specs)
    nin = 0
    cells = {}
    for p, q in grid12:
        key = "%s|%s" % (p, q)
        vals = R["enlace_k%d_%s" % (K_LINK, key)]["valores"][0]
        lo, hi = eps_interval(Bs[key], vals)
        x = frac(RB4[name]["celdas"][key]["por_k"][str(K_LINK)]["eps_comp"])
        ok = lo <= x <= hi
        nin += int(ok)
        cells[key] = dict(evidence(Bs[key], vals, lo, hi), contiene_B4=ok)
    check("SOLO_U enlace %s: eps_comp con k = %d contiene la fraccion de B4 en las 12 celdas (%d)"
          % (tag, K_LINK, nin), nin == 12, nin)
    out["enlace"] = {"instancia": tag, "k": K_LINK, "N": inst["N"], "M": inst["M"], "precision_bits": mot.PREC,
                     "entradas_comparadas": n_ent,
                     "celdas": cells, "catalogo": cat_data(cat), "trabajos_dp": job_data(R)}


def verify_output(path):
    """R2: recalcula sin DP, desde U_lo, U_hi, PREC y B guardados, los extremos exactos, la vista de 25 cifras,
    B = (1-p)^N (1-q)^M, las contenciones de B4 y D1-CE, la anchura relativa y R_3 de referencia."""
    d = json.load(open(path, encoding="ascii"))
    RB4 = json.load(open(B4D + "RESULTADO_B4.json", encoding="utf-8"))["instancias"]
    RC = json.load(open(RCE, encoding="utf-8"))["instancias"]
    bad, n = [], [0]

    def one(tag, key, c, NM, P, refs):
        S = 1 << P
        Ul, Uh = c["U_punto_fijo"]
        p, q = (F(s) for s in key.split("|"))
        B = F(c["B"])
        lo, hi = 1 - B * F(Uh, S), 1 - B * F(Ul, S)
        n[0] += 1
        if B != (1 - p) ** NM[0] * (1 - q) ** NM[1]:
            bad.append(tag + ": B")
        if [str(lo), str(hi)] != c["eps_comp_exacto"] or [mot.dec(lo, False), mot.dec(hi, True)] != c["eps_comp"]:
            bad.append(tag + ": extremos")
        if not all(lo <= x <= hi for x in refs):
            bad.append(tag + ": contencion")
        return lo, hi

    s2 = d["instancias"]["S2_34_6x8"]
    if s2["precision_bits"] != mot.PREC:
        bad.append("S2: precision")
    for key, oc in s2["celdas"].items():
        cb, cr = RB4["S2_34_6x8"]["celdas"][key]["por_k"], RC["S2_34_6x8"]["celdas"][key]["por_k"]
        for k in range(1, KMAX_S2 + 1):
            tag = "S2 %s k=%d" % (key, k)
            lo, hi = one(tag, key, oc[str(k)], (s2["N"], s2["M"]), s2["precision_bits"],
                         [frac(cb[str(k)]["eps_comp"]), frac(cr[str(k)]["eps_comp"])])
            if mot.relw(lo, hi) > mot.REL_W:
                bad.append(tag + ": anchura")
            if k == 3:
                ei = frac(cb["3"]["eps_ing"])
                if [str(ei / hi), str(ei / lo)] != oc["R_3_referencia_exacto"]:
                    bad.append(tag + ": R_3")
    ln = d["enlace"]
    for key, c in ln["celdas"].items():
        x = frac(RB4["TIRA_6x6"]["celdas"][key]["por_k"][str(ln["k"])]["eps_comp"])
        one("enlace %s" % key, key, c, (ln["N"], ln["M"]), ln["precision_bits"], [x])
    return {"comprobaciones": n[0], "fallos": bad}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output")
    ap.add_argument("--trabajo")
    ap.add_argument("--verificar", help="salida de una ejecucion anterior: comprobacion sin DP")
    a = ap.parse_args()
    if a.verificar:
        res = verify_output(a.verificar)
        print("VERIFICACION SIN DP: %d celdas y cortes; fallos: %s"
              % (res["comprobaciones"], res["fallos"] if res["fallos"] else "ninguno"), flush=True)
        raise SystemExit(1 if res["fallos"] or res["comprobaciones"] != 12 * (KMAX_S2 + 1) else 0)
    if not a.output or not a.trabajo:
        raise SystemExit("faltan --output y --trabajo")
    if os.path.exists(a.output) or os.path.exists(a.trabajo):
        raise SystemExit("--output o --trabajo ya existen")
    os.makedirs(a.trabajo)
    out = {"unidad": "T2-1 N2 SOLO_U", "programa_sha256": sha_file(os.path.abspath(__file__)), "dependencias": SHA,
           "instancias": {}}
    code = 1
    try:
        ce, rect, _ = b4.load_all(with_b1=True)
        RB4 = json.load(open(B4D + "RESULTADO_B4.json", encoding="utf-8"))["instancias"]
        RC = json.load(open(RCE, encoding="utf-8"))["instancias"]
        grid12 = [(p, q) for p in ce.GRID_P for q in (p, p / 2, F(0))]
        run_link(ce, rect, a.trabajo, RB4, grid12, out)
        run_s2(ce, rect, a.trabajo, RB4, RC, grid12, out)
        out["estado_N2"] = "CODIGO_ESCRITO; controles PASS en esta ejecucion"
        code = 0
    except (b4.CheckFail, b4.Stop) as ex:
        out["detenido_en"] = str(ex)
        out["estado_N2"] = "CODIGO_ESCRITO; detenido"
    except Exception as ex:  # se guarda lo que haya; nunca cuenta como PASS
        out["detenido_en"] = "error: %s: %s" % (type(ex).__name__, ex)
        out["estado_N2"] = "CODIGO_ESCRITO; detenido"
    out["coste"] = {"llamadas_decoder": _COST["llamadas"], "techo_llamadas": CALLS_MAX,
                    "segundos": round(time.time() - T0_WALL, 3), "techo_segundos": SEG_MAX}
    out["controles"] = b4.CHECKS
    out["segundos"] = round(time.perf_counter() - T0, 3)
    with open(a.output, "w", encoding="ascii", newline="\n") as fo:
        fo.write(json.dumps(out, ensure_ascii=True, indent=1) + "\n")
    print("SOLO_U:", "PASS" if code == 0 else "DETENIDO: " + out.get("detenido_en", ""), flush=True)
    raise SystemExit(code)


if __name__ == "__main__":
    main()
